/*****************************************************************************
*    <DataSim: VLBI data simulator>                                          *
*    Copyright (C) <2015> <Zheng Meyer-Zhao>                                 *
*                                                                            *
*    This file is part of DataSim.                                           *
                                                                             *
*    DataSim is free software: you can redistribute it and/or modify         *
*    it under the terms of the GNU General Public License as published by    *
*    the Free Software Foundation, either version 3 of the License, or       *
*    (at your option) any later version.                                     *
*                                                                            *
*    DataSim is distributed in the hope that it will be useful,              *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
*    GNU General Public License for more details.                            *
*                                                                            *
*    You should have received a copy of the GNU General Public License       *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include "vdifio.h"
#include "architecture.h"
#include "util.h"
#include "subband.h"
#include "datasim.h"
#include "model.h"
#include "vdifzipper.h"

using namespace std;

/*
 * Constructor
 */
Subband::Subband()
  : d_startIdx(0), d_blksize(0), d_length(0), d_antIdx(0), d_antSEFD(0), d_sbIdx(0),
    d_vpbytes(0), d_vpsamps(0), d_bandwidth(0), d_antname(0),
    d_mjd(0), d_seconds(0), d_freq(0), d_verbose(0), d_groupIdx(0)
{
  d_starttime = 0;

  // current pointer index for signal generation
  d_cptr = 0;

  // calculate the sample time and nearest sample
  d_sampletime = 0;
  d_vptime = 0;

  // use the time at the middle of a packet to calculate nearest sample
  d_nearestsample = 0;
  // fractinal sample error of the nearest sample (in us)
  // is nearestsampletime - starttime
  d_fracsamperror = 0;

  // initialize packet counter
  d_pkcounter = 0;
  // initialize shift counter
  d_shift = 0;

  // process pointer start at delay offset
  d_procptr =  0;

  // vdif file name to write data to
  d_filename = "";
  // open the output stream
  d_vdiffile.open("");

  try
  {
    //d_arr = vectorAlloc_cf32(d_length);
    d_arr = new cf32[d_length];
    // allocate memory for d_temp and d_tempt
    d_temp = vectorAlloc_cf32(d_blksize);
    d_tempt = vectorAlloc_cf32(d_blksize);
  }
  catch(bad_alloc& ba)
  {
    cout << "!!!!Failed allocating memory!!!!!" << endl;
    cout << "Exception caught: " << ba.what() << endl;
  }

  // initialize pDFTSpecC and bufsize
  vectorInitDFTC_cf32(&d_pDFTSpecCsig, d_blksize, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate,  &d_bufsigsize, &d_bufsig);
  vectorInitDFTC_cf32(&d_pDFTSpecCproc, d_vpsamps, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate,  &d_bufprocsize, &d_bufproc);
  vectorInitDFTC_cf32(&d_pDFTSpecCCToR, 2 * d_vpsamps, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate, &d_bufCToRsize, &d_bufCToR);


  // allocate memory for DFT buffer
  d_bufsig = vectorAlloc_u8(d_bufsigsize);
  d_bufproc = vectorAlloc_u8(d_bufprocsize);
  d_bufCToR = vectorAlloc_u8(d_bufCToRsize);

  // allocate memory for process buffer with size of vpsamps
  d_procbuffer = vectorAlloc_cf32(d_vpsamps);
  d_procbuffreq = vectorAlloc_cf32(d_vpsamps);
  d_procbufferrot = vectorAlloc_cf32(d_vpsamps);
  d_procbuffreqcorr = vectorAlloc_cf32(d_vpsamps);

  // allocate memory for complex and real signal array with size of 2*vpsamps
  // arrays used to convert complex to real
  d_buffreqtemp = vectorAlloc_cf32(2 * d_vpsamps);
  d_realC = vectorAlloc_cf32(2 * d_vpsamps);
  d_real = vectorAlloc_f32(2 * d_vpsamps);

  // allocate memory for d_vdifbuf
  d_vdifbuf = vectorAlloc_u8(d_vpbytes);  // initialize sample count and threshhold multiplier for quantization
  // allocate memory for d_delaycoeffs
  d_delaycoeffs = vectorAlloc_f64(2);
  // initialize delay coeffs
  d_delaycoeffs[0] = 0;
  d_delaycoeffs[1] = 0;
  // allocate memory for fractional sample error with size of vpsamps
  d_fracsamperrbuf = vectorAlloc_cf32(d_vpsamps);
  // allocate memory for fringe rotation buffer with size of vpsamps
  d_fringerotbuf = vectorAlloc_cf32(d_vpsamps);

  // initialize VDIF header
  createVDIFHeader((vdif_header *)d_vdifbuf, d_vpbytes - VDIF_HEADER_BYTES, d_sbIdx, BITS, 1, ISCOMPLEX, (char *)d_antname.c_str());
  setVDIFEpoch((vdif_header *)d_vdifbuf, d_mjd);
  setVDIFFrameMJDSec((vdif_header *)d_vdifbuf, d_mjd*86400 + d_seconds);

  d_sampcount = 0;
  d_tmul = 1.0;
  d_square = 0.0;
}

Subband::Subband(size_t const &startIdx, size_t const &blksize, size_t const &length, size_t const &antIdx, size_t const &antframespersec, unsigned int const &antSEFD,
             size_t const &sbIdx, size_t const &vpbytes, size_t const &vpsamps, f64* const &delaycoeffs, float const &bandwidth, string const &antname,
             int const &mjd, int const &seconds, float const &freq, size_t const &verbose, int groupIdx)
  : d_startIdx(startIdx), d_blksize(blksize), d_length(length), d_antIdx(antIdx), d_antframespersec(antframespersec), d_antSEFD(antSEFD),
    d_sbIdx(sbIdx), d_vpbytes(vpbytes), d_vpsamps(vpsamps), d_bandwidth(bandwidth), d_antname(antname),
    d_mjd(mjd), d_seconds(seconds), d_freq(freq), d_verbose(verbose), d_groupIdx(groupIdx)
{
  d_starttime = delaycoeffs[1];

  // current pointer index for signal generation
  d_cptr = d_length / 2;
  if(d_verbose >= 2) cout << "Start ptr index for copying data is " << d_cptr << endl;

  // calculate the sample time and nearest sample
  d_sampletime = 1.0 / static_cast<double>(d_bandwidth);
  d_vptime = d_vpsamps / static_cast<double>(d_bandwidth);

  // use the time at the middle of a packet to calculate nearest sample
  d_nearestsample = static_cast<int>((delaycoeffs[1] + delaycoeffs[0]*d_vptime*1e-6/2.0) * d_bandwidth + 0.5);
  // fractinal sample error of the nearest sample (in us)
  // is nearestsampletime - starttime
  double nearestsampletime = d_nearestsample * d_sampletime;
  d_fracsamperror = nearestsampletime - d_starttime;
  if(verbose >= 2)
  {
    cout << setprecision(20) << "Start time is " << d_starttime << endl;
    cout << "Nearest sample index is " << d_nearestsample << endl;
    cout << setprecision(20) << "bandwidth is " << d_bandwidth << ", nearestsample time is " << nearestsampletime << endl;
    cout << setprecision(20) << "Fractional sample error of the nearest sample is " << d_fracsamperror << endl;
  }

  // initialize packet counter
  d_pkcounter = 0;
  // initialize shift counter
  d_shift = 0;

  // process pointer start at delay offset
  d_procptr =  static_cast<int>(d_starttime * d_bandwidth + 0.5) + d_length / 2;
  if(d_verbose >= 2) cout << "Process pointer start at " << d_procptr << endl;

  // vdif file name to write data to
  stringstream ss, cc;
  ss << d_sbIdx;
  cc << d_groupIdx;
  d_filename = d_antname + "_" + ss.str() + "-" + cc.str() + ".vdif";

  // open the output stream
  d_vdiffile.open(d_filename.c_str());
  if(d_verbose >= 1) cout << "Open VDIF output file stream " << d_filename << endl;

  // allocate memory for d_arr
  if(d_verbose >= 1) cout << "Allocating memory for ant " << d_antIdx << " subband " << d_sbIdx << endl;
  try
  {
    //d_arr = vectorAlloc_cf32(d_length);
    d_arr = new cf32[d_length];
    // allocate memory for d_temp and d_tempt
    d_temp = vectorAlloc_cf32(d_blksize);
    d_tempt = vectorAlloc_cf32(d_blksize);
  }
  catch(bad_alloc& ba)
  {
    cout << "!!!!Failed allocating memory!!!!!" << endl;
    cout << "Exception caught: " << ba.what() << endl;
  }

  // initialize pDFTSpecC and bufsize
  vectorInitDFTC_cf32(&d_pDFTSpecCsig, d_blksize, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate,  &d_bufsigsize, &d_bufsig);
  vectorInitDFTC_cf32(&d_pDFTSpecCproc, d_vpsamps, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate,  &d_bufprocsize, &d_bufproc);
  vectorInitDFTC_cf32(&d_pDFTSpecCCToR, 2 * d_vpsamps, IPP_FFT_DIV_INV_BY_N, ippAlgHintAccurate, &d_bufCToRsize, &d_bufCToR);


  // allocate memory for DFT buffer
  d_bufsig = vectorAlloc_u8(d_bufsigsize);
  d_bufproc = vectorAlloc_u8(d_bufprocsize);
  d_bufCToR = vectorAlloc_u8(d_bufCToRsize);

  // allocate memory for process buffer with size of vpsamps
  d_procbuffer = vectorAlloc_cf32(d_vpsamps);
  d_procbuffreq = vectorAlloc_cf32(d_vpsamps);
  d_procbufferrot = vectorAlloc_cf32(d_vpsamps);
  d_procbuffreqcorr = vectorAlloc_cf32(d_vpsamps);

  // allocate memory for complex and real signal array with size of 2*vpsamps
  // arrays used to convert complex to real
  d_buffreqtemp = vectorAlloc_cf32(2 * d_vpsamps);
  d_realC = vectorAlloc_cf32(2 * d_vpsamps);
  d_real = vectorAlloc_f32(2 * d_vpsamps);

  // allocate memory for d_vdifbuf
  d_vdifbuf = vectorAlloc_u8(d_vpbytes);  // initialize sample count and threshhold multiplier for quantization
  // allocate memory for d_delaycoeffs
  d_delaycoeffs = vectorAlloc_f64(2);
  // initialize delay coeffs
  d_delaycoeffs[0] = delaycoeffs[0];
  d_delaycoeffs[1] = delaycoeffs[1];
  // allocate memory for fractional sample error with size of vpsamps
  d_fracsamperrbuf = vectorAlloc_cf32(d_vpsamps);
  // allocate memory for fringe rotation buffer with size of vpsamps
  d_fringerotbuf = vectorAlloc_cf32(d_vpsamps);

  // initialize VDIF header
  createVDIFHeader((vdif_header *)d_vdifbuf, d_vpbytes - VDIF_HEADER_BYTES, d_sbIdx, BITS, 1, ISCOMPLEX, (char *)d_antname.c_str());
  setVDIFEpochMJD((vdif_header *)d_vdifbuf, d_mjd);
  setVDIFFrameMJDSec((vdif_header *)d_vdifbuf, d_mjd*86400 + d_seconds);

  d_sampcount = 0;
  d_tmul = 1.0;
  d_square = 0.0;
}

/*
 * Destructor
 */
Subband::~Subband()
{
  //vectorFree(d_arr);
  delete [] d_arr;
  vectorFree(d_temp);
  vectorFree(d_tempt);
  vectorFree(d_bufsig);
  vectorFree(d_bufproc);
  vectorFree(d_bufCToR);
  vectorFree(d_procbuffer);
  vectorFree(d_procbuffreq);
  vectorFree(d_procbufferrot);
  vectorFree(d_procbuffreqcorr);
  vectorFree(d_buffreqtemp);
  vectorFree(d_realC);
  vectorFree(d_real);
  vectorFree(d_vdifbuf);
  vectorFree(d_delaycoeffs);
  vectorFree(d_fracsamperrbuf);
  vectorFree(d_fringerotbuf);
  if(d_verbose >= 1) cout << "Free memory for ant " << d_antIdx << " subband " << d_sbIdx << endl;
}

// copy constructor
Subband::Subband(Subband const &other)
{
  d_startIdx = other.d_startIdx;
  d_blksize = other.d_blksize;
  d_length = other.d_length;
  d_antIdx = other.d_antIdx;
  d_antframespersec = other.d_antframespersec;
  d_antSEFD = other.d_antSEFD;
  d_sbIdx = other.d_sbIdx;
  d_vpbytes = other.d_vpbytes;
  d_vpsamps = other.d_vpsamps;
  d_starttime = other.d_starttime;
  d_bandwidth = other.d_bandwidth;
  d_antname = other.d_antname;
  d_mjd = other.d_mjd;
  d_seconds = other.d_seconds;
  d_freq = other.d_freq;
  d_verbose = other.d_verbose;
  d_cptr = other.d_cptr;
  d_procptr = other.d_procptr;
  d_filename = other.d_filename;
  d_groupIdx = other.d_groupIdx;

  d_vdiffile.copyfmt(other.d_vdiffile);
  d_vdiffile.clear(other.d_vdiffile.rdstate());
  d_vdiffile.basic_ios<char>::rdbuf(other.d_vdiffile.rdbuf());

  d_sampletime = other.d_sampletime;
  d_vptime = other.d_vptime;
  d_nearestsample = other.d_nearestsample;
  d_fracsamperror = other.d_fracsamperror;
  d_pkcounter = other.d_pkcounter;
  d_shift = other.d_shift;
  d_vdifbuf = other.d_vdifbuf;
  d_delaycoeffs = other.d_delaycoeffs;
  d_fracsamperrbuf = other.d_fracsamperrbuf;
  d_fringerotbuf = other.d_fringerotbuf;
  d_arr = other.d_arr;
  d_temp = other.d_temp;
  d_tempt = other.d_tempt;
  d_procbuffer = other.d_procbuffer;
  d_procbuffreq = other.d_procbuffreq;
  d_procbufferrot = other.d_procbufferrot;
  d_procbuffreqcorr = other.d_procbuffreqcorr;
  d_buffreqtemp = other.d_buffreqtemp;
  d_realC = other.d_realC;
  d_real = other.d_real;
  d_bufsigsize = other.d_bufsigsize;
  d_bufprocsize = other.d_bufprocsize;
  d_bufCToRsize = other.d_bufCToRsize;
  d_bufsig = other.d_bufsig;
  d_bufproc = other.d_bufproc;
  d_bufCToR = other.d_bufCToR;
  d_pDFTSpecCsig = other.d_pDFTSpecCsig;
  d_pDFTSpecCproc = other.d_pDFTSpecCproc;
  d_pDFTSpecCCToR = other.d_pDFTSpecCCToR;
  d_sampcount = other.d_sampcount;
  d_tmul = other.d_tmul;
  d_square = other.d_square;
}

void Subband::swap(Subband& n2)
{
  using std::swap;
  swap(d_startIdx, n2.d_startIdx);
  swap(d_blksize, n2.d_blksize);
  swap(d_length, n2.d_length);
  swap(d_antIdx, n2.d_antIdx);
  swap(d_antframespersec, n2.d_antframespersec);
  swap(d_antSEFD, n2.d_antSEFD);
  swap(d_sbIdx, n2.d_sbIdx);
  swap(d_vpbytes, n2.d_vpbytes);
  swap(d_vpsamps, n2.d_vpsamps);
  swap(d_starttime, n2.d_starttime);
  swap(d_bandwidth, n2.d_bandwidth);
  swap(d_antname, n2.d_antname);
  swap(d_mjd, n2.d_mjd);
  swap(d_seconds, n2.d_seconds);
  swap(d_freq, n2.d_freq);
  swap(d_verbose, n2.d_verbose);
  swap(d_cptr, n2.d_cptr);
  swap(d_procptr, n2.d_procptr);
  swap(d_filename, n2.d_filename);

  std::ofstream d_tempvdif;
  d_tempvdif.copyfmt(d_vdiffile);
  d_tempvdif.clear(d_vdiffile.rdstate());
  d_tempvdif.basic_ios<char>::rdbuf(d_vdiffile.rdbuf());


  d_vdiffile.copyfmt(n2.d_vdiffile);
  d_vdiffile.clear(n2.d_vdiffile.rdstate());
  d_vdiffile.basic_ios<char>::rdbuf(n2.d_vdiffile.rdbuf());

  n2.d_vdiffile.copyfmt(d_tempvdif);
  n2.d_vdiffile.clear(d_tempvdif.rdstate());
  n2.d_vdiffile.basic_ios<char>::rdbuf(d_tempvdif.rdbuf());

  swap(d_sampletime, n2.d_sampletime);
  swap(d_vptime, n2.d_vptime);
  swap(d_nearestsample, n2.d_nearestsample);
  swap(d_fracsamperror, n2.d_fracsamperror);
  swap(d_pkcounter, n2.d_pkcounter);
  swap(d_shift, n2.d_shift);
  swap(d_vdifbuf, n2.d_vdifbuf);
  swap(d_delaycoeffs, n2.d_delaycoeffs);
  swap(d_fracsamperrbuf, n2.d_fracsamperrbuf);
  swap(d_fringerotbuf, n2.d_fringerotbuf);
  swap(d_arr, n2.d_arr);
  swap(d_temp, n2.d_temp);
  swap(d_tempt, n2.d_tempt);
  swap(d_procbuffer, n2.d_procbuffer);
  swap(d_procbuffreq, n2.d_procbuffreq);
  swap(d_procbufferrot, n2.d_procbufferrot);
  swap(d_procbuffreqcorr, n2.d_procbuffreqcorr);
  swap(d_buffreqtemp, n2.d_buffreqtemp);
  swap(d_realC, n2.d_realC);
  swap(d_real, n2.d_real);
  swap(d_bufsigsize, n2.d_bufsigsize);
  swap(d_bufprocsize, n2.d_bufprocsize);
  swap(d_bufCToRsize, n2.d_bufCToRsize);
  swap(d_bufsig, n2.d_bufsig);
  swap(d_bufproc, n2.d_bufproc);
  swap(d_bufCToR, n2.d_bufCToR);
  swap(d_pDFTSpecCsig, n2.d_pDFTSpecCsig);
  swap(d_pDFTSpecCproc, n2.d_pDFTSpecCproc);
  swap(d_pDFTSpecCCToR, n2.d_pDFTSpecCCToR);
  swap(d_sampcount, n2.d_sampcount);
  swap(d_tmul, n2.d_tmul);
  swap(d_square, n2.d_square);
}

// assignment operator
Subband const &Subband::operator=(Subband const &other)
{
  Subband tmp(other);
  swap(tmp);
  return *this;
}


/*
 * Public functions
 */
void Subband::fabricatedata(float* commFreqSig, gsl_rng *rng_inst, float sfluxdensity)
{
  // select data from commFreqSig
  // add station noise
  // apply Ormsby filter
  // inverse DFT/FFT
  // copy data to the subband array

  copyToTemp(commFreqSig);
  mulsfluxdensity(sfluxdensity);
  addstationnoise(rng_inst);
  normalizesignal(sfluxdensity);
  applyfilter();

  // value at DC should be zero
  assert(fabs(d_temp[0].re) < EPSILON);
  assert(fabs(d_temp[0].im) < EPSILON);

  inverseDFT(d_temp, d_tempt, d_pDFTSpecCsig, d_bufsig);
  copyToArr();
}

/*
 * Move data from the second half of the array to the first half
 * Reset the process pointer
 */
void Subband::movedata()
{
  for(size_t i = 0; i < d_length / 2; i++)
  {
    d_arr[i] = d_arr[i + d_length/2];
  }
  d_procptr -= d_length/2 ;
  if(d_verbose >= 2)
    cout << "Process pointer for Ant " << d_antIdx << " Subband " << d_sbIdx << " is at " << d_procptr << endl;
}

/*
 * Fill in the process buffer
 */
void Subband::fillprocbuffer()
{
  if(d_verbose >= 2)
    cout << "   fill process buffer for the current vdif packet ..." << endl;
  for(size_t i = 0; i < d_vpsamps; i++)
  {
    d_procbuffer[i] = d_arr[d_procptr+i];
  }
  d_procptr += d_vpsamps;
}

/*
 * Process data in the proc buffer
 * DFT
 * apply fractional sample correction
 * apply looffset correction
 * inverse DFT
 * apply fringe rotation
 * complex to real conversion
 */
void Subband::processdata()
{
  // apply fringe rotation
  //applyfringerotation();

  // apply DFT to proc data
  DFT(d_procbuffer, d_procbuffreq, d_pDFTSpecCproc, d_bufproc);

  // value at DC should be 0
  // manually set DC to 0, will this cause problem?
  d_procbuffreq[0].re = 0.0;
  d_procbuffreq[0].im = 0.0;

  // fill in fractional sample error buffer
  // and apply fractional sample error correction
  applyfracsamperrcorrection();

  // apply looffset correction

  if(d_verbose >= 2)
  {
    cout << "procbuf at DC is " << d_procbuffreq[0].re << " " << d_procbuffreq[0].im << endl;
    cout << "procbuf at freq 1 is " << d_procbuffreq[1].re << " " << d_procbuffreq[1].im << endl;
    cout << "procbuf at Nyquist is " << d_procbuffreq[d_vpsamps - 1].re << " " << d_procbuffreq[d_vpsamps - 1].im << endl;
  }

  // apply inverse DFT
  inverseDFT(d_procbuffreqcorr, d_procbuffer, d_pDFTSpecCproc, d_bufproc);

  // apply fringe rotation
  applyfringerotation();

  // complex to real conversion
  // 1. transform signal from time domain to frequency domain
  DFT(d_procbufferrot, d_procbuffreq, d_pDFTSpecCproc, d_bufproc);

  if(d_verbose >= 2)
    cout << "procbuf at DC after 2nd DFT is " << d_procbuffreq[0].re << " " << d_procbuffreq[0].im << endl;

  // 2. copy data from d_procbuffreq to the first half of d_buffreqtemp
  //    and fill in the second half of d_buffreqtemp with Hermitian property of the signal
  fillBuffreqtemp();

  // 3. transform signal from frequency domain back to time domain
  inverseDFT(d_buffreqtemp, d_realC, d_pDFTSpecCCToR, d_bufCToR);

  // 4. the imaginary part of values in d_realC should be 0
  //    copy the real part of the value to d_real
  complex_to_real();

  // update packet counter
  d_pkcounter++;
}

/*
 * apply phasecal
 */
void Subband::applyphasecal(int pcalinterval)
{
  for(size_t pcalfreq = 0; pcalfreq < (size_t)d_bandwidth/pcalinterval; pcalfreq++)
  {
    for(size_t sidx = 0; sidx < 2 * d_vpsamps; sidx++)
    {
      // sin(2*PI*index * pcalfreq / (2 * bandwitdh))
      // index = sampleindex mod (2 * bandwitdh)
      // since the signal is too strong, add a scaling factor of 1/500
      size_t modidx = sidx % (size_t)(2 * d_bandwidth);
      d_real[sidx] += 1.0/500.0 * sin(2*M_PI*modidx*pcalfreq*pcalinterval/(2*d_bandwidth));
    }
  }
  // assume the filter array has size of d_blksize
  // with values 0.0, 1/2, 4/5, 1, 1, ...., 1, 4/5, 1/2
  d_real[0] = 0.0;
  d_real[1] *= 1.0/2;
  d_real[2] *= 4.0/5;
  d_real[2*d_vpsamps-1] *= 1.0/2;
  d_real[2*d_vpsamps-2] *= 4.0/5;
}

void Subband::processdatawithpcal(int pcal)
{
  processdata();
  // 5. add phasecal
  applyphasecal(pcal);
}

/*
 * update pkcounter
 * update nearestsample
 * calculate fractional sample error for the next vdif packet and update procptr
 */
void Subband::updatevalues(Model* model)
{
  double prevdelay = d_delaycoeffs[1];
  double prevrate = d_delaycoeffs[0];
  // calculate delay interpolator for the next vdif packet
  // only consider scan 0 source 0
  // only consider first order
  // scanindex, offsettime, timespan in seconds, numincrements, antennaindex, scansourceindex, order, delaycoeffs
  model->calculateDelayInterpolator(0, (d_starttime+d_pkcounter*d_vptime + d_shift * d_sampletime)*1e-6, d_vptime*1e-6, 1, d_antIdx, 0, 1, d_delaycoeffs);
  if(d_verbose >= 2)
    cout << "delaycoeffs at offset " << (d_starttime+d_pkcounter*d_vptime + d_shift * d_sampletime)*1e-6 << " is " << d_delaycoeffs[0] << " " << d_delaycoeffs[1] << endl;

  // calculate fractional sample error for the next vdif packet
  d_fracsamperror += d_delaycoeffs[1] - (prevdelay + prevrate * d_vptime * 1e-6);

  // check if the fractional sample error is larger than 0.5*sampletime
  // or smaller than -0.5*sampletime
  // if so, move start position of the next vdif package accordingly

  if(d_verbose >= 2)
  {
    cout << "Antenna " << d_antIdx << " subband " << d_sbIdx << endl;
    cout << setprecision(20) << "  fractional sample error for packet " << d_pkcounter
         << " is " << d_fracsamperror << ", 0.5*sampletime is " << 0.5*d_sampletime << endl;
  }
  if(d_fracsamperror > 0.5*d_sampletime)
  {
    if(d_verbose >= 2)
      cout << "Antenna " << d_antIdx << " subband " << d_sbIdx << ":" << endl
           << "  fractional sample error is " << d_fracsamperror << ", larger than " << 0.5*d_sampletime << "!!" << endl
           << "  move process pointer to the next position!!" << endl;
    d_procptr++;
    d_shift++;
    d_fracsamperror -= d_sampletime;
    if(d_verbose >= 2)
      cout << "  process pointer now at position " << d_procptr << endl
           << "  new fractional sample error is " << d_fracsamperror << endl;
  }
  if(d_fracsamperror < -0.5*d_sampletime)
  {
    if(d_verbose >= 2)
      cout << "Antenna " << d_antIdx << " subband " << d_sbIdx << ":" << endl
           << "  fractional sample error is " << d_fracsamperror << ", smaller than " << -0.5*d_sampletime << "!!" << endl
           << "  move process pointer to the previous position!!" << endl;
    d_procptr--;
    d_shift--;
    d_fracsamperror += d_sampletime;
    if(d_verbose >= 2)
      cout << "  process pointer now at position " << d_procptr << endl
           << "  new fractional sample error is " << d_fracsamperror << endl;
  }
}

/*
 * Quantization
 */
void Subband::quantize()
{
  if(d_verbose >= 2)
    cout << "   Start quantization ..." << endl;
  size_t shift = 0;
  uint8_t bits, mask;
  uint8_t *optr = &d_vdifbuf[VDIF_HEADER_BYTES];
  float sample;
  //float thresh = sqrt(1 + csigma * csigma) * d_tmul;
  float thresh = d_tmul;

  for(size_t i = 0; i < 2 * d_vpsamps; i++)
  {
    sample = d_real[i];

    /* quantize it to 2 bits */
    if (sample > thresh) bits = 03;
    else if (sample < -thresh) bits = 00;
    else if (sample > 0.0) bits = 02;
    else if (sample < 0.0) bits = 01;
    else /* sample == 0.0 */ bits = 01; /* minor bias */

    /* update statistics up to limit */
    if (d_sampcount < BSMX) {
      d_sampcount++;
      d_square += sample*sample;
    }

    /* install the bits at proper location */
    bits <<= shift;
    mask = 0x3;
    mask <<= shift;
    (*optr) &= ~mask;
    (*optr) |= bits;

    /* move to the next bit location */
    shift += 2;
    while (shift > 7) { shift -= 8; optr++; } /* if/while */
  }
  d_tmul = sqrt(d_square / d_sampcount);
}

/*
 * Pack the processed data to VDIF packet
 * And update the vdif packet counter
 */
void Subband::writetovdif()
{
  if(d_verbose >= 2)
    cout << "   write to output vdif file ..." << endl;
  d_vdiffile.write((char *)d_vdifbuf, d_vpbytes);
}

/*
 * Close the output vdif stream
 */
void Subband::closevdif()
{
  d_vdiffile.close();
  if(d_verbose >= 3) cout << "Close output file stream " << d_filename << endl;
}

/*
 * Private functions
 */

/*
 * Copy data from common signal to temporary signal array
 */
void Subband::copyToTemp(float* commFreqSig)
{
  if(d_verbose >= 2)
    cout << "copy data to temp for ant " << d_antIdx << " subband " << d_sbIdx << endl;

  size_t index;
  for(size_t idx = 0; idx < d_blksize; idx++)
  {
    index = 2 * (d_startIdx+idx);
    d_temp[idx].re = commFreqSig[index];
    d_temp[idx].im = commFreqSig[index + 1];
    if((idx < 10) && (d_verbose >= 2))
    {
      cout << " commFreqSig[" << index << "] is " << commFreqSig[index]
           << " d_temp[" << idx << "] real is " << d_temp[idx].re<< endl;
    }
  }
}

/*
 * Change common signal amplitutde by multiplying square root of source flux density
 */
void Subband::mulsfluxdensity(float sfluxdensity)
{
  for(size_t idx = 0; idx < d_blksize; idx++)
  {
    d_temp[idx].re *= sqrt(sfluxdensity);
    d_temp[idx].im *= sqrt(sfluxdensity);
  }
}

/*
 * Add station noise
 */
void Subband::addstationnoise(gsl_rng *rng_inst)
{
  if(d_verbose >= 2)
    cout << "Adding station noise for ant " << d_antIdx << " subband " << d_sbIdx << endl;

  float* noise = new float [d_blksize*2];
  // generate station noise with the same variance as the common signal
  gencplx(noise, d_blksize*2, STDEV, rng_inst, d_verbose);

  // add station noise to d_temp
  for(size_t idx = 0; idx < d_blksize; idx++)
  {
    d_temp[idx].re += sqrt(d_antSEFD) * noise[2*idx];
    d_temp[idx].im += sqrt(d_antSEFD) * noise[2*idx+1];
  }

  delete noise;
}

/*
 * Signal normalization based on source flux density and antenna SEFD
 */
void Subband::normalizesignal(float sfluxdensity)
{
  for(size_t idx = 0; idx < d_blksize; idx++)
  {
    d_temp[idx].re /= sqrt(sfluxdensity + d_antSEFD);
    d_temp[idx].im /= sqrt(sfluxdensity + d_antSEFD);
  }
}

/*
 * Apply Ormsby filter to the temporary signal array
 */
void Subband::applyfilter()
{
  // assume the filter array has size of d_blksize
  // with values 0.0, 1/2, 4/5, 1, 1, ...., 1, 4/5, 1/2
  d_temp[0].re = 0.0;
  d_temp[1].re *= 1.0/2;
  d_temp[2].re *= 4.0/5;
  d_temp[d_blksize-1].re *= 1.0/2;
  d_temp[d_blksize-2].re *= 4.0/5;

  d_temp[0].im = 0.0;
  d_temp[1].im *= 1.0/2;
  d_temp[2].im *= 4.0/5;
  d_temp[d_blksize-1].im *= 1.0/2;
  d_temp[d_blksize-2].im *= 4.0/5;
  if(d_verbose >= 2)
  {
    cout << "Applied Ormsby filter to the temporary signal array for ant "
         << d_antIdx << " subband " << d_sbIdx  << endl;
  }
}

/*
 * Append the time domain signal in d_tempt to the signal array d_arr
 */
void Subband::copyToArr()
{
  if(d_verbose >= 2)
  {
    cout << "Copying time domain signal from temporary signal array d_tempt to signal array d_arr for ant "
         << d_antIdx << " subband " << d_sbIdx  << endl;
  }
  for(size_t i = 0; i < d_blksize; i++, d_cptr++)
  {
    d_arr[d_cptr] = d_tempt[i];

    if((i < 10) && (d_verbose >= 2))
    {
      cout << " current pointer index is " << d_cptr << " d_tempt[" << i << "] real is "
           << d_tempt[i].re << " d_arr[" << d_cptr << "] real is " << d_arr[d_cptr].re << endl;
    }
  }
}

/*
 * copy data from d_procbuffreq to the first half of d_buffreqtemp
 * fill in the second half of d_buffreqtemp with Hermitian property of the signal
 */
void Subband::fillBuffreqtemp()
{
  // copy data from d_procbuffreq to the first half of d_buffreqtemp
  for(size_t i = 0; i < d_vpsamps; i++)
  {
    // copy data
    d_buffreqtemp[i].re = d_procbuffreq[i].re;
    d_buffreqtemp[i].im = d_procbuffreq[i].im;

    // the second half of d_buffreqtemp possesses the Hermitian property of the signal
    d_buffreqtemp[2*d_vpsamps - i].re = d_buffreqtemp[i].re;
    d_buffreqtemp[2*d_vpsamps - i].im = -d_buffreqtemp[i].im;
  }
  // set the value at DC and Nyquist to 0
  // will this cause problem?
  d_buffreqtemp[d_vpsamps].re = d_buffreqtemp[0].re = 0.0;
  d_buffreqtemp[d_vpsamps].im = d_buffreqtemp[0].im = 0.0;
}

/*
 * copy the real part of d_realC to d_real
 */
void Subband::complex_to_real()
{
  for(size_t i = 0; i < 2 * d_vpsamps; i++)
  {
    // the imaginary part should be around zero
    if(i < 10 && d_verbose >= 2)
    {
      cout << "Absolute value of real part of d_realC at " << i << " is " << fabs(d_realC[i].re) << endl;
      cout << "Absolute value of imaginary part of d_realC at " << i << " is " << fabs(d_realC[i].im) << endl;
    }
    assert(fabs(d_realC[i].im) < EPSILON);
    d_real[i] = d_realC[i].re;
  }
}

/*
 * DFT
 */
void Subband::DFT(Ipp32fc* pSrc, Ipp32fc* pDst, vecDFTSpecC_cf32* pDFTSpecC, u8* buf)
{
  if(d_verbose >= 2) cout << "calling DFT ..." << endl;
  vectorDFT_CtoC_cf32(pSrc, pDst, pDFTSpecC, buf);
}

/*
 * inverse DFT
 */
void Subband::inverseDFT(Ipp32fc* pSrc, Ipp32fc* pDst, vecDFTSpecC_cf32* pDFTSpecC, u8* buf)
{
  if(d_verbose >= 2) cout << "calling inverse DFT ..." << endl;
  ippsDFTInv_CToC_32fc(pSrc, pDst, pDFTSpecC, buf);
}

/*
 * fill in fractional sample error buffer
 * and apply fractional sample correction
 */
void Subband::applyfracsamperrcorrection()
{

  int status;
  double arg = 2 * M_PI * d_bandwidth / (double)d_vpsamps;
  double fracerr;

  for(size_t idx = 0; idx < d_vpsamps; idx++)
  {
    // fractional sample error is calculated as 2*PI*fracerr/sampletime * idx/vpsamps
    fracerr = arg * idx * d_fracsamperror;
    d_fracsamperrbuf[idx].re = (float)cos(fracerr);
    d_fracsamperrbuf[idx].im = (float)sin(fracerr);
  }

  status = ippsMul_32fc(d_fracsamperrbuf, d_procbuffreq, d_procbuffreqcorr, d_vpsamps);
  if(status != vecNoErr)
    cout << "Error in application of fractional sample correction!!!" << endl;
}

/*
 * fill in fringe rotation buffer
 * and apply fringe rotation
 */
void Subband::applyfringerotation()
{
  int status;
  double phase;
  for(size_t idx = 0; idx < d_vpsamps; idx++)
  {
    // phase change is a function of time
    // use the same frequency for different subband
    // phase is 2*PI*freq*(delay+rate*idx/vpsamps)
    // rate is us/sec
    // delay is us
    // phase value should be opposite to what's in DiFX
    phase = 2 * M_PI * fraction_of(d_freq * (d_delaycoeffs[1] + d_delaycoeffs[0] * (double)idx / d_vpsamps));
    d_fringerotbuf[idx].re = cos(phase);
    d_fringerotbuf[idx].im = sin(phase);
  }

  status = ippsMul_32fc(d_fringerotbuf, d_procbuffer, d_procbufferrot, d_vpsamps);
  if(status != vecNoErr)
    cout << "Error in application of fringe rotation!!!" << endl;
}

/*
 * calculate the fractional part of a floating point number
 */
double Subband::fraction_of(double val)
{
  return val - rint(val - 0.5);
}
