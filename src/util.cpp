/*
 * datasim.cpp
 * The Non-zero baseline data simulator
 *
 * Author: Zheng Meyer-Zhao
 * 2013/10/24
 * To compile:
 * mpic++ `pkg-config --cflags --libs mpifxcorr mark5access difxmessage` -o datasim datasim.cpp -lmpifxcorr
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <gsl/gsl_randist.h>
#include "vdifio.h"
#include "util.h"
#include "datasim.h"
#include "configuration.h"

using namespace std;

/*
 * Check whether the fractional part of a floating point number is 0
 */
bool is_integer(float num)
{
  return floor(num)==num;
}

/*
 * Calculate the maximum spectrum resolution that can be used by all antennas
 */
int getSpecRes(Configuration* config, int configindex, float& specRes, size_t verbose)
{
  size_t cnt = 0;
  float tempSpecRes;
  bool isGCD = true; // Greatest Common Divisor
  vector<float> startfreq, freqdiff, bandwidth;
  vector<float>::iterator it, ity;
  for(int i = 0; i < config->getNumDataStreams(); i++)
  {
    startfreq.push_back(config->getDRecordedFreq(configindex, i, 0));
    bandwidth.push_back(config->getDRecordedBandwidth(configindex, i, 0));
  }
  for(it = startfreq.begin(); it != startfreq.end(); it++)
  {
    for(ity = it + 1; ity != startfreq.end(); ity++)
    {
      freqdiff.push_back(fabs(*it - *ity));
      if(verbose >= 2) cout << "Freqdiff is " << freqdiff.back() << endl;
    }
  }

  // loop through all frequency difference and bandwidth
  // we want to find the greatest common divisor
  // that all results after division are integers
  do
  {
    cnt++; // the highest chosen specRes is 0.5MHz, therefore, cnt starts from 1
    tempSpecRes = 1.0 / pow(2, cnt);
    for(it = freqdiff.begin(); it != freqdiff.end(); it++)
      isGCD &= is_integer(*it / tempSpecRes);
    for(it = bandwidth.begin(); it != bandwidth.end(); it++)
      isGCD &= is_integer(*it / tempSpecRes);

  // if one or more results are not integers
  // divide the value of tempSpecRes by 2
  // if GCD cannot be found at 1/(2^10), stop the calculation 
  }while(!isGCD && (cnt < 10));

  if(cnt == 10)
  {
    cout << "Failed to calculate spectral resoluation ..." << endl;
    return EXIT_FAILURE;
  }
  else
  {
    specRes = tempSpecRes;
    return EXIT_SUCCESS;
  }
}

/*
 * Calculate the maximum band covereage (bandwidth * nChans) among all antennas
 */
float getMaxChanFreq(Configuration* config, int configindex, size_t verbose)
{
  float chanFreq, maxChanFreq = 0;
  for(int i = 0; i < config->getNumDataStreams(); i++)
  {
    chanFreq = config->getDRecordedBandwidth(configindex, i, 0) * config->getDNumRecordedBands(configindex, i);
    maxChanFreq = (chanFreq > maxChanFreq) ? chanFreq : maxChanFreq;
  }
  if(verbose >= 1)
  {
    cout << "Max chan Freq is " << maxChanFreq << " MHz" << endl;
  }
  return maxChanFreq;
}

/*
 * Calculate the number of samples to be generated per time block for the common signal
 */
int getNumSamps(Configuration* config, int configindex, float specRes, size_t verbose)
{
  float maxChanFreq = getMaxChanFreq(config, configindex, verbose);
  assert(is_integer(maxChanFreq / specRes));
  return (int)maxChanFreq / specRes;
}

/*
 * Get the lowest start frequency among all antennas
 */
float getMinStartFreq(Configuration* config, int configindex, size_t verbose)
{
  float startFreq;
  float minStartFreq = config->getDRecordedFreq(configindex, 0, 0);
  for(int i = 0; i < config->getNumDataStreams(); i++)
  {
    startFreq = config->getDRecordedFreq(configindex, i, 0);
    minStartFreq = (startFreq < minStartFreq) ? startFreq : minStartFreq;
  }
  if(verbose >= 1)
  {
    cout << "Min start Freq is " << minStartFreq << " MHz" << endl;
  }
  return minStartFreq;
}

/*
 * Generate complex numbers using IPP Library
 * @cpDst Destination array where data will be store
 * @len   number of complex samples to generate
 * @mean  MEAN
 * @stdev standard deviation
 * @pSeed SEED
 */
void gencplx(cf32* cpDst, size_t len, f32 stdev, gsl_rng *rng_inst, size_t verbose)
{
  int status;
  f32* real = vectorAlloc_f32(len);
  f32* imag = vectorAlloc_f32(len);

  for(size_t idx = 0; idx < len; idx++)
  {
    real[idx] = gsl_ran_gaussian_ziggurat(rng_inst, stdev);
    imag[idx] = gsl_ran_gaussian_ziggurat(rng_inst, stdev);
  }
  
  status = vectorRealToComplex_f32(real, imag, cpDst, len);
  if(status != vecNoErr)
    cerr << "Error assembling random signals to complex!!!" << status << endl;

  if(verbose >= 3)
  {
    static size_t cnt = 1;
    if(cnt == 1)
    {
      for(size_t i = 0; i < len; i++)
        cout << "Real part value: " << real[i] << ", imaginary part value: " << imag[i] << endl;
      cnt--;
    }
  }
  for(size_t idx = 0; idx < len; idx++)
  {
    assert(cpDst[idx].re == real[idx]);
    assert(cpDst[idx].im == imag[idx]);
  }

  vectorFree(real);
  vectorFree(imag);
}

/*
 * Generate signal of time duration 'tdur' for all subband of all antennas
 */
void genSignal(size_t stdur, cf32* commFreqSig, vector<SBArr*>& sbVec, int numSamps, gsl_rng *rng_inst, float tdur, float csigma, size_t verbose)
{
  vector<SBArr*>::iterator it;
  if(verbose >= 1)
  {
    cout << "Generate " << tdur << " us signal" << endl;
  }
  for(size_t t = 0; t < stdur; t++)
  {
    gencplx(commFreqSig, numSamps, STDEV, rng_inst, verbose);
    // loop through each subband array of each antenna
    if(verbose >= 2) cout << "At time " << t << "us:" << endl;
    for(it = sbVec.begin(); it != sbVec.end(); ++it)
    {
      /* select data from commFreqSig 
       * add station noise
       * apply Ormsby filter
       * inverse DFT/FFT
       * copy data to the subband array
       */
      if(verbose >= 2)
      {
        cout << " Fabricate data for Antenna " << (*it)->getantIdx() << "subband " << (*it)->getsbIdx() << endl;
      }
      (*it)->fabricatedata(commFreqSig, rng_inst, csigma); 
    }
  }

  // after TDUR time signal is generated for each subband array
  // set the current pointer of each array back to the beginning of the second half
  for(it = sbVec.begin(); it != sbVec.end(); ++it)
  {
    (*it)->setcptr((*it)->getlength() / 2);
    if(verbose >= 2) 
    {
      cout << "Antenna " << (*it)->getantIdx() << " subband " << (*it)->getsbIdx()
                         << " set current pointer back to " << (*it)->getlength() / 2 << endl;
    }
  }
}

/*
 * loop through each subband of each antenna
 * move data from the second half of the array to the first half
 * set the process pointer to the proper location
 */
void movedata(vector<SBArr*>& sbVec, size_t verbose)
{
  if(verbose >= 2) cout << "Move data in each subband array forward" << endl;
  vector<SBArr*>::iterator it;
  for(it = sbVec.begin(); it != sbVec.end(); ++it)
  {
    (*it)->movedata();
  }
}

/*
 * loop through each subband of each antenna
 * select vdif packet size of data to process
 * process the data
 * quantization
 * pack to vdif
 */
int processAndPacketize(float csigma, size_t framespersec, vector<SBArr*>& sbVec, Model* model, size_t verbose)
{
  vector<SBArr*>::iterator it;
  for(it = sbVec.begin(); it != sbVec.end(); ++it)
  {
    if(verbose >= 2)
      cout << "Antenna " << (*it)->getantIdx() << " Subband "
           << (*it)->getsbIdx() << " process vdif packet" << endl;
    (*it)->fillprocbuffer(); 
    (*it)->processdata(); 
    (*it)->updatevalues(model);
 
    (*it)->quantize(csigma);
    (*it)->writetovdif();
    if(verbose >= 2)
      cout << "current seconds is " << ((vdif_header *)(*it)->getvdifbuf())->seconds
           << ", frame number is " << ((vdif_header *)(*it)->getvdifbuf())->frame << endl;
    //update VDIF Header for the next packet
    nextVDIFHeader((vdif_header *) (*it)->getvdifbuf(), (int) framespersec);    
  }
  return (EXIT_SUCCESS);
}

/*
 * calculate the lowest process pointer in terms of time among all subband arrays
 */
double getMinProcPtrTime(vector<SBArr*>& sbVec, size_t verbose)
{
  vector<SBArr*>::iterator it = sbVec.begin();
  double minprocptrtime = (*it)->getprocptr() * (1.0 / (*it)->getbandwidth());
  double temp;
  for(it = sbVec.begin(); it != sbVec.end(); ++it)
  {
    temp = (*it)->getprocptr() * (1.0 / (*it)->getbandwidth());
    if(verbose >= 2)
      cout << "Process Pointer in time for Antenna " << (*it)->getantIdx() << " Subband "
           << (*it)->getsbIdx() << " is " << temp << " us" <<endl;
    minprocptrtime = (minprocptrtime > temp) ? temp : minprocptrtime;
  }
  return minprocptrtime;
}
