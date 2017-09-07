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
 * Generate complex numbers using GSL
 * odd index is the real part, even index is the imaginary part
 */
void gencplx(float* cpDst, size_t len, f32 stdev, gsl_rng *rng_inst, size_t verbose)
{
  for(size_t idx = 0; idx < len; idx+=2)
  {
    cpDst[idx] = gsl_ran_gaussian_ziggurat(rng_inst, stdev);
    cpDst[idx+1] = gsl_ran_gaussian_ziggurat(rng_inst, stdev);
  }
}

/*
 * loop through each subband of each antenna
 * move data from the second half of the array to the first half
 * set the process pointer to the proper location
 */
void movedata(Subband* subband, size_t verbose)
{
  if(verbose >= 2) cout << "Move data in each subband array forward" << endl;
  subband->movedata();
}

/*
 * loop through each subband of each antenna
 * select vdif packet size of data to process
 * process the data
 * quantization
 * pack to vdif
 */
int processAndPacketize(size_t framespersec, Subband* subband, Model* model, size_t verbose)
{
  if(verbose >= 2)
    cout << "Antenna " << subband->getantIdx() << " Subband "
         << subband->getsbIdx() << " process vdif packet" << endl;
  subband->fillprocbuffer(); 
  subband->processdata(); 
  subband->updatevalues(model);

  subband->quantize();
  subband->writetovdif();
  if(verbose >= 2)
    cout << "current seconds is " << ((vdif_header *)subband->getvdifbuf())->seconds
         << ", frame number is " << ((vdif_header *)subband->getvdifbuf())->frame << endl;
  //update VDIF Header for the next packet
  nextVDIFHeader((vdif_header *) subband->getvdifbuf(), (int) framespersec);    

  return (EXIT_SUCCESS);
}

void prepareSignal(unsigned long stdur, int verbose, gsl_rng **rng_inst, int sfluxdensity,
  int numSamps, int myid, int numprocs, Subband* subband)
{
  float* commFreqSig1;                  // 0.5 seconds common frequency domain signal
  float* commFreqSig2;
  size_t workers = size_t(numprocs - 1);
  MPI_Request request1[workers], request2[workers];
  MPI_Status istatus1[workers], istatus2[workers];

  // allocate memory for the common frequency domain signal
  commFreqSig1 = new float [numSamps*2];
  commFreqSig2 = new float [numSamps*2];

  if(myid == MASTER)
  {
    gencplx(commFreqSig1, numSamps*2*stdur, STDEV, rng_inst[myid], verbose);
    for(size_t idx=0; idx < workers; idx++)
      MPI_Isend(&commFreqSig1[0], numSamps*2*stdur, MPI_FLOAT, idx, COMMSIG, MPI_COMM_WORLD, &request1[idx]);

  }
  else
  {
    MPI_Recv(&commFreqSig1[0], numSamps*2*stdur, MPI_FLOAT, MASTER, COMMSIG, MPI_COMM_WORLD, &istatus1[myid-1]);
    subband->fabricatedata(commFreqSig1, rng_inst[myid], sfluxdensity);
  }

  // free allocated common signal memory
  delete commFreqSig1;
  delete commFreqSig2;
}