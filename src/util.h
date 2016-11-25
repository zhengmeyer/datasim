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

#ifndef __UTIL_H__
#define __UTIL_H__

#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <mpi.h>
#include "configuration.h"
#include "subband.h"

/* parameters used for random number generation */
#define MEAN 1.0
#define STDEV 1.0
#define SEED 48573
#define MAXANT 20
#define MAXLEN 50

#define MASTER 0
#define COMMSIG 100
#define LOCK 200


typedef struct setup {
  int verbose;
  int test;                                 // test mode
  unsigned int seed;                        // random number generator seed
  float sfluxdensity;                       // source flux density in Jansky
  int antSEFDs[MAXANT] ;                    // antenna SEFD
  char inputfilename[MAXLEN];               // .input file name
  //vector<float> linesignal[2];              // line signal (sky)
  //vector<vector<float>> injectionsignal;    // injection signal (station)
} setup;

/*
 * Check whether the fractional part of a floating point number is 0
 */
bool is_integer(float num);
/*
 * Calculate the maximum spectrum resolution that can be used by all antennas
 */
int getSpecRes(Configuration* config, int configindex, float& specRes, size_t verbose);

/*
 * Calculate the maximum subband width among all antennas
 */
float getMaxChanFreq(Configuration* config, int configindex, size_t verbose);

/*
 * Calculate the number of samples to be generated per time block for the common signal
 */
int getNumSamps(Configuration* config, int configindex, float specRes, size_t verbose);

/*
 * Get the lowest start frequency among all antennas
 */
float getMinStartFreq(Configuration* config, int configindex, size_t verbose);

/*
 * Generate complex numbers using IPP Library
 * @cpDst Destination array where data will be store
 * @len   number of complex samples to generate
 * @stdev standard deviation
 * @rng_inst pointer to random number generator spec
 */
void gencplx(float* cpDst, size_t len, f32 stdev, gsl_rng *rng_inst, size_t verbose);

/*
 * loop through each subband of each antenna
 * move data from the second half of the array to the first half
 * set the process pointer to the proper location
 */
void movedata(Subband* subband, size_t verbose);

/*
 * loop through each subband of each antenna
 * select vdif packet size of data to process
 * process the data
 * quantization
 * pack to vdif
 */
int processAndPacketize(size_t framespersec, Subband* subband, Model* model, size_t verbose);

void genSignal(unsigned long stdur, int verbose, gsl_rng **rng_inst,
  float* commFreqSig, int numSamps, int lock, int myid, int numprocs, Subband* subband);
void copySignal(unsigned long stdur, int verbose, gsl_rng **rng_inst, int sfluxdensity, float* commFreqSig,
  int numSamps, int lock, int myid, int numprocs, Subband* subband);

#endif /* __UTIL_H__ */

/*
 * eof
 */
