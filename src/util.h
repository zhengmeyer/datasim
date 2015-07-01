/*
 * util.h
 * VLBI data simulator
 *
 * Author: Zheng Meyer-Zhao
 * 2013/10/24
 *
 * Utilities used by the data simulator
 */

#ifndef __UTIL_H__
#define __UTIL_H__

#include <cstdlib>
#include <gsl/gsl_randist.h>
#include "configuration.h"
#include "sbarr.h"

/* parameters used for random number generation */
#define MEAN 1.0
#define STDEV 1.0
#define SEED 48573

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
void gencplx(cf32* cpDst, size_t len, f32 stdev, gsl_rng *rng_inst, size_t verbose);
/*
 * Generate TDUR time signal for all subband of all antennas
 */
void genSignal(size_t stdur, cf32* commFreqSig, vector<SBArr*>& sbVec, int numSamps, gsl_rng *rng_inst, float tdur, float csigma, size_t verbose);
/*
 * loop through each subband of each antenna
 * move data from the second half of the array to the first half
 * set the process pointer to the proper location
 */
void movedata(vector<SBArr*>& sbVec, size_t verbose);
/*
 * loop through each subband of each antenna
 * select vdif packet size of data to process
 * process the data
 * quantization
 * pack to vdif
 */
int processAndPacketize(float csigma, size_t framespersec, vector<SBArr*>& sbVec, Model* model, size_t verbose);
/*
 * calculate the lowest process pointer among all subband arrays
 */
double getMinProcPtrTime(vector<SBArr*>& sbVec, size_t verbose);

#endif /* __UTIL_H__ */

/*
 * eof
 */