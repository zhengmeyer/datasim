/*
 * datasim.cpp
 * The Non-zero baseline data simulator
 *
 * Author: Zheng Meyer-Zhao
 * 2013/10/24
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <gsl/gsl_randist.h>
#include "configuration.h"
#include "datasim.h"
#include "util.h"
#include "sbarr.h"
#include "vdifio.h"
#include "model.h"
#include "vdifzipper.h"

using namespace std;

static void usage(int argc, char **argv)
{
  cout << endl;
  cout << "Usage:  " << argv[0] << " [<options>] <input file>" << endl;
  cout << endl;
  cout << "  options can include:" << endl;
  cout << "     -h" << endl;
  cout << "     --help        display this information and quit." << endl;
  cout << endl;
  cout << "     -v" << endl;
  cout << "     --verbose     increase the verbosity of the output; -v -v for more." << endl;
  cout << endl;
  cout << "     -t" << endl;
  cout << "     --test        run in test mode, generate 1 second data for each station no matter what's given in the configuration file." << endl;
  cout << endl;
}   
      
     
int initVecSBArr(Configuration* config, int configindex, float specRes, 
                  float minStartFreq, vector<SBArr*> &sbVec, Model* model,
                  float tdur, int verbose)
{
  size_t length, startIdx, blksize;
  size_t vpsamps;   // number of samples in a vdif packet
  size_t vpbytes;   // number of bytes in a single-thread vdif packet
  size_t framebytes, numrecordedbands;
  size_t numdatastreams = (size_t)config->getNumDataStreams();
  float freq, bw;
  string antname;
  int mjd, seconds;
  f64* tempcoeffs;
  f64* delaycoeffs;
  double vptime;
  size_t framespersec;
  
  mjd = config->getStartMJD();
  seconds = config->getStartSeconds();
  if(verbose >= 1)
  {
    cout << "MJD is " << mjd << ", start seconds is " << seconds << endl;
  }

  // allocate memory for delaycoeffs
  tempcoeffs = vectorAlloc_f64(2);
  delaycoeffs = vectorAlloc_f64(2);

  for(size_t i = 0; i < numdatastreams; i++)
  {
    framebytes = (size_t)config->getFrameBytes(configindex, i);
    numrecordedbands = (size_t)config->getDNumRecordedBands(configindex, i);
    framespersec = config->getFramesPerSecond(configindex, i);
    vptime = 1.0 * 1e6 / framespersec;
    antname = config->getTelescopeName(i);
    // change the last character of the output vdif name to lower case for fourfit postprocessing
    antname.back() = tolower(antname.back());
    if(verbose >= 1)
    { 
      cout << "Antenna " << i << endl;
      cout << " framebytes is " << framebytes << endl;
      cout << " numrecordedbands is " << numrecordedbands << endl;
      cout << " antenna name is " << antname << endl;
    }

    // only consider scan 0 sourc 0
    // scanindex, offsettime in seconds, timespan in seconds, numincrements, antennaindex, scansourceindex, order, delaycoeffs
    model->calculateDelayInterpolator(0, 0, vptime*1e-6, 1, i, 0, 1, tempcoeffs);
    model->calculateDelayInterpolator(0, tempcoeffs[1]*1e-6, vptime*1e-6, 1, i, 0, 1, delaycoeffs);
    if(verbose >= 2)
    {
      cout << "delay in us for datastream " << i << " at offsettime 0s is " << tempcoeffs[1] << endl;
      cout << "delay in us for datastream " << i << " at offsettime " << tempcoeffs[1]*1e-6 << "s " << " is " << delaycoeffs[1] << endl;
    }

    // calculate vdif packet size in terms of bytes and number of samples
    // each sample uses 4 bits, as the sample is complex and we use 2 bits sampling
    // therefore 2 samples per byte
    vpbytes = (framebytes - VDIF_HEADER_BYTES) / numrecordedbands + VDIF_HEADER_BYTES;
    vpsamps = (framebytes - VDIF_HEADER_BYTES) / numrecordedbands * 2;
  
    for(size_t j = 0; j < numrecordedbands; j++)
    {
      freq = config->getDRecordedFreq(configindex, i, j);
      bw = config->getDRecordedBandwidth(configindex, i, j);
      if(!is_integer((freq - minStartFreq) / specRes))
      {
        cout << "StartIndex position is not an integer ... " << endl
             << "Something is wrong here ... " << endl;
        return EXIT_FAILURE;
      }
      else
        startIdx = (freq - minStartFreq) / specRes;
      blksize = bw / specRes; // number of samples to copy from startIdx
      length = bw * tdur * 2; // size of the array twice of tdur

      if(verbose >= 1)
      {
        cout << "Ant " << i << " subband " << j << ":" << endl
             << "  start index is " << startIdx << endl 
             << "  block size is " << blksize << " length is " << length << endl
             << "  each vdif packet has " << vpsamps << " samples" << endl
             << "  VDIF_HEADER_BYTES is " << VDIF_HEADER_BYTES << " vpbytes is " << vpbytes << endl
             << "  number of samples in vdif packet is " << vpsamps << " framebytes is " << framebytes << endl
             << "  recorded freq is " << freq << endl;
      }
      sbVec.push_back(new SBArr(startIdx, blksize, length, i, j, vpbytes, vpsamps, delaycoeffs, bw, antname, mjd, seconds, freq, verbose));
    }
  } 

  vectorFree(tempcoeffs);
  vectorFree(delaycoeffs);
  return EXIT_SUCCESS;
}

void freeVecSBArr(vector<SBArr*> &sbVec)
{
  vector<SBArr*>::iterator it;
  for(it = sbVec.begin(); it != sbVec.end(); it++)
  {
    (*it)->closevdif();
  }
  // free the memory allocated
  for(size_t i = 0; i < sbVec.size(); i++)
  {
    delete sbVec[i];
  }
  sbVec.clear();
}

int main(int argc, char* argv[])
{
  size_t verbose = 0;
  size_t test = 0;
  float tdur = 0.5 * 1e6;
  float testdur; // time duration of test mode
  string inputFile;

  if(argc < 2)
  {
    usage(argc, argv);
    return EXIT_FAILURE;
  }

  for(size_t i = 1; i < (size_t)argc; i++)
  {
    if(argv[i][0] == '-')
    {
      if(strcmp(argv[i], "-h") == 0 ||
         strcmp(argv[i], "--help") == 0) 
      {
        usage(argc, argv);
        return EXIT_SUCCESS;
      }
      else if(strcmp(argv[i], "-v") == 0 ||
        strcmp(argv[i], "--verbose") == 0)
      {
        ++verbose;
      }
      else if(strcmp(argv[i], "-t") == 0 ||
        strcmp(argv[i], "--test") == 0)
      {
        test = 1;
        testdur = 1;
        cout << "Entering test mode ..." << endl;
      }
    }
    else
    {
      if(!inputFile.empty())
      {
        cout << "inputFile value is " << inputFile << endl;
        cerr << "Multiple input files provided, only one is expected" << endl;
        cerr << "Run with -h or --help for help information" << endl;

        return EXIT_FAILURE;
      }
      inputFile = argv[i];   
    }
  }

  Configuration* config;
  Model* model;
  float dur;                            // duration of the simulated data in seconds
  float durus;                          // duration of the simulated data in microseconds
  float specRes, minStartFreq;
  int numdatastreams, numrecordedbands, freqindex;
  int numSamps;                         // number of samples to be generated per time block for the common signal
  int stime;                            // step time in microsecond == time block
  int configindex = 0;
  cf32* commFreqSig;                    // 0.5 seconds common frequency domain signal
  vector<SBArr*> sbVec;                 // vector of subband arrays

  double timer = 0.0, tt, vptime;       // vptime is the time duration of the samples within the packet
  float csigma = sqrt(CORR / (1-CORR));
  size_t framespersec;

  // initialize random number generator
  gsl_rng *rng_inst;
  gsl_rng_env_setup();
  rng_inst = gsl_rng_alloc(gsl_rng_ranlux389);
  gsl_rng_set(rng_inst, SEED);

  config = new Configuration(inputFile.c_str(), 0);
  model = config->getModel();

  // retrieve general information from config
  // if in test mode, only generate data for 1 second
  dur = test ? testdur : config->getExecuteSeconds();

  numdatastreams = config->getNumDataStreams();
  cout << "Generate " << dur << " seconds data for " << numdatastreams << " stations ..." << "\n"
       << "csigma is " << csigma << endl;

  for(int i = 0; i < numdatastreams; i++)
  {
    framespersec = config->getFramesPerSecond(configindex, i);
    vptime = 1.0 * 1e6 / framespersec;
    if(!is_integer(vptime))
    {
      cout << "VDIF packet time in microsecond is not an integer!! Something is wrong here ..." << endl;
      return EXIT_FAILURE;
    }
    numrecordedbands = config->getDNumRecordedBands(configindex, i);
    if(verbose >= 1)
    {
      cout << "Telescope " << config->getTelescopeName(i) << " has " << numrecordedbands << " recorded bands" << endl;
      cout << " Number of frames per second is " << framespersec << endl;
    }
    for(int j = 0; j < numrecordedbands; j++)
    {
      freqindex = config->getDRecordedFreqIndex(configindex, i, j); 
      if(verbose >= 1)
      {
        cout << "Subband " << j << ":" << "\n"
             << "  Frequency " << config->getFreqTableFreq(freqindex) << "\n"
             << "  Bandwidth " << config->getFreqTableBandwidth(freqindex) << endl; 
      }
    }
  }

  // general information from model
  if(verbose >= 1)
  {
    cout << "Number of scans is " << model->getNumScans() << endl;
    cout << "Scan start second is " << model->getScanStartSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
    cout << "Scan end second is " << model->getScanEndSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
    cout << "vptime in seconds is " << vptime/1e6 << endl;
  }

  // calculate specRes, number of samples per time block, step time
  if(getSpecRes(config, configindex, specRes, verbose) != EXIT_SUCCESS)
  {
    cout << "Failed to calculate spectral resolution ..." << endl;
    return EXIT_FAILURE;
  };
  numSamps = getNumSamps(config, configindex, specRes, verbose);
  stime = static_cast<int>(1 / specRes); // step time in microsecond
  if(verbose >= 1)
  {
    cout << "SpecRes is " << specRes << " MHz" << endl;
    cout << "number of samples per time block is " << numSamps << "\n"
         << "step time is " << stime << " us" << endl;
    cout << "tdur is " << tdur << endl;
  }
  minStartFreq = getMinStartFreq(config, configindex, verbose);

  // create and initialize array for each subband of each antenna
  if(initVecSBArr(config, configindex, specRes, minStartFreq, sbVec, model, tdur, verbose)
      != EXIT_SUCCESS)
  {
    cout << "Failed to initialize Subband Array ..." << endl;
    return EXIT_FAILURE;
  };
  // allocate memory for the common frequency domain signal
  commFreqSig = vectorAlloc_cf32(numSamps);

  cout << "Start generating data ...\n"
          "The process may take a while, please be patient!" << endl;

  // generate tdur time (e.g. 500000 microseconds) common frequency domain signal
  // and copy them to the right subband of each antenna
  genSignal(tdur/stime, commFreqSig, sbVec, numSamps, rng_inst, tdur, csigma, verbose);

  // loop though the simulation time dur
  do
  {
    // move data in each array from the second half to the first half
    // and set the process pointer to the proper location
    // i.e. data is moved half array ahead, therefore process pointer 
    // is moved half array ahead
    movedata(sbVec, verbose);
    // fill in the second half of each array, and
    // reset the current pointer of each array after data is generated
    genSignal(tdur/stime, commFreqSig, sbVec, numSamps, rng_inst, tdur, csigma, verbose);

    // set tt to the lowest process pointer time
    tt = getMinProcPtrTime(sbVec, verbose);
    if(verbose >= 2) cout << "the lowest process pointer is at time " << tt << " us" << endl;

    if(tt >= tdur)
    {
      cout << "**the lowest process pointer time is larger that tdur!!!\n"
              "**Something is wrong here!!!" << endl;
      return (EXIT_FAILURE);
    }

    // dur is in second
    // durus is in microsecond
    durus = dur * 1e6;
    while((tt < tdur) && (timer < durus))
    {
      // process and packetize one vdif packet for each subband array
      if(processAndPacketize(csigma, framespersec, sbVec, model, verbose))
        return(EXIT_FAILURE);
      tt += vptime;
      timer += vptime; 
      if(verbose >= 2) cout << "tt is " << tt << ", timer is " << timer << endl;
    }

  }while(timer < durus);
  
  freeVecSBArr(sbVec);  
  if(verbose >= 2) cout << "free memory for commFreq" << endl;
  vectorFree(commFreqSig); 

  // combine VDIF files
  vdifzipper(config, configindex, durus, verbose);

  cout << "All data has been generated successfully, bye!" << endl;

  // free random number generator
  gsl_rng_free(rng_inst);

  return(EXIT_SUCCESS);
}
