/*****************************************************************************
*    <DataSim: VLBI data simulator>                                          * 
*    Copyright (C) <2015> <Zheng Meyer-Zhao>                                 *
*                                                                            *
*    This file is part of DataSim.                                           *
*                                                                            *
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
#include <unistd.h>
#include <vector>
#include <string>
#include <gsl/gsl_randist.h>
#include <getopt.h>
#include <mpi.h>
#include "configuration.h"
#include "datasim.h"
#include "util.h"
#include "subband.h"
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
  cout << "     -f" << endl;
  cout << "     --flux        source flux density in Jansky (default is 1 Jy)." << endl;
  cout << endl;
  cout << "     -s" << endl;
  cout << "     --sefd        antenna SEFDs in a comma-seperated list in Jansky.\n"
       << "                   If there are more antennas than provided SEFDs,\n"
       << "                   SEFD of the remaining antennas is set to 1000 Jansky."
       << "                   Maximum number of antennas supported is 20."<< endl;
  cout << endl;
  cout << "     -d" << endl;
  cout << "     --seed        random number generator seed." << endl;
  cout << endl;
  cout << "     -v" << endl;
  cout << "     --verbose     increase the verbosity of the output; -v -v for more." << endl;
  cout << endl;
  cout << "     -t" << endl;
  cout << "     --test        run in test mode, generate 1 second data for each station,\n"
       << "                   no matter what's given in the configuration file." << endl;
  cout << endl;
  /*
  cout << "     -l" << endl;
  cout << "     --line        line signal in the form of frequency,amplitude." << endl;
  cout << endl;
  cout << "     -i" << endl;
  cout << "     --injection   injection signal in the form of frequency,amplitude." << endl;
  cout << endl;
  */
}

static void cmdparser(int argc, char* argv[], setup &setupinfo)
{
  char tmp;
  static struct option long_options[] = {
    {"help",      no_argument,        0,  'h'},
    {"flux",      required_argument,  0,  'f'},
    {"sefd",      required_argument,  0,  's'},
    {"seed",      required_argument,  0,  'd'},
    {"verbose",   no_argument,        0,  'v'},
    {"test",      no_argument,        0,  't'},
    //{"line",      required_argument,  0,  'l'},
    //{"injection", required_argument,  0,  'i'},
    {0,           0,                  0,   0 }
  };
  int long_index = 0;
  while((tmp=getopt_long(argc,argv,"hf:s:d:vt",
              long_options, &long_index)) != -1)
  {
    switch(tmp)
    {
      case 'h':
        usage(argc, argv);
        exit (EXIT_SUCCESS);
        break;
      case 'f':
        if(*optarg == '-' || (atof(optarg) - 0.0) < EPSILON)
        {
          cerr << "Option -f requires a non-zero floating-point number as argument." << endl;
          exit (EXIT_FAILURE);
        }
        else
          setupinfo.sfluxdensity = atof(optarg);
        break;
      case 's':
        if(*optarg == '-' || *optarg == ' ')
        {
          cerr << "Option -s requires a comma-separated list of unsigned integer as argument." << endl;
          exit (EXIT_FAILURE);
        }
        else
        {
          istringstream ss(optarg);
          string token;
          int cnt = 0;
          while(getline(ss, token, ','))
          {
            setupinfo.antSEFDs[cnt] = atoi(token.c_str());
            cnt++;
          }
        }
        break;
      case 'd':
        if(*optarg == '-' || *optarg == ' ')
        {
          cerr << "Option -d requires an unsigned integer as argument." << endl;
          exit (EXIT_FAILURE);
        }
        else
        {
          setupinfo.seed = atoi(optarg);
        }
        break;
      case 'v':
        setupinfo.verbose++;
        break;
      case 't':
      setupinfo.test = 1;
      break;
      /*
      case 'l':
        if(*optarg == '-' || *optarg == ' ')
        {
          cerr << "Option -l requires a comma-separated list of float as argument." << endl;
          exit (EXIT_FAILURE);
        }
        else
        {
          istringstream ss(optarg);
          string token;
          for(size_t i = 0; i < setupinfo.linesignal.size(); i++)
          {
            getline(ss, token, ',');
            setupinfo.linesignal[i].push_back(atof(token.c_str()));
          }
        }
        break;
      case 'i':
        if(*optarg == '-' || *optarg == ' ')
        {
          cerr << "Option -i requires a comma-separated list of float as argument." << endl;
          exit (EXIT_FAILURE);
        }
        else
        {
          istringstream ss(optarg);
          string token;
          for(size_t i = 0; i < setupinfo.injectionsignal.size())
          {
            getline(ss, token, ',');
            setupinfo.injectionsignal[i].push_back(atof(token.c_str()));
          }
        }
        break;
      */
      default:
        usage(argc, argv);
        exit (EXIT_FAILURE);
        break;
    }
  }
  // check if input VDIF file name is given
  if(optind > argc - 1)
  {
    cerr << "No .input file provided, nothing to do ..." << endl;
    usage(argc, argv);
    exit (EXIT_FAILURE);
  }
  else if(optind < argc - 1)
  {
    cerr << "Multiple input files provided, only one is expected." << endl;
    usage(argc, argv);
    exit (EXIT_FAILURE);
  }
  else
  {
    strcpy(setupinfo.inputfilename, argv[optind]);
  }
}

int main(int argc, char* argv[])
{
  // initialize MPI
  MPI_Init(&argc, &argv);

  int numprocs, myid;
  double start, end, elapse;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  setup setupinfo;

  // master parse command line argument
  // and send struct setupinfo to all worker processes
  if(myid == MASTER)
  {
    setupinfo.verbose = 0;
    setupinfo.test = 0;
    setupinfo.seed = SEED;
    setupinfo.sfluxdensity = 10;     // source flux density in Jansky
    for(size_t idx = 0; idx < MAXANT; idx++)
      setupinfo.antSEFDs[idx] = 200;
    /*
    for(size_t i = 0; i < setupinfo.linesignal.size(); i++)
      setupinfo.linesignal[i] = 0;
    for(size_t i = 0; i < setupinfo.injectionsignal.size(); i++)
      setupinfo.injectionsignal[i] = 0;
    */

    // parse command line argument
    cmdparser(argc, argv, setupinfo);
  }

  // MPI_Type_create_struct
  // Create struct for setupinfo
  const int nitems = 6;
  int block[6] = {1, 1, 1, 1, MAXANT, MAXLEN};
  MPI_Datatype type[6] = {MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_FLOAT, MPI_INT, MPI_CHAR};
  MPI_Aint disp[6];
  MPI_Datatype structtype;

  disp[0] = offsetof(setup, verbose);
  disp[1] = offsetof(setup, test);
  disp[2] = offsetof(setup, seed);
  disp[3] = offsetof(setup, sfluxdensity);
  disp[4] = offsetof(setup, antSEFDs);
  disp[5] = offsetof(setup, inputfilename);

  MPI_Type_create_struct(nitems, block, disp, type, &structtype);
  MPI_Type_commit(&structtype);

  // broadcast setupinfo
  MPI_Bcast(&setupinfo, 1, structtype, MASTER, MPI_COMM_WORLD);

  //cout << "Process " << myid << ": " << setupinfo.seed << " " << setupinfo.inputfilename << endl;

  float tdur = 0.2 * 1e6;
  float testdur = 1; // time duration of test mode
  Configuration* config;
  Model* model;
  float dur;                            // duration of the simulated data in seconds
  float durus;                          // duration of the simulated data in microseconds
  float specRes, minStartFreq;
  int numdatastreams, numrecordedbands, freqindex;
  int numSamps;                         // number of samples to be generated per time block for the common signal
  int stime;                            // step time in microsecond == time block
  int configindex = 0;

  double timer = 0.0;
  double tt = 0.0;
  double procptrtime = 0.0;               
  double refvptime;                     // vptime is the time duration of the samples within the packet
  size_t framespersec;

  // initialize random number generator
  gsl_rng **rng_inst;
  gsl_rng_env_setup();

  rng_inst = (gsl_rng **) malloc(numprocs * sizeof(gsl_rng *));
  rng_inst[myid] = gsl_rng_alloc (gsl_rng_ranlux389);
  gsl_rng_set(rng_inst[myid],setupinfo.seed+myid);

  config = new Configuration(setupinfo.inputfilename, 0);
  model = config->getModel();

  // retrieve general information from config
  // if in test mode, only generate data for 1 second
  dur = setupinfo.test ? testdur : config->getExecuteSeconds();
  // dur is in second
  // durus is in microsecond
  durus = dur * 1e6;

  numdatastreams = config->getNumDataStreams();

  // use framespersec of the first antenna to calculate vptime
  // use this vptime as time reference
  framespersec = config->getFramesPerSecond(configindex, 0);
  refvptime = 1.0 * 1e6 / framespersec;
 
  bool sbbased = false;
  bool timebased = false;

  // print out information of the simulation
  // and retrieve information of the number of antenna and subbands

  // subbandsinfo is only used by MASTER to store sbinfo
  // and distribute the information to all processes using MPI_Scatter(...)
  // hope this will work

  // create a 2D array with the totoal number of subbands
  // containting antenna index and subband index of each subband
  int* subbandsinfo;
  int sbinfo[2] = {0, 0};
  // total subbands counter
  int sbcount = 0;
  if(myid == MASTER)
  {
    start = MPI_Wtime();
    // array with num-antenna number of elements
    // the value of which is the number of subbands the corresponding antenna has
    vector<int> antsb;
    if(!is_integer(refvptime))
    {
      cout << "VDIF packet time in microsecond is not an integer!! Something is wrong here ..." << endl;
      return EXIT_FAILURE;
    }

    cout << "Generate " << dur << " seconds data for " << numdatastreams << " stations ..." << "\n"
      << "source flux density is " << setupinfo.sfluxdensity << ".\n"
      << "random number generator seed is " << setupinfo.seed << endl;

    for(int i = 0; i < numdatastreams; i++)
    {
      framespersec = (size_t)config->getFramesPerSecond(configindex, i);
        
      numrecordedbands = config->getDNumRecordedBands(configindex, i);

      cout << "Telescope " << config->getTelescopeName(i) << "\n"
           << " Number of recorded band(s) is " << numrecordedbands << "\n"
           << " Antenna SEFD is " << setupinfo.antSEFDs[i] << "\n"
           << " Number of frames per second is " << framespersec << endl;
      
      // antenna subbands counter
      int antsbcnt = 0; 
      for(int j = 0; j < numrecordedbands; j++)
      {
        sbcount++;
        antsbcnt++;
        freqindex = config->getDRecordedFreqIndex(configindex, i, j); 
        cout << "Subband " << j << ":" << "\n"
             << "  Frequency " << config->getFreqTableFreq(freqindex) << "\n"
             << "  Bandwidth " << config->getFreqTableBandwidth(freqindex) << endl; 
      }
      antsb.push_back(antsbcnt);
    }

    // general information from model
    if(setupinfo.verbose >= 1)
    {
      cout << "Total number of subbands is  " << sbcount << endl;
      cout << "Number of scans is " << model->getNumScans() << endl;
      cout << "Scan start second is " << model->getScanStartSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
      cout << "Scan end second is " << model->getScanEndSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
    }

    subbandsinfo = new int [(sbcount+1) * 2];

    int numprocessed = 1;
    subbandsinfo[0] = 0;
    subbandsinfo[1] = 0;
    for(int i = 0; i < numdatastreams; i++)
    {
      numrecordedbands = config->getDNumRecordedBands(configindex, i);
      for(int j = 0; j < numrecordedbands; j ++)
      {
        int idx = 2 * (numprocessed + j);
        subbandsinfo[idx] = i;
        subbandsinfo[idx+1] = j;
      }
      numprocessed += numrecordedbands;
    }
    for(size_t idx = 0; idx < (size_t) (sbcount+1); idx++)
    {
      cout << subbandsinfo[idx*2] << " " << subbandsinfo[idx*2+1] << endl;
    }
  }
  // scatter sbinfo to each process
  // disbribute (antidx, sbidx) information to each process

  MPI_Scatter(&subbandsinfo[0], 2, MPI_INT, &sbinfo[0], 2, MPI_INT, MASTER, MPI_COMM_WORLD);
  
//==================================================
  // calculate specRes, number of samples per time block, step time
  if(getSpecRes(config, configindex, specRes, setupinfo.verbose) != EXIT_SUCCESS)
  {
    cout << "Process " << myid << ": Failed to calculate spectral resolution ..." << endl;
    return EXIT_FAILURE;
  };
  numSamps = getNumSamps(config, configindex, specRes, setupinfo.verbose);
  stime = static_cast<int>(1 / specRes); // step time in microsecond
  if(setupinfo.verbose >= 0 && myid == MASTER)
  {
    cout << "SpecRes is " << specRes << " MHz" << endl;
    cout << "number of samples per time block is " << numSamps << "\n"
         << "step time is " << stime << " us" << endl;
    cout << "tdur is " << tdur << endl;
  }
  minStartFreq = getMinStartFreq(config, configindex, setupinfo.verbose);

  // 'div' is the number of pieces the simulation time is divided into
  size_t div = 1;
  // use color for later implementation using MPI Shared-memeroy
  // Set it to zero for the moment
  int color = 0; 
  //int color = myid / (sbcount + 1);
  MPI_Comm local_comm;
  // If there are more subbands then number of cores, use time-baseed parallelization.
  // Otherwise, use a combination of subband- and time-based parallelization.

  if((sbcount+1) > numprocs)
  {
    cout << "Total number of cores is less than number of subbands plus 1." << endl;
    div = numprocs;
    //color = myid;
  }
  else
  {
    cout << "Use subband-based parallelization." << endl;
    div = numprocs / (sbcount + 1);
    //color = myid / (sbcount + 1);
    sbbased = true;
  }

  if(div != 1)
  {
    cout << "Use time-based parallelization." << endl;
    timebased = true;
    durus = durus/div;
    size_t remaining = (size_t)durus % div;
    if(remaining != 0)
      cout << "Cannot divide " << durus << " ms into " << div << " pieces.\n"
           << "Omitting the last " << remaining << " ms." << endl;
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, myid, &local_comm);
  int local_rank, local_size;
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);

  // master generates common signal
  // each subband copies its data to the right place
  // if use time-based parallelization
  // master is local_rank 0 of each group
  
  size_t stdur = tdur/stime;
  
  Subband* subband;
  int antidx = sbinfo[0];
  int sbidx = sbinfo[1];

  // each process initializes its corresponding subband
  size_t length, startIdx, blksize;
  size_t vpsamps;   // number of samples in a vdif packet
  size_t vpbytes;   // number of bytes in a single-thread vdif packet
  size_t framebytes;
  float freq, bw;
  string antname;
  int mjd, seconds;
  f64* tempcoeffs;
  f64* delaycoeffs;
  double antvptime;
  size_t antframespersec;

  framebytes = (size_t)config->getFrameBytes(configindex, antidx);
  numrecordedbands = config->getDNumRecordedBands(configindex, antidx);
  antframespersec = (size_t)config->getFramesPerSecond(configindex, antidx);
  antvptime = 1.0 * 1e6 / antframespersec;
  antname = config->getTelescopeName(antidx);
  // change the last character of the output vdif name to lower case for fourfit postprocessing
  //antname.back() = tolower(antname.back());
  antname.at(antname.size()-1) = tolower(antname.at(antname.size()-1));

  if(local_rank != MASTER)
  {
    // every antenna initialize its own subbands for antenna-based parallelization

    mjd = config->getStartMJD();
    seconds = config->getStartSeconds();
    if(setupinfo.verbose >= 1)
    {
      cout << "MJD is " << mjd << ", start seconds is " << seconds << endl;
    }

    // allocate memory for delaycoeffs
    tempcoeffs = vectorAlloc_f64(2);
    delaycoeffs = vectorAlloc_f64(2);

    if(setupinfo.verbose >= 1)
    { 
      cout << "Antenna " << antidx << endl;
      cout << " framebytes is " << framebytes << endl;
      cout << " numrecordedbands is " << numrecordedbands << endl;
      cout << " antenna name is " << antname << endl;
    }

    // only consider scan 0 source 0
    // scanindex, offsettime in seconds, timespan in seconds, numincrements, antennaindex, scansourceindex, order, delaycoeffs
    model->calculateDelayInterpolator(0, 0, antvptime*1e-6, 1, antidx, 0, 1, tempcoeffs);
    model->calculateDelayInterpolator(0, tempcoeffs[1]*1e-6, antvptime*1e-6, 1, antidx, 0, 1, delaycoeffs);
    if(setupinfo.verbose >= 2)
    {
      cout << "delay in us for datastream " << antidx << " at offsettime 0s is " << tempcoeffs[1] << endl;
      cout << "delay in us for datastream " << antidx << " at offsettime " << tempcoeffs[1]*1e-6 << "s " << " is " << delaycoeffs[1] << endl;
    }

    // calculate vdif packet size in terms of bytes and number of samples
    // each sample uses 4 bits, as the sample is complex and we use 2 bits sampling
    // therefore 2 samples per byte
    vpbytes = (framebytes - VDIF_HEADER_BYTES) / numrecordedbands + VDIF_HEADER_BYTES;
    vpsamps = (framebytes - VDIF_HEADER_BYTES) / numrecordedbands * 2;

    freq = config->getDRecordedFreq(configindex, antidx, sbidx);
    bw = config->getDRecordedBandwidth(configindex, antidx, sbidx);
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

    if(setupinfo.verbose >= 1)
    {
      cout << "Process: " << myid << ": Ant " << antidx << " subband " << sbidx << ":" << endl
           << "  start index is " << startIdx << endl 
           << "  block size is " << blksize << " length is " << length << endl
           << "  each vdif packet has " << vpsamps << " samples" << endl
           << "  VDIF_HEADER_BYTES is " << VDIF_HEADER_BYTES << " vpbytes is " << vpbytes << endl
           << "  number of samples in vdif packet is " << vpsamps << " framebytes is " << framebytes << endl
           << "  recorded freq is " << freq << endl;
    }

    subband = new Subband(startIdx, blksize, length, antidx, setupinfo.antSEFDs[antidx], sbidx, 
                       vpbytes, vpsamps, delaycoeffs, bw, antname, mjd, seconds, freq, setupinfo.verbose, color);

    // finish initializing subband
    // free memories for temporary allacated arrays
    vectorFree(tempcoeffs);
    vectorFree(delaycoeffs);

    if(setupinfo.verbose == 1)
      cout << "Process " << myid
           << ", local rank " << local_rank << " of group " << color << endl;
  }

  /*
   * generate common frequency domain signal
   * and send it to each subband of each antenna
   */
  if(local_rank == MASTER)
  {
    cout << "Master process of group " << color << " starts generating data ...\n"
    "The process may take a while, please be patient!" << endl;
    if(setupinfo.verbose >= 1)
        cout << "Master process of group " << color << " generates " << tdur << " us signal" << endl;
  }

  float* commFreqSig1;                  // 0.5 seconds common frequency domain signal
  float* commFreqSig2;
  //size_t workers = size_t(numprocs - 1);
  MPI_Request request[2];
  MPI_Status status[2];

  // allocate memory for the common frequency domain signal
  int sampsize = numSamps*2*stdur;
  commFreqSig1 = new float [sampsize];
  commFreqSig2 = new float [sampsize];

  if(local_rank == MASTER)
  {
    gencplx(commFreqSig1, sampsize, STDEV, rng_inst[myid], setupinfo.verbose);
    MPI_Ibcast(commFreqSig1, sampsize, MPI_FLOAT, MASTER, local_comm, &request[0]);
    // User tdur + 1 as reference to calculate the smallest process pointer time 
    procptrtime = tdur + 1;
    gencplx(commFreqSig2, sampsize, STDEV, rng_inst[myid], setupinfo.verbose);
  }
  else
    MPI_Ibcast(commFreqSig1, sampsize, MPI_FLOAT, MASTER, local_comm, &request[0]);
  MPI_Wait(&request[0], &status[0]);

  if(local_rank != MASTER)
  {
    subband->fabricatedata(commFreqSig1, rng_inst[myid], setupinfo.sfluxdensity);

    // move data in each array from the second half to the first half
    // and set the process pointer to the proper location
    // i.e. data is moved half array ahead, therefore process pointer 
    // is moved half array ahead
    movedata(subband, setupinfo.verbose);

    // after TDUR time signal is generated for each subband array
    // set the current pointer of each array back to the beginning of the second half

    subband->setcptr(subband->getlength() / 2);
    if(setupinfo.verbose >= 2) 
      cout << "Antenna " << subband->getantIdx() << " subband " << subband->getsbIdx()
                         << " set current pointer back to " << subband->getlength() / 2 << endl;
    // each subband calculates its own process pointer time
    procptrtime = subband->getprocptr() * (1.0 / subband->getbandwidth());
    if(setupinfo.verbose >= 2) cout << "Process " << myid << ": procptrtime is " << procptrtime << endl;
  }

  MPI_Allreduce(&procptrtime, &tt, 1, MPI_DOUBLE, MPI_MIN, local_comm);
  if((local_rank == MASTER) && (tt >= tdur))
  {
    cout << "**the lowest process pointer time is larger than tdur!!!\n"
            "**Something is wrong here!!!" << endl;
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    return (EXIT_FAILURE);
  }


  do
  {
    MPI_Ibcast(commFreqSig2, sampsize, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &request[1]);
    while((tt < tdur) && (timer < durus))
    {
      if(local_rank != MASTER)
      {
        // process and packetize one vdif packet for each subband array
        int rc = processAndPacketize(antframespersec, subband, model, setupinfo.verbose);
        if(rc)
        {
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          return(EXIT_FAILURE);
        }
      }
      // update tt
      tt += refvptime;
      timer += refvptime;
    }
    MPI_Wait(&request[1], &status[1]);

    if(local_rank != MASTER)
    {
      subband->fabricatedata(commFreqSig2, rng_inst[myid], setupinfo.sfluxdensity);

      // move data in each array from the second half to the first half
      // and set the process pointer to the proper location
      // i.e. data is moved half array ahead, therefore process pointer 
      // is moved half array ahead
      movedata(subband, setupinfo.verbose);

      // after TDUR time signal is generated for each subband array
      // set the current pointer of each array back to the beginning of the second half

      subband->setcptr(subband->getlength() / 2);
      if(setupinfo.verbose >= 2) 
        cout << "Antenna " << subband->getantIdx() << " subband " << subband->getsbIdx()
                           << " set current pointer back to " << subband->getlength() / 2 << endl;
      // each subband calculates its own process pointer time
      procptrtime = subband->getprocptr() * (1.0 / subband->getbandwidth());
      if(setupinfo.verbose >= 2) cout << "Process " << myid << ": procptrtime is " << procptrtime << endl;

      // receive data in the second half of the array
      // process and packetize at the same time
      MPI_Ibcast(commFreqSig1, sampsize, MPI_FLOAT, MASTER, local_comm, &request[0]);
    }
    else
    {
      gencplx(commFreqSig1, sampsize, STDEV, rng_inst[myid], setupinfo.verbose);
      MPI_Ibcast(commFreqSig1, sampsize, MPI_FLOAT, MASTER, local_comm, &request[0]);
      procptrtime = tdur + 1;
      gencplx(commFreqSig2, sampsize, STDEV, rng_inst[myid], setupinfo.verbose);
    }

    MPI_Allreduce(&procptrtime, &tt, 1, MPI_DOUBLE, MPI_MIN, local_comm);
    if((local_rank == MASTER) && (tt >= tdur))
    {
      cout << "**the lowest process pointer time is larger than tdur!!!\n"
              "**Something is wrong here!!!" << endl;
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      return (EXIT_FAILURE);
    }

    while((tt < tdur) && (timer < durus))
    {
      if(local_rank != MASTER)
      {
        // process and packetize one vdif packet for each subband array
        int rc = processAndPacketize(antframespersec, subband, model, setupinfo.verbose);
        if(rc)
        {
          MPI_Abort(MPI_COMM_WORLD, ERROR);
          return(EXIT_FAILURE);
        }
      }
      // update tt
      tt += refvptime;
      timer += refvptime;
    }
    MPI_Wait(&request[0], &status[0]);

    if((local_rank == MASTER) && (setupinfo.verbose >=2))
      cout << "tt is " << tt << ", timer is " << timer << endl;

  } while(timer < durus);


  // free allocated common signal memory
  delete commFreqSig1;
  delete commFreqSig2;

  if(local_rank != MASTER)
  {
    subband->closevdif();
    // free allocated subband memory
    delete subband;
  }

  if(myid== MASTER)
  {
    end = MPI_Wtime();
    elapse = (end - start)/60.0; // convert time from seconds to minutes
  
    cout << "All data has been generated successfully, bye!" << endl;
    cout << "Total duration is " << elapse << " minutes." << endl;
  }

  if(local_rank < numdatastreams)
  {
    // combine VDIF files
    vdifzipper(config, configindex, durus, setupinfo.verbose, local_rank, color);
  }


  MPI_Type_free(&structtype);
  MPI_Comm_free(&local_comm);


  // free random number generator
  gsl_rng_free(rng_inst[myid]);
  MPI_Finalize();

  return(EXIT_SUCCESS);
}
