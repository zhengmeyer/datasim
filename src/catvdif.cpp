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
#include <fstream>
#include <string>
#include <sstream>
#include <stdint.h>
#include "architecture.h"
#include "configuration.h"
#include "vdifio.h"
#include "catvdif.h"
#include <mpi.h>

using namespace std;


// Read in vdif header produced after vdifzipper

void catvdif(Configuration* config, int configindex, float durus, size_t verbose, int myid, size_t div)
{
  int mjd, seconds;
  
  cout << "Combine VDIF files of each time-segment into a single-thread multi-channel VDIF file ..." << endl;
 
  mjd = config->getStartMJD();
  seconds = config->getStartSeconds();
  if(verbose >= 1)
  {
    cout << "mjd is " << mjd << ", seconds is " << seconds << endl;
  }

  // each antenna combines its time-segment into one VDIF file

  size_t framebytes;                    // data frame size of the output VDIF file
  size_t numrecordedbands;              // number of channels
  float bw;                             // bandwidth
  size_t framespersec, totalnumframes;
  string antname;
  ofstream outputvdif;

  if(verbose >= 1)
  {
    cout << "Antenna " << myid << ":" << endl;
  }
  // retrieve number of subbands, total framebytes, antenna name
  // and the vdif packet bytes for each singal channel vdif file
  framebytes = (size_t)config->getFrameBytes(configindex, myid);
  numrecordedbands = (size_t)config->getDNumRecordedBands(configindex, myid);
  antname = config->getTelescopeName(myid);
  // change the last character of the output vdif name to lower case for fourfit postprocessing
  //antname.back() = tolower(antname.back());
  antname.at(antname.size()-1) = tolower(antname.at(antname.size()-1));

  ifstream tsegfile[div];

  // retrieve bandwidth information and number of frames per second of each antenna
  bw = config->getDRecordedBandwidth(configindex, myid, 0);
  framespersec = config->getFramesPerSecond(configindex, myid);
  totalnumframes = durus / 1e6 * framespersec;

  if(verbose >= 1)
  {
    cout << " framebyte is " << framebytes << "bytes, number of channels is " << numrecordedbands<< "\n"
         << " bandwitdh is " << bw << "MHz, framespersec is " << framespersec<< "\n"
         << " total number of frames is " << totalnumframes << endl;
  }
  // allocate memory for vdif packet buffer and input file streams
  uint8_t* outputvdifbuf;
  uint8_t* inputvdifbuf;
  inputvdifbuf = vectorAlloc_u8(framebytes);
  outputvdifbuf = vectorAlloc_u8(framebytes);
  fill_n(inputvdifbuf, framebytes, 0);
  fill_n(outputvdifbuf, framebytes, 0);
  if(verbose >= 2)
  {
    cout << " Allocated memory for vdif packet buffer for antenna " << myid << endl;  
  }

  // initialize VDIF header of the output buffer
  createVDIFHeader((vdif_header *)outputvdifbuf, framebytes - VDIF_HEADER_BYTES, myid, BITS, numrecordedbands, ISCOMPLEX, (char *)antname.c_str());
  setVDIFEpoch((vdif_header *)outputvdifbuf, mjd);
  setVDIFFrameMJDSec((vdif_header *)outputvdifbuf, mjd*86400 + seconds);
  if(verbose >= 2)
  {
    cout << " VDIF header initialized" << endl; 
  }

  try
  {
    // open multi-channel vdif file to write to
    outputvdif.open((antname + ".vdif").c_str(), ios::binary);
    if(verbose >= 2)
    {
      cout << " Opened " << antname << ".vdif for writing ..." << endl;
    }
    // open vdif file of each time segment
    for(size_t  tseg = 0; tseg < div; tseg++)
    {
      stringstream ss;  
      ss << tseg;
      ostringstream filename;
      filename << antname << "-" << ss.str() << ".vdif";
      tsegfile[tseg].open(filename.str().c_str(), ios::binary);
      if(verbose >= 2)
      {
        cout << " Opened input file " << filename.str() << endl;         
      }

	    for(size_t idx = 0; idx < totalnumframes; idx++)
	    {
	      // reset the value of the input vdif buffer
	      fill_n(inputvdifbuf, framebytes, 0);

	      tsegfile[tseg].read((char *)inputvdifbuf, framebytes * sizeof(uint8_t));
	      for(size_t byte = VDIF_HEADER_BYTES; byte < framebytes; byte++)
	      	outputvdifbuf[byte] = inputvdifbuf[byte];

	      // write to output vdif file
	      outputvdif.write((char *)outputvdifbuf, framebytes * sizeof(uint8_t));

	      //update VDIF Header for the next packet
	      if(idx != totalnumframes - 1)
	        nextVDIFHeader((vdif_header *)outputvdifbuf, (int) framespersec);     
	    }
	    // close input vdif files
	    for(size_t ch = 0; ch < numrecordedbands; ch++)
	    {
	      tsegfile[tseg].close();
	      if(verbose >= 1)
	      { 
	        cout << " Closed input file for time-segment " << tseg << endl;
	      }
	    }
	  }
    // close multi-channel vdif file after writing
    outputvdif.close();
    if(verbose >= 1)
    {
      cout << " Closed " << antname << ".vdif ..." << endl;
    }
  } catch (ofstream::failure e) {                                                                              
    cerr << "Exception opening/closinging input or output vdif file" << endl;                                              
  }
  // free memory of input and output vdif buffer
  vectorFree(inputvdifbuf);
  vectorFree(outputvdifbuf);
  if(verbose >= 2)
  {
    cout << "Freed memory allocated for input and output vdif buffers" << endl;
  }
}
