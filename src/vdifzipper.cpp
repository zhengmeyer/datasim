/*
 * vdifzipper.cpp
 * The Non-zero baseline data simulator
 * Combine sigle-channel VDIF files to a single thread multi-channel VDIF file
 *
 * Author: Zheng Meyer-Zhao
 * 2014/02/21
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdint>
#include "architecture.h"
#include "configuration.h"
#include "vdifio.h"
#include "vdifzipper.h"

using namespace std;
/**
 * @config Configuration file
 * @configindex Configuration index
 * @durus Observation time in microseconds
 */
void vdifzipper(Configuration* config, int configindex, float durus, size_t verbose)
{
  size_t numdatastreams = (size_t)config->getNumDataStreams();
  int mjd, seconds;
  
  cout << "Combine VDIF files of each channel into a single-thread multi-channel VDIF file ..." << endl;
 
  mjd = config->getStartMJD();
  seconds = config->getStartSeconds();
  if(verbose >= 1)
  {
    cout << "mjd is " << mjd << ", seconds is " << seconds << endl;
  }

  // loop through all the antennas
  for(size_t i = 0; i < numdatastreams; i++)
  {
    size_t framebytes, numrecordedbands;
    size_t chvpbytes, sampsperbyte = 8 / BITS;
    size_t osampsperbyte;
    float bw;
    size_t framespersec, totalnumframes;
    size_t shift, ishift, oshift;
    uint8_t *optr, *iptr;
    uint8_t bits, imask, omask;
    string antname;
    ofstream outputvdif;

    if(verbose >= 1)
    {
      cout << "Antenna " << i << ":" << endl;
    }
    // retrieve number of subbands, total framebytes, antenna name
    // and the vdif packet bytes for each singal channel vdif file
    framebytes = (size_t)config->getFrameBytes(configindex, i);
    numrecordedbands = (size_t)config->getDNumRecordedBands(configindex, i);
    antname = config->getTelescopeName(i);
    // change the last character of the output vdif name to lower case for fourfit postprocessing
    antname.back() = tolower(antname.back());

    chvpbytes = (framebytes - VDIF_HEADER_BYTES) / numrecordedbands + VDIF_HEADER_BYTES;

    ifstream chfile[numrecordedbands];
    // retrieve bandwidth information of each antenna
    // calculate vdif packet time and number of frames per second
    bw = config->getDRecordedBandwidth(configindex, i, 0);
    framespersec = config->getFramesPerSecond(configindex, i);
    totalnumframes = durus / 1e6 * framespersec;
    osampsperbyte = sampsperbyte / numrecordedbands;

    if(verbose >= 1)
    {
      cout << " framebyte is " << framebytes << "bytes, number of channels is " << numrecordedbands<< "\n"
           << " bandwitdh is " << bw << "MHz, framespersec is " << framespersec<< "\n"
           << " total number of frames is " << totalnumframes << endl;
    }
    // allocate memory for vdif packet buffer and input file streams
    uint8_t* outputvdifbuf;
    uint8_t* inputvdifbuf;
    outputvdifbuf = vectorAlloc_u8(framebytes);
    fill_n(outputvdifbuf, framebytes, 0);
    inputvdifbuf = vectorAlloc_u8(chvpbytes);
    fill_n(inputvdifbuf, chvpbytes, 0);
    if(verbose >= 2)
    {
      cout << " Allocated memory for vdif packet buffer for antenna " << i << endl;  
    }

    // initialize VDIF header of the output buffer
    createVDIFHeader((vdif_header *)outputvdifbuf, framebytes - VDIF_HEADER_BYTES, i, BITS, numrecordedbands, ISCOMPLEX, (char *)antname.c_str());
    setVDIFEpoch((vdif_header *)outputvdifbuf, mjd);
    setVDIFFrameMJDSec((vdif_header *)outputvdifbuf, mjd*86400 + seconds);
    if(verbose >= 2)
    {
      cout << " VDIF header initialized" << endl; 
    }

    try {
      // open multi-channel vdif file to write to
      outputvdif.open(antname + ".vdif", ios::binary);
      if(verbose >= 2)
      {
        cout << " Opened " << antname << ".vdif for writing ..." << endl;
      }
      // open vdif file of each channel 
      for(size_t ch = 0; ch < numrecordedbands; ch++)
      {
        stringstream ss;  
        ss << ch;
        ostringstream filename;
        filename << antname << "_" << ss.str() << ".vdif";
        chfile[ch].open(filename.str(), ios::binary);
        if(verbose >= 2)
        {
          cout << " Opened input file " << filename.str() << endl;         
        }
      }

      for(size_t idx = 0; idx < totalnumframes; idx++)
      {
        optr = &outputvdifbuf[VDIF_HEADER_BYTES];
        // loop through all the channels
        for(size_t ch = 0; ch < numrecordedbands; ch++)
        {
          optr = &outputvdifbuf[VDIF_HEADER_BYTES];
          // reset the value of the input vdif buffer
          fill_n(inputvdifbuf, chvpbytes, 0);

          chfile[ch].read((char *)inputvdifbuf, chvpbytes * sizeof(uint8_t));

          // set input vdif buffer pointer at proper location
          iptr = &inputvdifbuf[VDIF_HEADER_BYTES];

          // loop through the data bytes of the channel
          for(size_t bytes = VDIF_HEADER_BYTES; bytes < chvpbytes; bytes++)
          {
            if((idx == 0) && (bytes < VDIF_HEADER_BYTES + 4) && (verbose >= 2))
              cout << "input bytes " << bytes <<": value " << static_cast<unsigned>(*iptr) << endl;
            //*optr = (*iptr);
            size_t osampcnt = 0;
            for(size_t samp = 0; samp < sampsperbyte; samp++)
            {
              // shift to the proper input bit
              ishift = samp * BITS;
              imask = 03;
              imask <<= ishift;
              bits = (*iptr) & imask;
              bits >>= ishift;
              // store the value at the proper output bits
              shift = ishift * numrecordedbands + ch * 2;
              oshift = shift % 8;
              bits <<= oshift;
              omask = 03;
              omask <<= oshift;
              (*optr) &= ~omask;
              (*optr) |= bits;
              if((idx == 0) && (bytes < VDIF_HEADER_BYTES + 4) && (verbose >= 2)) {
                cout << "ishift is " << static_cast<unsigned>(ishift) << ", imask is " << static_cast<unsigned>(imask) << "\n"
                     << "bits is " << static_cast<unsigned>(bits) << ", shift is " << static_cast<unsigned>(shift) << "\n"
                     << "oshift is " << static_cast<unsigned>(oshift) << ", omask is " << static_cast<unsigned>(omask) << endl;
              }
              // update output sample counter
              osampcnt++;
              // if the current output sample byte is full
              // move the pointer to the next byte
              if(osampcnt == osampsperbyte) {
                if((idx == 0) && (bytes < VDIF_HEADER_BYTES + 4) && (verbose >= 2))
                  cout << "output bytes value is " << static_cast<unsigned>(*optr) << endl; 
                optr++;
                osampcnt = 0;
              }
            }
            iptr++;
            //optr++;
          }
        }
        // write to output vdif file
        outputvdif.write((char *)outputvdifbuf, framebytes * sizeof(uint8_t));

        //update VDIF Header for the next packet
        if(idx != totalnumframes - 1)
          nextVDIFHeader((vdif_header *)outputvdifbuf, (int) framespersec);     
      }
      // close input vdif files
      for(size_t ch = 0; ch < numrecordedbands; ch++)
      {
        chfile[ch].close();
        if(verbose >= 1)
        { 
          cout << " Closed input file for channel " << ch << endl;
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
    vectorFree(outputvdifbuf);
    vectorFree(inputvdifbuf);
    if(verbose >= 2)
    {
      cout << "Freed memory allocated for input and output vdif buffers" << endl;
    }
  } 
}

//int main(int argc, char* argv[])
//{
//  Configuration* config;
//  int configindex = 0;
//  float dur, durus;
//  size_t verbose = 0;
//
//  config = new Configuration(argv[1], configindex);
//  // retrieve observation time in seconds
//  // and convert it to microseconds
//  dur = config->getExecuteSeconds();
//  durus = dur * 1e6;
//  vdifzipper(config, configindex, durus, verbose);
//
//  return 0;
//}
