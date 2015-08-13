/*
 * datasim.h
 * VLBI data simulator
 *
 * Author: Zheng Meyer-Zhao
 * Date: 2013/11/04
 */

#ifndef __DATASIM_H__
#define __DATASIM_H__
  
  #include <cstdint>
  #include <string>

  #define EPSILON 1e-6

  /* parameters used for quantization and vdif */
  #define CORR 0.05               // amplitude of common signal relative to the station receiver noise
  #define BSMX 1024 * 1024
  #define BITS 2
  #define ISCOMPLEX 0
  #define BITSPERBYTE 8

  typedef struct setup {
    int verbose;
    int test;                     // test mode
    float sfluxdensity;           // source flux density in Jansky
    std::string inputfilename;    // .input file name
  } setup;

#endif /* __DATASIM_H__ */

/*
 * eof
 */
