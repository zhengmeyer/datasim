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
  
  #define EPSILON 1e-6

  /* parameters used for quantization and vdif */
  #define CORR 0.05               // amplitude of common signal relative to the station receiver noise
  #define BSMX 1024 * 1024
  #define BITS 2
  #define ISCOMPLEX 0
  #define BITSPERBYTE 8

#endif /* __DATASIM_H__ */

/*
 * eof
 */
