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
  #define CORR 0.05
  #define BSMX 1024 * 1024
  #ifndef BITS
    #define BITS 2
  #endif
  
  #ifndef ISCOMPLEX
    #define ISCOMPLEX 0
  #endif

#endif /* __DATASIM_H__ */

/*
 * eof
 */
