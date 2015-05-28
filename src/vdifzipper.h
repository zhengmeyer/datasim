/*
 * vdifzipper.h
 * VLBI data simulator
 *
 * Author: Zheng Meyer-Zhao
 * 2014/02/27
 *
 * Utilities used by the data simulator
 */

#ifndef __VDIFZIPPER_H__
#define __VDIFZIPPER_H__

#include "configuration.h"

#ifndef BITS
  #define BITS 2
#endif

#ifndef ISCOMPLEX
  #define ISCOMPLEX 0
#endif

void vdifzipper(Configuration* config, int configindex, float durus, size_t verbose);

#endif /* __VDIFZIPPER_H__ */

/*
 * eof
 */
