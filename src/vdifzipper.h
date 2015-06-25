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
#include "datasim.h"

void vdifzipper(Configuration* config, int configindex, float durus, size_t verbose);

#endif /* __VDIFZIPPER_H__ */

/*
 * eof
 */
