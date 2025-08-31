/* ------------------------------------------------------------------------- *\

    Declarations used in both eye tracking and calibraiton
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>


#ifndef min
#define min(x,y) ((x)<(y) ? (x):(y))
#endif
#ifndef max
#define max(x,y) ((x)>(y) ? (x):(y))
#endif
#define F_PI ((float)M_PI)

#define MAX_EL_RATIO    0.9995f     // <1 to help with phi direction in 1st frame

// utility functions used both in tracking and calibration
float norm_euc2(float *v1, float *v2);
void  rotationMat(float *angles, float mat[3][3], int inv);
void  rotcam(float crd[2], float ang, float out[2]);
void  compensatePupilElongation(float *xy, int xylen, float el_r, float pup_ang, float pup_stretch);
float camRayProjToMinor(float xyc[2], float xye[2], float iris_z);
void  compensatePCA(float vec_to_eyeb[2], float iris_ratio, float pca_angle[2], float *el_r, float *phi_minor);
