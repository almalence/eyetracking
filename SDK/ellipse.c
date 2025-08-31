/* ------------------------------------------------------------------------- *\

    Ellipse fitting functions
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

    Notes:
    - EllipseDirectFit() function code created from Matlab version provided under the license below

\* ------------------------------------------------------------------------- */

#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>

#include "ellipse.h"
#include "eigen.h"

#define F_PI ((float)M_PI)

//  Direct ellipse fit, proposed in article
//    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
//     "Direct Least Squares Fitting of Ellipses"
//     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
//
//  Our code is based on a numerically stable version
//  of this fit published by R. Halir and J. Flusser
//
//     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
//
//     Output: A = [a b c d e f]' is the vector of algebraic 
//             parameters of the fitting ellipse:
//             ax^2 + bxy + cy^2 +dx + ey + f = 0
//             the vector A is normed, so that ||A||=1
//
//  This is a fast non-iterative ellipse fit.
//
//  It returns ellipses only, even if points are
//  better approximated by a hyperbola.
//  It is somewhat biased toward smaller ellipses.
//
// Copyright (c) 2009, Nikolai Chernov
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in
//       the documentation and/or other materials provided with the distribution
//     * Neither the name of the University of Alabama at Birmingham nor the names
//       of its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//
int EllipseDirectFit(float *xy, int xylen, int skip_st, int skip_en, float A[6])
{
    int i,j,k;

    int npix = xylen - (skip_en-skip_st);

    // the centroid of the data set
    float cx=0, cy=0;
    for (i=0; i<xylen; ++i)
    {
        if ((i>=skip_st) && (i<skip_en)) continue;
        cx += xy[i*2];
        cy += xy[i*2+1];
    }
    cx /= npix;
    cy /= npix;

    float S1[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    float S2[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    float S3[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};

    // ToDo: S1, S3 and iS3 seem to be symmetric, reduce the number of operations
    for (i=0; i<xylen; ++i)
    {
        float D1[3];
        float D2[3];

        if ((i>=skip_st) && (i<skip_en)) continue;

        D1[0] = (xy[i*2]  -cx) * (xy[i*2]  -cx);
        D1[1] = (xy[i*2]  -cx) * (xy[i*2+1]-cy);
        D1[2] = (xy[i*2+1]-cy) * (xy[i*2+1]-cy);
        D2[0] = xy[i*2]  -cx;
        D2[1] = xy[i*2+1]-cy;
        D2[2] = 1.0f;

        for (j=0; j<3; ++j)
            for (k=0; k<3; ++k)
            {
                S1[j][k] += D1[j]*D1[k];
                S2[j][k] += D1[j]*D2[k];
                S3[j][k] += D2[j]*D2[k];
            }
    }

    // inversion of matrix S3
    float det = S3[0][0]*S3[1][1]*S3[2][2] + S3[0][1]*S3[1][2]*S3[2][0] + S3[0][2]*S3[1][0]*S3[2][1] -
                S3[0][0]*S3[1][2]*S3[2][1] - S3[0][1]*S3[1][0]*S3[2][2] - S3[0][2]*S3[1][1]*S3[2][0];
    if (det==0) return 0; // can not invert matrix
    float iS3[3][3];      // inverse of S3: adjoint matrix divided by determinant
    iS3[0][0]= (S3[1][1]*S3[2][2]-S3[2][1]*S3[1][2])/det;  iS3[1][0]=-(S3[1][0]*S3[2][2]-S3[1][2]*S3[2][0])/det;  iS3[2][0]= (S3[1][0]*S3[2][1]-S3[1][1]*S3[2][0])/det;
    iS3[0][1]=-(S3[0][1]*S3[2][2]-S3[2][1]*S3[0][2])/det;  iS3[1][1]= (S3[0][0]*S3[2][2]-S3[0][2]*S3[2][0])/det;  iS3[2][1]=-(S3[0][0]*S3[2][1]-S3[0][1]*S3[2][0])/det;
    iS3[0][2]= (S3[0][1]*S3[1][2]-S3[1][1]*S3[0][2])/det;  iS3[1][2]=-(S3[0][0]*S3[1][2]-S3[0][2]*S3[1][0])/det;  iS3[2][2]= (S3[0][0]*S3[1][1]-S3[0][1]*S3[1][0])/det;

    float T[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    for (j=0; j<3; ++j)
        for (k=0; k<3; ++k)
            T[k][j] += -iS3[k][0]*S2[j][0] -iS3[k][1]*S2[j][1] -iS3[k][2]*S2[j][2];

    // reusing S1 to keep M matrix calculations (also symmetric here)
    for (j=0; j<3; ++j)
        for (k=0; k<3; ++k)
            S1[k][j] += S2[k][0]*T[0][j] + S2[k][1]*T[1][j] + S2[k][2]*T[2][j];

    float M[3][3];
    M[0][0] =  S1[2][0]/2;  M[0][1] =  S1[2][1]/2;  M[0][2] =  S1[2][2]/2;
    M[1][0] = -S1[1][0];    M[1][1] = -S1[1][1];    M[1][2] = -S1[1][2];
    M[2][0] =  S1[0][0]/2;  M[2][1] =  S1[0][1]/2;  M[2][2] =  S1[0][2]/2;

    float eval[3];
    float evec[3][3];
    eigen_solve(M, 0.0001f, 8, eval, evec);

    float cond[3];
    for (i=0; i<3; ++i)
        cond[i] = 4*evec[0][i]*evec[2][i]-evec[2][i]*evec[1][i];

    int poscond=0;
    if (cond[1]>cond[0]) poscond=1;
    if (cond[2]>cond[poscond]) poscond=2;
    if (cond[poscond]<=0) return 0;

    A[0] = evec[0][poscond];
    A[1] = evec[1][poscond];
    A[2] = evec[2][poscond];

    A[3] = evec[0][poscond]*T[0][0] + evec[1][poscond]*T[0][1] + evec[2][poscond]*T[0][2];
    A[4] = evec[0][poscond]*T[1][0] + evec[1][poscond]*T[1][1] + evec[2][poscond]*T[1][2];
    A[5] = evec[0][poscond]*T[2][0] + evec[1][poscond]*T[2][1] + evec[2][poscond]*T[2][2];

    A[5] += A[0]*cx*cx + A[1]*cx*cy + A[2]*cy*cy - A[3]*cx - A[4]*cy;
    A[3] -= 2*A[0]*cx + A[1]*cy;
    A[4] -= 2*A[2]*cy + A[1]*cx;

    float normA = sqrtf(A[0]*A[0] + A[1]*A[1] + A[2]*A[2] + A[3]*A[3] + A[4]*A[4] + A[5]*A[5]);
    for (i=0; i<6; ++i)
        A[i] = A[i]/normA;

    return 1;
}


int EllipseGeomFromAlgebraic(float A[6], float *cx, float *cy, float *semx, float *semy, float *phi)
{
    float a = A[0];
    float b = A[1]/2;
    float c = A[2];
    float d = A[3]/2;
    float f = A[4]/2;
    float g = A[5];

    *cx = (c*d - b*f)/(b*b-a*c);
    *cy = (a*f - b*d)/(b*b-a*c);
    
    float s1 = (a-c)*(a-c)+4*b*b;
    if (s1<0) return 0;
    float s2 = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g);
    float s3 = (b*b-a*c)*( sqrtf(s1)-(a+c));
    float s4 = (b*b-a*c)*(-sqrtf(s1)-(a+c));
    if (s2/s3<0 || s2/s4<0) return 0;
    *semx = sqrtf( s2 / s3);
    *semy = sqrtf( s2 / s4);

    *phi = (atan2f(2*b, a-c)+F_PI)/2;
    if (*phi>3*F_PI/4) *phi -= F_PI;

    return 1;
}
