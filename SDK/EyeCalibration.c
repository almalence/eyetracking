/* ------------------------------------------------------------------------- *\

    Calibration of eye parameters
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <memory.h>

#include "EyeTracking.h"
#include "EyeUtil.h"
#include "ellipse.h"

#define N_DIAGS      0 // 4              // set to 0 for 5-point calibration, to 4 for 9-point calibration
#define N            (5+N_DIAGS)    // number of calibration locations, sequence: Left Right Top Bottom Diags_if_any Center

#if N_DIAGS==0
#define DIAGS_FILL0
#define DIAGS_FILL1
#else
#define DIAGS_FILL0   0,0,0,0,
#define DIAGS_FILL1   1,1,1,1,
#endif

#define MAX_PUPIL_XY 256            // pupil edge can be defined by that many [x,y] pairs (haven't seen over 100 in practice)


void estimateCtr(float p2d[N][2], float gt[N][2], float ctrp[N][2])
{
    for (int i=0; i<N; ++i)
    {
        float d_ref[2] = {0};
        float d_cam[2] = {0};
        for (int j=0; j<N; ++j)
        {
            d_ref[0] += fabsf( gt[i][0] - gt[j][0] );
            d_ref[1] += fabsf( gt[i][1] - gt[j][1] );
            d_cam[0] += fabsf( p2d[i][0] - p2d[j][0] );
            d_cam[1] += fabsf( p2d[i][1] - p2d[j][1] );
        }
    
        ctrp[i][0] = p2d[i][0] - gt[i][0] * d_cam[0] / d_ref[0];
        ctrp[i][1] = p2d[i][1] - gt[i][1] * d_cam[1] / d_ref[1];
    }
}

// find eyeball center on image plane and slippage of each location
// Note: camera perspective projection not accounted for here
void getCtrAndSlip(float p2d[N][2], float gt[N][2], float ctr[2], float slippage[N][2])
{
    float ctr_gt[2] = {0};
    float ctr_p2d[2] = {0};
    for (int i=0; i<N; ++i)
    {
        ctr_gt[0]  += gt[i][0];
        ctr_gt[1]  += gt[i][1];
        ctr_p2d[0] += p2d[i][0];
        ctr_p2d[1] += p2d[i][1];
    }
    ctr_gt[0]  /= N; ctr_gt[1]  /= N;
    ctr_p2d[0] /= N; ctr_p2d[1] /= N;

    float scl_gt[2] = {0};
    float scl_p2d[2] = {0};
    for (int i=0; i<N; ++i)
    {
        scl_gt[0]  += fabsf(gt[i][0] -ctr_gt[0]);
        scl_gt[1]  += fabsf(gt[i][1] -ctr_gt[1]);
        scl_p2d[0] += fabsf(p2d[i][0]-ctr_p2d[0]);
        scl_p2d[1] += fabsf(p2d[i][1]-ctr_p2d[1]);
    }

    float p2d_ideal[N][2];
    for (int i=0; i<N; ++i)
    {
        p2d_ideal[i][0] = ctr_p2d[0] + (gt[i][0]-ctr_gt[0]) * scl_p2d[0] / scl_gt[0];
        p2d_ideal[i][1] = ctr_p2d[1] + (gt[i][1]-ctr_gt[1]) * scl_p2d[1] / scl_gt[1];

        slippage[i][0] = p2d_ideal[i][0] - p2d[i][0];
        slippage[i][1] = p2d_ideal[i][1] - p2d[i][1];
    }

    float ctrp[N][2];
    estimateCtr(p2d_ideal, gt, ctrp);

    ctr[0] = ctr[1] = 0;
    for (int i=0; i<N; ++i)
    {
        ctr[0] += ctrp[i][0];
        ctr[1] += ctrp[i][1];
    }
    ctr[0] /= N;
    ctr[1] /= N;
}


void encl_rect(float semx, float semy, float phi, float *rw, float *rh)
{
    float ux = semx * cosf(phi);
    float uy = semx * sinf(phi);
    float vx = semy * cosf(phi + F_PI/2);
    float vy = semy * sinf(phi + F_PI/2);
    *rw = 2*sqrtf(ux*ux + vx*vx);
    *rh = 2*sqrtf(uy*uy + vy*vy);
}


// ToDo: for off-axis camera, ratio at the central position will be not 1, and we may need to take direction-specific component of it
// Looks like it will be an issue for camera_rot ~= 0 situations only
float diffFromPCA(float (*pxy)[MAX_PUPIL_XY][2], int *nxy, float (*cxy)[2], int anc_st[N], int anc_en[N], float gv_cam[N][3], int use[N], float ctr[2], float img_ctr[2], float iris_z, float ax, float ay)
{
    float ratio[N] = {0};

    for (int i=0; i<N; ++i)
        if (use[i])
        {
            float rect_true_n = 0;
            float rect_obs_n = 0;

            for (int f=anc_st[i]; f<anc_en[i]; ++f)
            {
                float rw, rh;
                float A[6];

                float *xy = (float *)(pxy[f]);
                if (!EllipseDirectFit(xy, nxy[f], 0, 0, A)) return 1e10f;
                float tmpcx, tmpcy, semx, semy, phi;
                EllipseGeomFromAlgebraic(A, &tmpcx, &tmpcy, &semx, &semy, &phi);

                // find enclosing rectangle for the observed ellipse
                encl_rect(semx, semy, phi, &rw, &rh);
                // observed bounding rectangle ratio at given position
                rect_obs_n = rect_obs_n + rw/rh;

                // reverse-compensate PCA
                float pca_ang[3] = {-ay, -ax, 0};
                float pca_mat[3][3];
                rotationMat(pca_ang, pca_mat, 0);
                float gv_pca[3];
                gv_pca[0] = pca_mat[0][0]*gv_cam[i][0] + pca_mat[0][1]*gv_cam[i][1] + pca_mat[0][2]*gv_cam[i][2];
                gv_pca[1] = pca_mat[1][0]*gv_cam[i][0] + pca_mat[1][1]*gv_cam[i][1] + pca_mat[1][2]*gv_cam[i][2];
                gv_pca[2] = pca_mat[2][0]*gv_cam[i][0] + pca_mat[2][1]*gv_cam[i][1] + pca_mat[2][2]*gv_cam[i][2];

                // reverse-compensate refraction and camera perspective elongation of the true ellipse ratio
                float xyc[2] = { cxy[f][0]-img_ctr[0], cxy[f][1]-img_ctr[1] };
                float xye[2] = { cxy[f][0]-ctr[0],     cxy[f][1]-ctr[1] };
                float cra_proj = camRayProjToMinor(xyc, xye, iris_z);
                float el_r = cosf((acosf(gv_pca[2])+cra_proj)/1.121f)/cosf(cra_proj);
                phi = atan2f(gv_pca[1], gv_pca[0]) + F_PI/2;
                encl_rect(1, el_r, phi, &rw, &rh);
                // expected bounding rectangle ratio for given position
                rect_true_n = rect_true_n + rw/rh;
            }

            ratio[i] = rect_obs_n / rect_true_n;
        }

    // goal is to reach the same relation to true ratio at all positions
    float meanratio = 0;
    int nuse = 0;
    for (int i=0; i<N; ++i)
        if (use[i])
        {
            meanratio += ratio[i];
            ++nuse;
        }
    meanratio /= nuse;

    float diff = 0;
    for (int i=0; i<N; ++i)
        if (use[i])
            diff += (ratio[i]-meanratio)*(ratio[i]-meanratio);

    return diff;
}


float ctrDiffFromPupilRatio(float (*pxy)[MAX_PUPIL_XY][2], int *nxy, float (*cxy)[2], int anc_st[N], int anc_en[N], int use[N],
    float *el_ratio0, float *iris_ratio, float ctr[2], float img_ctr[2], float pxl2mm[N], float iris_z, float iris2eye_mm, float iris2hrot_mm, float iris2vrot_mm,
    float pup_ang, float pup_stretch, float ax, float ay, float pup_refr, float ctrn[N][2])
{
    float pca_angle[2] = {ax, ay};

    memset(ctrn, 0, N*2*sizeof(float));

    float diff = 0;
    for (int i=0; i<N; ++i)
    {
        if (use[i])
        {
            for (int f=anc_st[i]; f<anc_en[i]; ++f)
            {
                float xy[MAX_PUPIL_XY*2];

                // elongation will overwrite xy values, need to keep pxy for future passes
                memcpy(xy, pxy[f], nxy[f]*2*sizeof(float));
                // pupil elongation compensation
                compensatePupilElongation(xy, nxy[f], el_ratio0[f], pup_ang, pup_stretch);
                float A[6];
                if (!EllipseDirectFit(xy, nxy[f], 0, 0, A)) return 1e10f;
                float tmpcx, tmpcy, semx, semy, phi;
                EllipseGeomFromAlgebraic(A, &tmpcx, &tmpcy, &semx, &semy, &phi);
                float el_r = fminf(semx,semy) / fmaxf(semx,semy);

                float phi_minor = phi + (semx<semy)*F_PI/2;

                // apply refraction/perspective correction
                float xyc[2] = { cxy[f][0]-img_ctr[0], cxy[f][1]-img_ctr[1] };
                float xye[2] = { cxy[f][0]-ctr[0],     cxy[f][1]-ctr[1] };
                float cra_proj = camRayProjToMinor(xyc, xye, iris_z);
                el_r = cosf(acosf(el_r*cosf(cra_proj))*pup_refr - cra_proj);

                // pupillary circular axis angle compensation
                float vec_to_eyeb[2] = { ctr[0]-cxy[f][0], ctr[1]-cxy[f][1] };
                compensatePCA(vec_to_eyeb, iris_ratio[f], pca_angle, &el_r, &phi_minor);
                // phi_mior points toward 2d eyeball center, range: [-pi/2 .. 3*pi/2]

                // Note: unlike in Estimate2dEyeCenter(), estimation is in unrotated camera pixel coordinates here
                float dv[2] = { sinf(phi_minor), -cosf(phi_minor) };
                float proj_ang = acosf(el_r);
                float len2dy = iris2eye_mm / pxl2mm[i] * sinf(proj_ang);
                float len2dx = iris2hrot_mm / iris2vrot_mm * len2dy;
                float e1h = cxy[f][0] + len2dx*dv[0];
                float e1v = cxy[f][1] + len2dy*dv[1];

                ctrn[i][0] += e1h; ctrn[i][1] += e1v;
            }

            ctrn[i][0] /= anc_en[i]-anc_st[i];
            ctrn[i][1] /= anc_en[i]-anc_st[i];
            diff += sqrtf( norm_euc2(ctrn[i], ctr) );
        }
    }

    return diff;
}


int calibrateRotAxesPupilDecenter(eye_cfg *ei, uint8_t ** frames, uint8_t * frame, int M, float gv_cam[N][3], int primary_anc, float ctr[2], float *iris2pupil, float pxl2mm[N])
{
    int status = EYECAL_ALL_OK;

    ctr[0] = ctr[1] = 0;
    *iris2pupil = 0.3f;

    ei->smooth_2deye  = 1; // eyeball center detection jitter may cause center 'flipping' to erroneous cluster without smoothing
    ei->smooth_2dtail = 0;
    ei->smooth_3deye  = 0;
    ei->smooth_iris   = 3; // erroneous iris radius detecitons with low confidence will get in the way without iris smoothing

    // Note: in case of poor iris detection - pupil decenter calibraiton will fail but pupil decentration can be ignored
    // since even typically high decenter of 0.9mm will only result in
    // 0.9*sin(15) = 0.22mm systematic error in rotation axis distance in the very worst case

    memset(pxl2mm, 0, N*sizeof(float));

    float pupil2d[N][2] = {0};
    float iris2dx[N] = {0};
    float iris2dy[N] = {0};
    int   nirisx[N] = {0};
    int   nirisy[N] = {0};
    int   npupil[N] = {0};

    float *iris_r1[N]={NULL};
    float *iris_r2[N]={NULL};
    for (int i=0; i<N; ++i)
    {
        iris_r1[i] = (float*)malloc(M*sizeof(float));
        iris_r2[i] = (float*)malloc(M*sizeof(float));
        if ((iris_r1[i]==NULL) || (iris_r2[i]==NULL)) { status=EYECAL_NO_MEMORY; goto free_iris; }
    }
    int n_iris_r1[N]={0};
    int n_iris_r2[N]={0};

    ei->el_known = 1;
    ei->iris_recenter = 2; // re-estimate iris both horizontally and vertically (for pupil decenter estimation)

    for (int i=0; i<M; ++i)
    {
        int ianc = i*N/M;
        int use = (i+0.5f)*N/M - ianc >= 0.5f;
        // pre-fill iris ratio of previous frame with ground-truth value
        ei->iris_ratio = gv_cam[ianc][2];
        // pre-fill ground truth ratio and phi (will overwrite values detected from pupil)
        ei->el_ratio =  gv_cam[ianc][2];
        ei->phi_minor = atan2f(gv_cam[ianc][1], gv_cam[ianc][0])+F_PI/2;

        memcpy(frame, frames[i], ei->w*ei->h);
        getEyePose(ei, frame, NULL, NULL);

        if (use && (ei->blink < 0.5))
        {
            pxl2mm[ianc]   += ei->pxl2mm;
            pupil2d[ianc][0] += ei->iris2d_p[0];
            pupil2d[ianc][1] += ei->iris2d_p[1];
            ++npupil[ianc];

            // for decenter - use locations where iris2d was estimated from 2 radii
            if (ei->iris_r1 && ei->iris_r2)
            {
                iris2dx[ianc] += ei->iris2d[0];
                ++nirisx[ianc];

                iris2dy[ianc] += ei->iris2d[1];
                ++nirisy[ianc];
            }

            // record separate estimations of iris radius from left and right edge
            // to correct later on with pupil decenter and adjust resulting pxl2mm values
            if (ei->iris_r1>0) iris_r1[ianc][n_iris_r1[ianc]++] = ei->iris_r1;
            if (ei->iris_r2>0) iris_r2[ianc][n_iris_r2[ianc]++] = ei->iris_r2;
        }
    }

    // set back default values
    ei->el_known = 0;
    ei->iris_recenter = 1;

    float iris2d[N][2];
    for (int i=0; i<N; ++i)
    {
        if (npupil[i]<3) { status=EYECAL_NO_PUPIL; goto free_iris; }

        // averaged per-position pxl2mm, iris2d
        pxl2mm[i] /= npupil[i];
        pupil2d[i][0] /= npupil[i];
        pupil2d[i][1] /= npupil[i];
        if (nirisx[i]) iris2d[i][0] = iris2dx[i]/nirisx[i];
            else  iris2d[i][0] = pupil2d[i][0];
        if (nirisy[i]) iris2d[i][1] = iris2dy[i]/nirisy[i];
            else  iris2d[i][1] = pupil2d[i][1];
    }

    // averaged separate estimations of iris radius from left and right edge
    // excluding outliers
    float ir_r1[N] = {0};
    float ir_r2[N] = {0};
    int n_r1[N] = {0};
    int n_r2[N] = {0};

    for (int i=0; i<N; ++i)
    {
        if (n_iris_r1[i])
        {
            float mean_r1 = 0;
            for (int j=0; j<n_iris_r1[i]; ++j)
                mean_r1 += iris_r1[i][j];
            mean_r1 /= n_iris_r1[i];

            for (int j=0; j<n_iris_r1[i]; ++j)
                if (fabsf(iris_r1[i][j]-mean_r1) < 0.2f*mean_r1)
                {
                    ir_r1[i] += iris_r1[i][j];
                    ++n_r1[i];
                }
            ir_r1[i] /= n_r1[i];
        }

        if (n_iris_r2[i])
        {
            float mean_r2 = 0;
            for (int j=0; j<n_iris_r2[i]; ++j)
                mean_r2 += iris_r2[i][j];
            mean_r2 /= n_iris_r2[i];

            for (int j=0; j<n_iris_r2[i]; ++j)
                if (fabsf(iris_r2[i][j]-mean_r2) < 0.2f*mean_r2)
                {
                    ir_r2[i] += iris_r2[i][j];
                    ++n_r2[i];
                }
            ir_r2[i] /= n_r2[i];
        }
    }

    // estimate pxl2mm for all locations from location closest to the camera
    // using single base value of pxl2mm for better consistency
    float z_pri = ei->iris_z + ei->iris2eye_mm * (1-gv_cam[primary_anc][2]) / pxl2mm[primary_anc];
    for (int i=0; i<N; ++i)
    {
        if (i!=primary_anc)
        {
            float z = ei->iris_z + ei->iris2eye_mm * (1-gv_cam[i][2]) / pxl2mm[primary_anc];
            pxl2mm[i] = pxl2mm[primary_anc] * z/z_pri;
        }
    }

    // compute pupil decenter
    float pup_dec_h=0, pup_dec_v=0;
    int nh=0, nv=0;
    for (int i=0; i<N; ++i)
    {
        if (nirisx[i] > 2)
        {
            float sh = sqrtf( 1 - gv_cam[i][0]*gv_cam[i][0] );
            if (sh>0.5f)
            {
                pup_dec_h += (pupil2d[i][0]-iris2d[i][0]) / sh * pxl2mm[i];
                ++nh;
            }
        }

        if (nirisy[i] > 2)
        {
            float sv = sqrtf( 1 - gv_cam[i][1]*gv_cam[i][1] );
            if (sv>0.5)
            {
                pup_dec_v += (pupil2d[i][1]-iris2d[i][1]) / sv * pxl2mm[i];
                ++nv;
            }
        }
    }

    // pupil decenter is of low importance for off-axis camera, so it can be skipped
    if (nh) pup_dec_h /= nh;
    if (nv) pup_dec_v /= nv;

    float pxl2mm_obs[N]; // observed pxl2mm, used for 3d slippage correction
    memcpy(pxl2mm_obs, pxl2mm, N*sizeof(float));

    // re-estimate pxl2mm accounting for pupil decenter
    for (int i=0; i<N; ++i)
    {
        ir_r1[i] -= pup_dec_h/pxl2mm[i];
        ir_r2[i] += pup_dec_h/pxl2mm[i];
        if (n_r1[i]+n_r2[i] >= 3)
        {
            float iris_r = (ir_r1[i]*n_r1[i] + ir_r2[i]*n_r2[i]) / (n_r1[i]+n_r2[i]);
            pxl2mm_obs[i] = ei->hmn_iris_mm/2/iris_r;
        }
    }
    float pxl_scl = pxl2mm_obs[primary_anc]/pxl2mm[primary_anc];
    for (int i=0; i<N; ++i)
        pxl2mm[i] *= pxl_scl;

    // compute distance to axis (as gaze vector length) from the length between
    // reference vector projections and distances between iris2d's
    // assuming angle order: L R T B ...
    float d_ref_h = fabsf(gv_cam[0][0] - gv_cam[1][0]);
    float dist_mm_p_h = fabsf( pupil2d[0][0]*pxl2mm[0] -  pupil2d[1][0]*pxl2mm[1]) / d_ref_h;
    float dist_mm_h = dist_mm_p_h;
    float d_ref_v = fabsf(gv_cam[2][1] - gv_cam[3][1]);
    float dist_mm_p_v = fabsf( pupil2d[2][1]*pxl2mm[2] -  pupil2d[3][1]*pxl2mm[3]) / d_ref_v;
    float dist_mm_v = dist_mm_p_v;

    if (min(nirisx[0], nirisx[1]) >= 2)
        dist_mm_h = fabsf( iris2d[0][0]*pxl2mm[0] - iris2d[1][0]*pxl2mm[1]) / d_ref_h;

    if (min(nirisy[2], nirisy[3]) >= 2)
        dist_mm_v = fabsf( iris2d[2][1]*pxl2mm[2] - iris2d[3][1]*pxl2mm[3]) / d_ref_v;

    float iris2pupil_h = dist_mm_h - dist_mm_p_h + 0.3f;
    float iris2pupil_v = dist_mm_v - dist_mm_p_v + 0.3f;
    // horizontal iris ellipse estimation is more reliable
    // iris2pupil = ( iris2pupil_h + iris2pupil_v ) /2;
    // not expecting pupil-iris distance to be above 1mm
    *iris2pupil = fmaxf(-1, fminf(1, iris2pupil_h));

    // very small and very large values are possible with prescription glasses
    if ( (fminf(dist_mm_h,dist_mm_v)<4) || (fmaxf(dist_mm_h,dist_mm_v)>25) )  { status=EYECAL_AXES_RANGE; goto free_iris; }
    else
    {
        ei->iris2hrot_mm = dist_mm_h;
        ei->iris2vrot_mm = dist_mm_v;
        ei->iris2eye_mm  = ei->iris2vrot_mm;
        ei->pupil2eye_mm = ei->iris2vrot_mm - *iris2pupil;
        ei->pupil_decenter_mm[0] = pup_dec_h;
        ei->pupil_decenter_mm[1] = pup_dec_v;
    }

    // estimate eyeball center on image plane
    float gt[N][2];
    for (int i=0; i<N; ++i)
    {
        gt[i][0] = -gv_cam[i][0];
        gt[i][1] = -gv_cam[i][1];
    }
    float tmp[N][2];
    getCtrAndSlip(iris2d, gt, ctr, tmp);

free_iris:
    for (int i=0; i<N; ++i)
    {
        if (iris_r1[i]) free(iris_r1[i]);
        if (iris_r2[i]) free(iris_r2[i]);
    }

    return status;
}


int calibratePupilElongationPCA(eye_cfg *ei, uint8_t ** frames, uint8_t * frame, int M, float gv_cam[N][3], float ctr[2], float pxl2mm[N], float ctrn[N][2])
{
    int status = EYECAL_ALL_OK;

    memset(ctrn, 0, N*2*sizeof(float));

    ei->smooth_2deye  = -1;      // use eyeball center found in the first pass

    ei->eyeball2d[0] = ei->ctr_tail[0] = ei->eyeball_ctr[0] = ctr[0];
    ei->eyeball2d[1] = ei->ctr_tail[1] = ei->eyeball_ctr[1] = ctr[1];

    float (*pxy)[MAX_PUPIL_XY][2] = (float (*)[MAX_PUPIL_XY][2])malloc(M*MAX_PUPIL_XY*2*sizeof(float));
    int   *nxy = (int *)malloc(M*sizeof(int));
    float (*cxy)[2] = (float (*)[2])malloc(M*2*sizeof(float));
    float *el_ratio0 = (float *)malloc(M*sizeof(float));
    float *iris_ratio = (float *)malloc(M*sizeof(float));

    if ( (pxy == NULL) ||  (nxy == NULL) || (cxy == NULL)  || (el_ratio0 == NULL)  || (iris_ratio == NULL) )
    {
        status = EYECAL_NO_MEMORY; goto free_pass2;
    }

    int numpos[N] = {0};
    for (int i=0, j=0; i<M; ++i)
    {
        int ianc = i*N/M;
        int use = (i+0.5f)*N/M - ianc >= 0.5f;

        memcpy(frame, frames[i], ei->w*ei->h);
        getEyePose(ei, frame, NULL, NULL);

        if ( use && (ei->blink < 0.5) )
        {
            nxy[j] = min(ei->xylen, MAX_PUPIL_XY);
            memcpy(pxy[j], ei->xy, nxy[j]*2*sizeof(float));
            cxy[j][0] = ei->cxy[0];
            cxy[j][1] = ei->cxy[1];
            el_ratio0[j] = ei->el_ratio0;
            iris_ratio[j] = ei->iris_ratio;
            
            ++numpos[ianc];
            ++j;
        }
    }

    for (int i=0; i<N; ++i)
        if (((i<=3) || (i==N)) && (numpos[i]<3)) { status = EYECAL_NO_PUPIL; goto free_pass2; }

    // start and end of particular calibration location in arrays
    int anc_st[N], anc_en[N];
    anc_st[0] = 0;
    anc_en[0] = numpos[0];
    for (int i=1; i<N; ++i)
    {
        anc_st[i] = anc_en[i-1];
        anc_en[i] = anc_en[i-1] + numpos[i];
    }

    // compute el_ratio/phi_minor corrections (pca, elongation) from reference given by gaze vectors in camera space vs measured

    // there might be slight differences between eyeball centers found in this pass and the next one
    // since iris centers are used here for camRayProjToMinor() and compensatePCA()
    // and pupil centers are used in the next pass (and in normal tracking)

    // find the best pup_ang and pup_stretch via iterative minimization
    float best_pup_refr = ei->pupil_refr;
    float best_pup_ang = 0;
    float best_pup_stretch = 1.0f;
    float best_ax = 0;
    float best_ay = 0;

    // first run - find pupilaary circular axis (pca)
    float img_ctr[2] = {(float)ei->w/2, (float)ei->h/2};
    float ctr_unrot[2];
    rotcam(ctr, -ei->camera_rot, ctr_unrot);
    ctr_unrot[0] += img_ctr[0];
    ctr_unrot[1] += img_ctr[1];

    // only calibrate for pca if positions above 9 degree from both sides relative to the eye center are present
    float pcax_st=0, pcax_en=0, pcax_step=0.25;
    float pcay_st=0, pcay_en=0, pcay_step=0.25;
    float gv_minx=gv_cam[0][0], gv_miny=gv_cam[0][1];
    float gv_maxx=gv_cam[0][0], gv_maxy=gv_cam[0][1];
    for (int i=1; i<N; ++i)
    {
        gv_minx = fminf(gv_minx, gv_cam[i][0]);
        gv_miny = fminf(gv_miny, gv_cam[i][1]);
        gv_maxx = fmaxf(gv_maxx, gv_cam[i][0]);
        gv_maxy = fmaxf(gv_maxy, gv_cam[i][1]);
    }
    if ( (gv_minx < -0.15f) && (gv_maxx > 0.15f) )
    {
        pcax_st = -15; pcax_en = 15;
    }
    if ( (gv_miny < -0.15f) && (gv_maxy > 0.15f) )
    {
        pcay_st = -15; pcay_en = 15;
    }

    // assuming angle order: L R T B ... C

    float best_diff = 1e10f;
    for (float ax=pcax_st; ax<=pcax_en; ax+=pcax_step)
    {
        int use[N] = {1, 1, 0, 0, DIAGS_FILL0 1};
        float diff = diffFromPCA(pxy, nxy, cxy, anc_st, anc_en, gv_cam, use, ctr_unrot, img_ctr, ei->iris_z, ax, 0);
        if (diff < best_diff)
        {
            best_ax = ax;
            best_diff = diff;
        }
    }

    best_diff = 1e10f;
    for (float ay=pcay_st; ay<=pcay_en; ay+=pcay_step)
    {
        int use[N] = {0, 0, 1, 1, DIAGS_FILL0 1};
        float diff = diffFromPCA(pxy, nxy, cxy, anc_st, anc_en, gv_cam, use, ctr_unrot, img_ctr, ei->iris_z, best_ax, ay);
        if (diff < best_diff)
        {
            best_ay = ay;
            best_diff = diff;
        }
    }

    // re-run horizontal pass with best vertical pca estimation
    best_diff = 1e10f;
    for (float ax=pcax_st; ax<=pcax_en; ax+=pcax_step)
    {
        int use[N] = {1, 1, 0, 0, DIAGS_FILL0 1};
        float diff = diffFromPCA(pxy, nxy, cxy, anc_st, anc_en, gv_cam, use, ctr_unrot, img_ctr, ei->iris_z, ax, best_ay);
        if (diff < best_diff)
        {
            best_ax = ax;
            best_diff = diff;
        }
    }

    ei->pca_angle[0] = best_ax;
    ei->pca_angle[1] = best_ay;

    // second run - find optimal parameters to best match direction towards center (but not scale) via pupil elongation:
    // - ignore camera perspective projection  ... cra-proj has heavy influence on angles ??? ... perform same correction as to gv_cam in first run ... ???
    // - ignore pupil refraction
    // - parameters affecting: pup_ang, pup_stretch

    // only calibrate for pupil elongation if sufficienctly front-looking position is present
    float pup_ang_st1=0, pup_ang_en1=0, pup_ang_step1=F_PI/8;
    float pup_ang_st2=0, pup_ang_en2=0, pup_ang_step2=F_PI/32;
    float pup_stretch_st1=1, pup_stretch_en1=1, pup_stretch_step1=0.05f;
    float pup_stretch_st2=0, pup_stretch_en2=0, pup_stretch_step2=0.005f;
    for (int i=0; i<N; ++i)
        if (gv_cam[i][2]>= 0.98f)
        {
            pup_ang_st1=0;          pup_ang_en1=F_PI;
            pup_ang_st2=-F_PI/8;    pup_ang_en2=F_PI/8;
            pup_stretch_st1=0.8f;   pup_stretch_en1=1.0f;
            pup_stretch_st2=-0.05f; pup_stretch_en2=0.05f;
            break;
        }

    best_diff = 1e10f;
    int use[N] = {1, 1, 1, 1, DIAGS_FILL1 0};
    float tmp[N][2];
    float best1_pup_ang = 0;
    float best1_pup_stretch = 1;
    for (float pup_ang=pup_ang_st1; pup_ang<=pup_ang_en1; pup_ang+=pup_ang_step1)
        for (float pup_stretch=pup_stretch_st1; pup_stretch<=pup_stretch_en1; pup_stretch+=pup_stretch_step1)
        {
            float diff = ctrDiffFromPupilRatio(pxy, nxy, cxy, anc_st, anc_en, use, el_ratio0, iris_ratio, ctr_unrot,
                img_ctr, pxl2mm, ei->iris_z, ei->iris2eye_mm, ei->iris2hrot_mm, ei->iris2vrot_mm,
                pup_ang, pup_stretch, best_ax, best_ay, best_pup_refr, tmp);
            if (diff<best_diff)
            {
                best1_pup_ang = pup_ang;
                best1_pup_stretch = pup_stretch;
                best_diff = diff;
            }
        }

    best_diff = 1e10f;
    best_pup_ang = best1_pup_ang;
    best_pup_stretch = best1_pup_stretch;
    for (float pup_ang=best1_pup_ang+pup_ang_st2; pup_ang<=best1_pup_ang+pup_ang_en2; pup_ang+=pup_ang_step2)
        for (float pup_stretch=best1_pup_stretch+pup_stretch_st2; pup_stretch<=best1_pup_stretch+pup_stretch_en2; pup_stretch+=pup_stretch_step2)
        {
            float diff = ctrDiffFromPupilRatio(pxy, nxy, cxy, anc_st, anc_en, use, el_ratio0, iris_ratio, ctr_unrot,
                img_ctr, pxl2mm, ei->iris_z, ei->iris2eye_mm, ei->iris2hrot_mm, ei->iris2vrot_mm,
                pup_ang, pup_stretch, best_ax, best_ay, best_pup_refr, tmp);
            if (diff<best_diff)
            {
                best_pup_ang = pup_ang;
                best_pup_stretch = pup_stretch;
                best_diff = diff;
            }
        }

    ctrDiffFromPupilRatio(pxy, nxy, cxy, anc_st, anc_en, use, el_ratio0, iris_ratio, ctr_unrot,
        img_ctr, pxl2mm, ei->iris_z, ei->iris2eye_mm, ei->iris2hrot_mm, ei->iris2vrot_mm,
        best_pup_ang, best_pup_stretch, best_ax, best_ay, best_pup_refr, ctrn);

    for (int i=0; i<N; ++i)
    {
        ctrn[i][0] -= img_ctr[0];
        ctrn[i][1] -= img_ctr[1];
        rotcam(ctrn[i], ei->camera_rot, ctrn[i]);
    }

    ei->pupil_refr = best_pup_refr;
    ei->pupil_angle = best_pup_ang*180/F_PI;
    ei->pupil_stretch = best_pup_stretch;

free_pass2:
    if (pxy) free(pxy);
    if (nxy) free(nxy);
    if (cxy) free(cxy);
    if (el_ratio0) free(el_ratio0);
    if (iris_ratio) free(iris_ratio);

    return status;
}


void calibrateAlphaAngle(eye_cfg *ei, uint8_t ** frames, uint8_t * frame, int M, float anc[N][2], float ctr[2], float pxl2mm[N], float alpha[2])
{
    // Note: alpha angle is due to the angle between optical and visual eye axis
    //       alpha stretching (scale) - either due to the pupil swim effects of HMD or eyewear

    // using this to start center estimations afresh (removes trash from the first pass)
    memset(ei->valid, 0, N_EYE_OBSERVATIONS*sizeof(int));
    ei->smooth_2deye  = 1;
    ei->smooth_2dtail = 4;
    ei->smooth_3deye  = 0;
    ei->smooth_iris   = 3;
    ei->eyeball2d[0] = ei->ctr_tail[0] = ei->eyeball_ctr[0] = ctr[0];
    ei->eyeball2d[1] = ei->ctr_tail[1] = ei->eyeball_ctr[1] = ctr[1];

    float rotrad      = ei->camera_rot * F_PI/180;

    alpha[0] = alpha[1] = 0;
    int npos = 0;
    for (int i=0; i<M; ++i)
    {
        int ianc = i*N/M;
        int use = (i+0.5f)*N/M - ianc >= 0.5f;

        memcpy(frame, frames[i], ei->w*ei->h);
        getEyePose(ei, frame, NULL, NULL);

        if ( use && (ei->blink < 0.5) )
        {
            float gaze_deg[2] = { 180/F_PI * asinf(-ei->gaze_vector[0]), 180/F_PI * asinf(-ei->gaze_vector[1]) };

            alpha[0] += anc[ianc][0] - gaze_deg[0];
            alpha[1] += anc[ianc][1] - gaze_deg[1];
            ++npos;
        }
    }

    if (npos)
    {
        alpha[0] /= npos;
        alpha[1] /= npos;
    }

    float alpha3[3] = {alpha[1], -alpha[0], 0};
    rotationMat(alpha3, ei->alpha_mat, 0);
}


// Obtain eye calibration parameters from a sequence of captured frames
//
// anc - calibration locations in degree from center (Left Right Top Bottom Center)
//       recommended anc = [-15 0; 15 0; 0 -15; 0 15; 0 0];
// M   - number of input frames
//
// ToDo: return status at which step the failure happened
//
eye_calib * calibrateEye(eye_cfg *ei, uint8_t ** frames, int M, float anc[N][2], int *err)
{
    eye_calib *ec = NULL;

    // frames are modified during getEyePose(), but multiple calibration passes should work on the same data
    // so, create separate storage where data from frames[] is copied and processed
    uint8_t *frame = malloc(ei->w*ei->h);
    *err = EYECAL_NO_MEMORY;
    if (frame == NULL) return NULL;

    // assume equal number of frames per calibration location
    // and use only second half of each frame interval for measurements

    // calculate projection of reference gaze vectors onto camera image plane
    float gv_cam[N][3];
    for (int i=0; i<N; ++i)
    {
        // reference gaze vector in HMD space
        float gv[3] = { -sinf(anc[i][0]*F_PI/180), -sinf(anc[i][1]*F_PI/180), 0 };
        gv[2] = sqrtf( 1-(gv[0]*gv[0]+gv[1]*gv[1]) );

        gv_cam[i][0] = ei->hmd2cam[0][0]*gv[0] + ei->hmd2cam[0][1]*gv[1] + ei->hmd2cam[0][2]*gv[2];
        gv_cam[i][1] = ei->hmd2cam[1][0]*gv[0] + ei->hmd2cam[1][1]*gv[1] + ei->hmd2cam[1][2]*gv[2];
        gv_cam[i][2] = ei->hmd2cam[2][0]*gv[0] + ei->hmd2cam[2][1]*gv[1] + ei->hmd2cam[2][2]*gv[2];
    }

    // location closest to the camera (most direct look to the camera)
    float closest = gv_cam[0][2];
    int primary_anc = 0;
    for (int i=1; i<N; ++i)
        if (gv_cam[i][2] > closest)
        {
            primary_anc = i;
            closest = gv_cam[i][2];
        }

    // reset eye instance calibrations
    memset(ei->alpha_mat, 0, 3*3*sizeof(float));
    ei->alpha_mat[0][0] = ei->alpha_mat[1][1] = ei->alpha_mat[2][2] = 1;

    // -------- first pass - initial estimation of position of eyeball rotation axes and pupil decentration

    float iris2pupil;
    float pxl2mm[N];
    float ctr[2];
    *err = calibrateRotAxesPupilDecenter(ei, frames, frame, M, gv_cam, primary_anc, ctr, &iris2pupil, pxl2mm);
    if (*err!=EYECAL_ALL_OK) goto free_frame; // Rotation axes calibration fail

    // -------- second pass - estimation of pupil elongation and pupillary circular axis angles

    float ctrn[N][2];
    *err = calibratePupilElongationPCA(ei, frames, frame, M, gv_cam, ctr, pxl2mm, ctrn);
    if (*err!=EYECAL_ALL_OK) goto free_frame; // Pupil calibration fail

    // -------- third pass - figure alpha angle and rotation axes based on ellipse ratio

    float alpha[2];
    calibrateAlphaAngle(ei, frames, frame, M, anc, ctr, pxl2mm, alpha);


    ec = (eye_calib *)malloc(sizeof(eye_calib));
    if (ec==NULL) { *err = EYECAL_NO_MEMORY; goto free_frame; }
    memset(ec, 0, sizeof(eye_calib));

    ec->iris2hrot         = ei->iris2hrot_mm;
    ec->iris2vrot         = ei->iris2vrot_mm;
    ec->iris2pupil        = iris2pupil;
    ec->pupil_decenter[0] = ei->pupil_decenter_mm[0];
    ec->pupil_decenter[1] = ei->pupil_decenter_mm[1];
    ec->eyeball_ctr[0]    = ei->eyeball_ctr[0];
    ec->eyeball_ctr[1]    = ei->eyeball_ctr[1];
    ec->pca_angle[0]      = ei->pca_angle[0];
    ec->pca_angle[1]      = ei->pca_angle[1];
    ec->pupil_refr        = ei->pupil_refr;
    ec->pupil_angle       = ei->pupil_angle;
    ec->pupil_stretch     = ei->pupil_stretch;
    ec->gaze_alpha[0]     = alpha[0];
    ec->gaze_alpha[1]     = alpha[1];

free_frame:
    free(frame);
    return ec;
}


void releaseCalibration(eye_calib * ec)
{
    if (ec) free(ec);
}
