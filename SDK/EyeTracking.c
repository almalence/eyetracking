/* ------------------------------------------------------------------------- *\

    Estimate eye pose and gaze vector from camera frame
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <memory.h>

#include "EyeUtil.h"
#include "EyeTracking.h"
#include "ellipse.h"

#ifdef HAVE_NEON
#include <arm_neon.h>
#endif

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif


// macro to access pixel value of monochrome camera frame
#define PX(im, x, y)    (im[(y)*w+(x)])
#define PXI(im, x, y)   ((int)(im[(y)*w+(x)]))
// same for 2x downscaled (=half width) image
#define PX2(im, x, y)    (im[(y)*(w/2)+(x)])

// angle from horizontal up to which iris is usually visible (15 degree up, 30 degree down)
#define IRIS_SCTR_UP     (F_PI/12)
#define IRIS_SCTR_DN     (F_PI/6)
#define IRIS_SCTR        (IRIS_SCTR_UP+IRIS_SCTR_DN)
#define IRIM_MIDH        ((int)(IRIM_H*IRIS_SCTR_UP/IRIS_SCTR))    // 26 iris scan middle line


// Squared Eucleadean norm of 2d vectors
float norm_euc2(float *v1, float *v2)
{
    return (v1[0]-v2[0]) * (v1[0]-v2[0]) + (v1[1]-v2[1]) * (v1[1]-v2[1]);
}

// soncos function from mingwex library. put here to make MSVC linker happy
#ifdef WIN32
void sincosf (float __x, float *p_sin, float *p_cos)
{
  long double c, s;

  __asm__ __volatile__ ("fsincos\n\t"
    "fnstsw    %%ax\n\t"
    "testl     $0x400, %%eax\n\t"
    "jz        1f\n\t"
    "fldpi\n\t"
    "fadd      %%st(0)\n\t"
    "fxch      %%st(1)\n\t"
    "2: fprem1\n\t"
    "fnstsw    %%ax\n\t"
    "testl     $0x400, %%eax\n\t"
    "jnz       2b\n\t"
    "fstp      %%st(1)\n\t"
    "fsincos\n\t"
    "1:" : "=t" (c), "=u" (s) : "0" (__x) : "eax");
  *p_sin = (float) s;
  *p_cos = (float) c;
}
#else
void sincosf (float __x, float *p_sin, float *p_cos);   // somehow didn't found this in math.h with gcc 16 under linux, but it is present in libm
#endif

// using volatile here to prevent sinf/cosf optimization into sincosf, which is problematic with MSVC linkage
void rotcam(float crd[2], float ang, float out[2])
{
    ang *= F_PI/180;

    // temporary storage to prevent overwrite if crd and out are the same
    float tmp[2];
    float sang,cang; sincosf(ang, &sang, &cang);
    tmp[0] =  crd[0]*cang - crd[1]*sang;
    tmp[1] =  crd[0]*sang + crd[1]*cang;
    out[0] = tmp[0]; out[1] = tmp[1];
}

// Rotation matrix from rotation angles in degree
// X-Y-Z Euler rotation, right-hand rule
void rotationMat(float *angles, float mat[3][3], int inv)
{
    float rx = angles[0]*F_PI/180;
    float ry = angles[1]*F_PI/180;
    float rz = angles[2]*F_PI/180;

    if (inv)
    {
        rx=-rx; ry=-ry; rz=-rz;
    }

    float srx, crx; sincosf(rx, &srx, &crx);
    float sry, cry; sincosf(ry, &sry, &cry);
    float srz, crz; sincosf(rz, &srz, &crz);

    float RX[3][3] = {{    1,    0,    0 },
                      {    0,  crx, -srx },
                      {    0,  srx,  crx }};
    float RY[3][3] = {{  cry,    0,  sry },
                      {    0,    1,    0 },
                      { -sry,    0,  cry }};
    float RZ[3][3] = {{  crz, -srz,    0 },
                      {  srz,  crz,    0 },
                      {    0,    0,    1 }};
    float RXY[3][3];

    if (inv)
    {
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                RXY[i][j] = RY[i][0]*RZ[0][j] + RY[i][1]*RZ[1][j] + RY[i][2]*RZ[2][j];
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                mat[i][j] = RX[i][0]*RXY[0][j] + RX[i][1]*RXY[1][j] + RX[i][2]*RXY[2][j];
    }
    else
    {
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                RXY[i][j] = RY[i][0]*RX[0][j] + RY[i][1]*RX[1][j] + RY[i][2]*RX[2][j];
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                mat[i][j] = RZ[i][0]*RXY[0][j] + RZ[i][1]*RXY[1][j] + RZ[i][2]*RXY[2][j];
    }
}

// pre-computed iris raster angular values: IRIS_SCTR_UP-IRIS_SCTR*(0:IRIM_H-1)/(IRIM_H-1)
float iris_d[IRIM_H] = {
    0.26180458f, 0.25155067f, 0.24174213f, 0.23196411f, 0.22200012f, 0.21202278f, 0.20225430f, 0.19224548f, 0.18222523f, 0.17222786f,
    0.16246414f, 0.15245914f, 0.14236450f, 0.13243008f, 0.12260723f, 0.11268616f, 0.10269213f, 0.09282494f, 0.08281565f, 0.07292223f,
    0.06291771f, 0.05300260f, 0.04305601f, 0.03310466f, 0.02320540f, 0.01324654f, 0.00331225f, -0.00662562f, -0.01655233f, -0.02649311f,
    -0.03644753f, -0.04641035f, -0.05634525f, -0.06617737f, -0.07620430f, -0.08611155f, -0.09611943f, -0.10601237f, -0.11602116f, -0.12584782f,
    -0.13577271f, -0.14583111f, -0.15561676f, -0.16563129f, -0.17564869f, -0.18565178f, -0.19541454f, -0.20536613f, -0.21536255f, -0.22540665f,
    -0.23521328f, -0.24514771f, -0.25503159f, -0.26486015f, -0.27512360f, -0.28472900f, -0.29508018f, -0.30490875f, -0.31464577f, -0.32489204f,
    -0.33467865f, -0.34432602f, -0.35455322f, -0.36434364f, -0.37410927f, -0.38436317f, -0.39404297f, -0.40429688f, -0.41406250f, -0.42382813f,
    -0.43408203f, -0.44384766f, -0.45361328f, -0.46386719f, -0.47363281f, -0.48339844f, -0.49365234f, -0.50292969f, -0.51269531f, -0.52343750f
};

eye_cfg * initEyeConfig(int w, int h, int min_iris_r, int max_iris_r, int valid_area[6], int pup_thr,
                         float camera_hfov, float camera_rot, int fps, float cam2hmd[6])
{
    // -------- allocate eye tracking instance

    eye_cfg *ei = (eye_cfg *)malloc(sizeof(eye_cfg));
    if (ei==NULL) return NULL;
    memset(ei, 0, sizeof(eye_cfg));

    int irim_w = (int)(max_iris_r-min_iris_r+1);

    // +32 bytes to keep vectorized avx / neon code happy
    ei->tmp = (uint8_t*)malloc(w*h+32);
    ei->im2x = (uint8_t*)malloc((w/2)*(h/2)+32);
    ei->digest = (int16_t*)malloc((w/4)*(h/4)*sizeof(int16_t));
    ei->dweight = (int16_t*)malloc((w/4)*(h/4)*sizeof(int16_t));
    ei->iris_raster = (uint8_t*)malloc(IRIM_H*irim_w*sizeof(uint8_t));
    if ( (ei->tmp==NULL) || (ei->im2x==NULL) || (ei->digest==NULL) || (ei->dweight==NULL) || (ei->iris_raster==NULL) )
    {
        releaseEyeConfig(ei);
        return NULL;
    }
    ei->xy = (float*)ei->digest;    // use same storage for flating-point [x,y] pupil edge coordinates

    // -------- Camera-related constants

    ei->fps = fps;
    ei->w = w;
    ei->h = h;
    ei->irim_w = irim_w;
    ei->min_iris_r  = min_iris_r;                                     // minimum and maximum iris radius, in pixels
    ei->max_iris_r  = max_iris_r;
    ei->iris_z = (w/2)/tanf(camera_hfov/2*F_PI/180);                  // z-axis distance between camera and pupil (or image plane) in pixels

    ei->fps30       = (int)ceilf(fps/30.0f);                          // how many frames per 1/30th of a second
    ei->fps_ec2d    = (int)(fps/ei->fps30+0.5f);                      // fps at which 2d eyeball center estimations are performed

    ei->pup_thr     = pup_thr;
    ei->camera_rot  = camera_rot;                                     // how much to rotate (clockwise) camera frame to make eye nearly-horizontal
    memcpy(ei->valid_area, valid_area, sizeof(int[6]));

    // Camera-to-HMD coordinate space: x/y/z rotations (degree), x/y/z displacement (mm)
    // X-Y-Z Euler rotations used (default in eg Blender) 
    // Note: camera_zrot above is excluded from Z rotation below
    float tmpmat[3][3];
    rotationMat(cam2hmd, tmpmat, 0);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            ei->cam2hmd[i][j] = tmpmat[i][j];
    ei->cam2hmd[0][3] = cam2hmd[3];
    ei->cam2hmd[1][3] = cam2hmd[4];
    ei->cam2hmd[2][3] = cam2hmd[5];

    // Inverse (HMD to Camera) rotation matrix, used in calibration
    rotationMat(cam2hmd, tmpmat, 1);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            ei->hmd2cam[i][j] = tmpmat[i][j];
    ei->hmd2cam[0][3] = -cam2hmd[3];
    ei->hmd2cam[1][3] = -cam2hmd[4];
    ei->hmd2cam[2][3] = -cam2hmd[5];


    // --------  Default temporal filtering strength
    // - zero means no adjustment or filtering
    // - the higher the number - the stronger the smoothing (proportional to 2^n)
    // - the higher the number - the weaker the adjustment (proportional to 2^-n)
    int f = (int)(log2f(ei->fps_ec2d/30.f));
    ei->smooth_2deye  = max(1, 3+f);    // 2d eyeball coordinate
    ei->smooth_2dtail = max(2, 4+f);    // 2d eyeball coordinate tail (about 1-2 secs)
    f = (int)(log2f(fps/30.f));
    ei->smooth_3deye  = max(2, 3+f);    // 3d eyeball coordinate
    ei->smooth_iris   = max(1, 5+f);    // iris radius in pixels

    ei->gv_smooth_diff = 0.035f;        // apply smoothing to gaze vector if readings are within ~2 degree

    // -------- Anatomic constants (may be adjusted in calibration)

    // Iris to rotation axes
    // 15mm and 12mm from cornea apex; minus roughly 3mm from cornea apex to plane where iris edge is
    // A.Ohlendorf etc. 2022. Positions of the horizontal and vertical centre of rotation in eyes with different refractive errors    
    // Horizontal rotation has axis farther from the iris than the true eyeball center
    // Vertical rotation has axis at about the true eyeball center
    // D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 10
    ei->iris2hrot_mm = 11.5f;
    ei->iris2vrot_mm = 9.4f;

    // Iris to eyeball center
    // 7.92 mm 
    // Model of the human eye developed for MCNPX (model A) based on the dimensions provided in
    // NCRP Report no 130 (NCRP 1999) and the eye model from Charles and Brown (1975)
    // https://www.researchgate.net/figure/Model-of-the-human-eye-developed-for-MCNPX-model-A-based-on-the-dimensions-provided-in_fig2_51703220
    // Estimate2dEyeCenter and FilterEyeballPosition assume that eyeball center coincide with vertical rotation center
    ei->iris2eye_mm = ei->iris2vrot_mm;

    // Pupil to eyeball center, more or less the same as iris-to-eye, but:
    // - moves back and forth with accommodation by up to 0.4mm (accommodated lens push entrance pupil forward)
    // - can be different in synthetic eyes
    ei->pupil2eye_mm = ei->iris2eye_mm - 0.3f;

    // Matrix for angle alpha (between optical and visual axes)
    // typical values: 4 degree nasal, 3 degree upwards
    // D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 64
    // horizontal angle, then vertical; for the left eye coordinate system is starting at top-left corner (looking at HMD screen)
    float gaze_alpha3[3] = {4, 2, 0};
    rotationMat(gaze_alpha3, ei->alpha_mat, 0);

    // Iris diameter
    // 12 mm
    // H. Gross. 2008. Handbook of Optical Systems: Vol. 4 Survey of Optical Instruments.
    // Wiley-VCH Verlag GmbH+Co. KGaA, Weinheim.
    //
    // 11.6mm
    // Guillaume Francois, Pascal Gautron, Gaspard Breton, and Kadi Bouatouch
    // Image-Based Modeling of the Human Eye
    // IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 15, NO. 5, SEPTEMBER/OCTOBER 2009
    ei->hmn_iris_mm = 11.6f;

    ei->pupil_refr = 1.121f;                                    // pupil widening due to refraction
    ei->pupil_stretch = 1;                                      // pupil elongation (non-circularity)

    // rate at which HMD slippage (eyeball center shift) may occur [0..1]
    // may need to be set higher for synthetic datasets where center shifts abruptly
    ei->slippage_rate = 0.3f;

    // -------- state data

    // most fields are already pre-inited to 0 with memset above

    ei->gaze_vector[2] = 1.0f;                                  // forward center for initial gaze vector
    ei->gaze_vector_smooth[2] = 1.0f;
    ei->blink = 1;                                              // assume fully-closed eye initially
    ei->pxl2mm = ei->hmn_iris_mm/(min_iris_r+max_iris_r);
    ei->iris_ratio = MAX_EL_RATIO;                              // stabilized ratio of iris ellipse

    ei->iris_recenter = 1;                                      // re-estimate iris center from iris edges: 0=off 1=horizontal only 2=both 3=both + ellipse angle(phi_minor)

    ei->eyeball2d[0] = ei->eyeball_ctr[0];                      // assuming initial eyeball 2d position at designed position
    ei->eyeball2d[1] = ei->eyeball_ctr[1];

    ei->iris3d_cam[2] = ei->iris_z*ei->pxl2mm;                  // estimated 3d iris and eye center locations in camera space
    ei->eyeball3d_cam[2] = ei->iris3d_cam[2] + ei->iris2eye_mm;
    ei->eyeball3d0x[2] = ei->eyeball3d_cam[2];                  // stabilized eye center in camera space, assuming 0-degree vertical rotation
    ei->iris3d[2] = ei->iris2eye_mm;                            // estimated 3d iris center location in hmd space

    // values used during calibration
    ei->cxy[0] = w/2.0f;                                        // pupil center in camera frame
    ei->cxy[1] = h/2.0f;
    ei->el_ratio0 = MAX_EL_RATIO;                               // pupil ellipse ratio before any corrections
    ei->el_ratio = MAX_EL_RATIO;                                // pupil ellipse ratio (after refraction/projection/etc corrections)
    
    // -------- state data

    ei->ctr_tail[0] = ei->eyeball_ctr[0];
    ei->ctr_tail[1] = ei->eyeball_ctr[1];
    ei->last_bin = -1;
    
    // -------- prefill iris scan direction angles

    // look to the left / right of the pupil center, at up to [-30 30] degree from horizontal line,
    float rot = camera_rot * F_PI/180;
    for (int i=0; i<IRIM_H; ++i)
    {
        ei->iris_d1[i] = F_PI + iris_d[i] - rot;
        ei->iris_d2[i] =      - iris_d[i] - rot;
    }

    return ei;
}


void releaseEyeConfig(eye_cfg *ei)
{
    if (ei)
    {
        if (ei->tmp) free(ei->tmp);
        if (ei->im2x) free(ei->im2x);
        if (ei->digest) free(ei->digest);
        if (ei->dweight) free(ei->dweight);
        if (ei->iris_raster) free(ei->iris_raster);
        free(ei);
    }
}


// Default values:
// iris2hrot       11
// iris2vrot       9
// iris2pupil      0.3
// pupil_decenter  [0 0]
// eyeball_ctr     [0 0]
// pca_angle       [0 0]
// pupil_refr      1.121
// pupil_angle     0
// pupil_stretch   1
// gaze_alpha      [4 -2]
void applyCalibration(eye_cfg *ei, eye_calib *ec)
{

    ei->iris2hrot_mm         = ec->iris2hrot;
    ei->iris2vrot_mm         = ec->iris2vrot;
    ei->iris2eye_mm          = ei->iris2vrot_mm;
    ei->pupil2eye_mm         = ei->iris2vrot_mm - ec->iris2pupil;

    ei->pupil_decenter_mm[0] = ec->pupil_decenter[0];
    ei->pupil_decenter_mm[1] = ec->pupil_decenter[1];

    ei->eyeball2d[0] = ei->ctr_tail[0] = ei->eyeball_ctr[0]  = ec->eyeball_ctr[0];
    ei->eyeball2d[1] = ei->ctr_tail[1] = ei->eyeball_ctr[1]  = ec->eyeball_ctr[1];

    ei->pca_angle[0]         = ec->pca_angle[0];
    ei->pca_angle[1]         = ec->pca_angle[1];
    ei->pupil_refr           = ec->pupil_refr;
    ei->pupil_angle          = ec->pupil_angle;
    ei->pupil_stretch        = ec->pupil_stretch;

    float alpha3[3] = {ec->gaze_alpha[1], -ec->gaze_alpha[0], 0};
    rotationMat(alpha3, ei->alpha_mat, 0);
}


#define HORZ_BODY(x)                                          \
    uint8_t pmax = max(lt, rt);                               \
    uint8_t pmin = min(lt, rt);                               \
    int pxi = min(pmax, max(pmin, PX(tmp, (x), y)));          \
    /* image brightening towards corners */                   \
    int xc = (x)-valid[0];                                    \
    int yc = y-valid[1];                                      \
    pxi += ( (((pxi * (xc*xc+yc*yc))>>16) * brc + 128) >>8);  \
    pxi = min(255, pxi);

#ifdef HAVE_NEON
void removeSpeckles1_neon(uint8_t *restrict im, uint8_t *restrict tmp, uint8_t *restrict im2x, int valid[6], int w, int h)
{
    // vertical pass
    memcpy(&PX(tmp, 0, 0), &PX(im, 0, 0), w);
    memcpy(&PX(tmp, 0, 2), &PX(im, 0, 2), w);

    for(int y = 4; y < h-4; y += 2)
    {
        uint8x16_t v1, v2, v3, v4, v5, v6, v7, v8, v9;
        uint8_t *p = &PX(im, 0, y-4);
        uint8_t *po = &PX(tmp,0,y);

        for(int x=0; x<w; x += 16)
        {
            v1 = vld1q_u8(p); v2 = vld1q_u8(p+w); v3 = vld1q_u8(p+2*w); v4 = vld1q_u8(p+3*w); v5 = vld1q_u8(p+4*w);
            v6 = vld1q_u8(p+5*w); v7 = vld1q_u8(p+6*w); v8 = vld1q_u8(p+7*w); v9 = vld1q_u8(p+8*w); p += 16;
            // ToDo: may need +1 (and in horizontal pass too) to match AVX version (AVX computes (a+b+1)/2)
            v1 = vrhaddq_u8(vrhaddq_u8(v1, v2), vrhaddq_u8(v3, v4));
            v2 = vrhaddq_u8(vrhaddq_u8(v6, v7), vrhaddq_u8(v8, v9));
            v1 = vminq_u8(vmaxq_u8(v1, v2), vmaxq_u8(vminq_u8(v1, v2), v5));
            vst1q_u8(po, v1); po += 16;
        }
    }

    memcpy(&PX(tmp, 0, h-4), &PX(im, 0, h-4), w);
    memcpy(&PX(tmp, 0, h-2), &PX(im, 0, h-2), w);

    int r2 = valid[2]*valid[2]+valid[3]*valid[3];
    int brc = (valid[5]<<16) / r2;    // brightening coeff

    // horizontal pass
    for(int y = 0; y < h; y += 2)
    {
        int pxi, xc, yc = y-valid[1];

        xc = -valid[0]; pxi = PX(tmp, 0, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,0,y/2) = min(255, pxi);

        xc = 2-valid[0]; pxi = PX(tmp, 2, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,1,y/2) = min(255, pxi);

        uint8_t *p = &PX(tmp, 0, y);
        uint8_t *po = &PX2(im2x,2,y/2);

        for(int x = 4; x < w-4; x += 16)
        {
            uint8x16_t w1 = vld1q_u8(p);
            uint8x16_t w2 = vld1q_u8(p+16);
            p += 16;

            uint8x16x2_t ww = vuzpq_u8(w1, w2);
            w1 = ww.val[0]; w2 = ww.val[1];
            
            uint8x8_t v1l = vget_low_u8(w1);
            uint8x8_t v1h = vget_high_u8(w1);
            uint8x8_t v2l = vget_low_u8(w2);
            uint8x8_t v2h = vget_high_u8(w2);
            
            uint8x8_t v1 = vrhadd_u8(vrhadd_u8(v1l, v2l), vrhadd_u8(vext_u8(v1l, v1h, 1), vext_u8(v2l, v2h, 1)));
            uint8x8_t v2 = vrhadd_u8(vrhadd_u8(vext_u8(v2l, v2h, 2), vext_u8(v1l, v1h, 3)), vrhadd_u8(vext_u8(v2l, v2h, 3), vext_u8(v1l, v1h, 4)));
            v1 = vmin_u8(vmax_u8(v1, v2), vmax_u8(vmin_u8(v1, v2), vext_u8(v1l, v1h, 2)));

            xc = x - valid[0];
            int32x4_t vxc = vmovq_n_s32(xc);
            int32x4_t vyc = vmovq_n_s32(yc*yc);
            int32x4_t vxc1 = vaddq_s32(vxc, vmovl_s16(vcreate_s16(0x0006000400020000L)));
            int32x4_t vxc2 = vaddq_s32(vxc, vmovl_s16(vcreate_s16(0x000E000C000A0008L)));
            uint32x4_t vc1 = vreinterpretq_u32_s32(vaddq_s32(vmulq_s32(vxc1, vxc1), vyc));
            uint32x4_t vc2 = vreinterpretq_u32_s32(vaddq_s32(vmulq_s32(vxc2, vxc2), vyc));

            uint16x8_t v16x8 = vmovl_u8(v1);
            uint32x4_t v32x4l = vmovl_u16(vget_low_u16(v16x8));
            uint32x4_t v32x4h = vmovl_u16(vget_high_u16(v16x8));

            int32x4_t vbrc = vmovq_n_s32(brc);
            v1 = vmovn_u16(vcombine_u16(
                vmovn_u32(vminq_u32(vaddq_u32(v32x4l, vrshrq_n_u32(vmulq_u32(vshrq_n_u32(vmulq_u32(v32x4l, vc1), 16), vbrc), 8)), vmovq_n_u32(255))),
                vmovn_u32(vminq_u32(vaddq_u32(v32x4h, vrshrq_n_u32(vmulq_u32(vshrq_n_u32(vmulq_u32(v32x4h, vc2), 16), vbrc), 8)), vmovq_n_u32(255)))));

            vst1_u8(po, v1);
            po += 8;
        }

        xc = w-4-valid[0]; pxi = PX(tmp, w-4, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,w/2-2,y/2) = min(255, pxi);

        xc = w-2-valid[0]; pxi = PX(tmp, w-2, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,w/2-1,y/2) = min(255, pxi);
    }
}
#endif

#ifdef HAVE_AVX2
void removeSpeckles1_avx2(uint8_t *restrict im, uint8_t *restrict tmp, uint8_t *restrict im2x, int valid[6], int w, int h)
{
    // vertical pass
    memcpy(&PX(tmp, 0, 0), &PX(im, 0, 0), w);
    memcpy(&PX(tmp, 0, 2), &PX(im, 0, 2), w);

    for(int y = 4; y < h-4; y += 2)
    {
        __m256i v1, v2, v3, v4, v5, v6, v7, v8, v9;
        uint8_t *p = &PX(im, 0, y-4);
        uint8_t *po = &PX(tmp,0,y);

        for(int x = 0; x < w; x += 32)
        {
            v1 = _mm256_loadu_si256((__m256i *)p);         v2 = _mm256_loadu_si256((__m256i *)(p+w));     v3 = _mm256_loadu_si256((__m256i *)(p+2*w));
            v4 = _mm256_loadu_si256((__m256i *)(p+3*w));   v5 = _mm256_loadu_si256((__m256i *)(p+4*w));   v6 = _mm256_loadu_si256((__m256i *)(p+5*w));
            v7 = _mm256_loadu_si256((__m256i *)(p+6*w));   v8 = _mm256_loadu_si256((__m256i *)(p+7*w));   v9 = _mm256_loadu_si256((__m256i *)(p+8*w)); p += 32;
            v1 = _mm256_avg_epu8(_mm256_avg_epu8(v1, v2), _mm256_avg_epu8(v3, v4));
            v2 = _mm256_avg_epu8(_mm256_avg_epu8(v6, v7), _mm256_avg_epu8(v8, v9));
            v1 = _mm256_min_epu8(_mm256_max_epu8(v1, v2), _mm256_max_epu8(_mm256_min_epu8(v1, v2), v5));
            _mm256_storeu_si256((__m256i *)po, v1); po += 32;
        }
    }

    memcpy(&PX(tmp, 0, h-4), &PX(im, 0, h-4), w);
    memcpy(&PX(tmp, 0, h-2), &PX(im, 0, h-2), w);

    int r2 = valid[2]*valid[2]+valid[3]*valid[3];
    int brc = (valid[5]<<16) / r2;    // brightening coeff

    // horizontal pass
    for(int y = 0; y < h; y += 2)
    {
        int pxi, xc, yc = y-valid[1];

        xc = -valid[0]; pxi = PX(tmp, 0, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,0,y/2) = min(255, pxi);

        xc = 2-valid[0]; pxi = PX(tmp, 2, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,1,y/2) = min(255, pxi);

        __m256i v1, v2, v3, v4, v5, v6, v7, v8, v9; __m128i vv;
        uint8_t *p = &PX(tmp, 0, y);
        uint8_t *po = &PX2(im2x,2,y/2);

        __m256i vc0 = _mm256_cvtepu16_epi32(_mm_set_epi16(14, 12, 10, 8, 6, 4, 2, 0));

        for(int x = 4; x < w-4; x += 64)
        {
            v1 = _mm256_loadu_si256((__m256i *)p); v2 = _mm256_loadu_si256((__m256i *)(p+32)); p += 64; vv = _mm_loadu_si128((__m128i *)p);
            v3 = _mm256_permute2x128_si256(v1, v2, 0x20); v4 = _mm256_permute2x128_si256(v1, v2, 0x31);

            v1 = _mm256_packus_epi16(_mm256_and_si256(v3, _mm256_set1_epi16(255)), _mm256_and_si256(v4, _mm256_set1_epi16(255)));
            v2 = _mm256_packus_epi16(_mm256_srli_epi16(v3, 8), _mm256_srli_epi16(v4, 8));
            __m256i vv1 = _mm256_permute2x128_si256(v1, _mm256_castsi128_si256(_mm_packus_epi16(_mm_and_si128(vv, _mm_set1_epi16(255)), _mm_set1_epi8(0))), 0x21);
            __m256i vv2 = _mm256_permute2x128_si256(v2, _mm256_castsi128_si256(_mm_packus_epi16(_mm_srli_epi16(vv, 8), _mm_set1_epi8(0))), 0x21);

            v3 = _mm256_alignr_epi8(vv1, v1, 1); v4 = _mm256_alignr_epi8(vv2, v2, 1); v5 = _mm256_alignr_epi8(vv1, v1, 2); v6 = _mm256_alignr_epi8(vv2, v2, 2);
            v7 = _mm256_alignr_epi8(vv1, v1, 3); v8 = _mm256_alignr_epi8(vv2, v2, 3); v9 = _mm256_alignr_epi8(vv1, v1, 4);

            v1 = _mm256_avg_epu8(_mm256_avg_epu8(v1, v2), _mm256_avg_epu8(v3, v4));
            v2 = _mm256_avg_epu8(_mm256_avg_epu8(v6, v7), _mm256_avg_epu8(v8, v9));
            v1 = _mm256_min_epu8(_mm256_max_epu8(v1, v2), _mm256_max_epu8(_mm256_min_epu8(v1, v2), v5));

            xc = x - valid[0];
            __m256i vxc = _mm256_set1_epi32(xc);
            __m256i vyc = _mm256_set1_epi32(yc*yc);
            
            __m256i vxc1 = _mm256_add_epi32(vxc, vc0);
            __m256i vxc2 = _mm256_add_epi32(vxc, _mm256_add_epi32(vc0, _mm256_set1_epi32(16)));
            __m256i vxc3 = _mm256_add_epi32(vxc, _mm256_add_epi32(vc0, _mm256_set1_epi32(32)));
            __m256i vxc4 = _mm256_add_epi32(vxc, _mm256_add_epi32(vc0, _mm256_set1_epi32(48)));
            __m256i vc1 = _mm256_add_epi32(_mm256_mullo_epi32(vxc1, vxc1), vyc);
            __m256i vc2 = _mm256_add_epi32(_mm256_mullo_epi32(vxc2, vxc2), vyc);
            __m256i vc3 = _mm256_add_epi32(_mm256_mullo_epi32(vxc3, vxc3), vyc);
            __m256i vc4 = _mm256_add_epi32(_mm256_mullo_epi32(vxc4, vxc4), vyc);

            __m256i v16x16l = _mm256_cvtepu8_epi16(_mm256_extractf128_si256(v1, 0));
            __m256i v16x16h = _mm256_cvtepu8_epi16(_mm256_extractf128_si256(v1, 1));
            __m256i v32x8n1 = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(v16x16l, 0));
            __m256i v32x8n2 = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(v16x16l, 1));
            __m256i v32x8n3 = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(v16x16h, 0));
            __m256i v32x8n4 = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(v16x16h, 1));

            __m256i vbrc = _mm256_set1_epi32(brc);
            v32x8n1 = _mm256_min_epu32(_mm256_add_epi32(v32x8n1, _mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(_mm256_mullo_epi32(v32x8n1, vc1), 16), vbrc), _mm256_set1_epi32(128)), 8)), _mm256_set1_epi32(255));
            v32x8n2 = _mm256_min_epu32(_mm256_add_epi32(v32x8n2, _mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(_mm256_mullo_epi32(v32x8n2, vc2), 16), vbrc), _mm256_set1_epi32(128)), 8)), _mm256_set1_epi32(255));
            v32x8n3 = _mm256_min_epu32(_mm256_add_epi32(v32x8n3, _mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(_mm256_mullo_epi32(v32x8n3, vc3), 16), vbrc), _mm256_set1_epi32(128)), 8)), _mm256_set1_epi32(255));
            v32x8n4 = _mm256_min_epu32(_mm256_add_epi32(v32x8n4, _mm256_srli_epi32(_mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(_mm256_mullo_epi32(v32x8n4, vc4), 16), vbrc), _mm256_set1_epi32(128)), 8)), _mm256_set1_epi32(255));

            v32x8n1 = _mm256_packus_epi32(_mm256_permute2x128_si256(v32x8n1, v32x8n2, 0x20), _mm256_permute2x128_si256(v32x8n1, v32x8n2, 0x31));
            v32x8n2 = _mm256_packus_epi32(_mm256_permute2x128_si256(v32x8n3, v32x8n4, 0x20), _mm256_permute2x128_si256(v32x8n3, v32x8n4, 0x31));
            v1 = _mm256_packus_epi16(_mm256_permute2x128_si256(v32x8n1, v32x8n2, 0x20), _mm256_permute2x128_si256(v32x8n1, v32x8n2, 0x31));

            _mm256_storeu_si256((__m256i *)po, v1);
            po += 32;
        }

        xc = w-4-valid[0]; pxi = PX(tmp, w-4, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,w/2-2,y/2) = min(255, pxi);

        xc = w-2-valid[0]; pxi = PX(tmp, w-2, y);
        pxi += (((pxi * (xc*xc + yc*yc)) >> 16) * brc + 128) >> 8;
        PX2(im2x,w/2-1,y/2) = min(255, pxi);
        
    }
}
#endif

// speckle removal filter
// will not keep thin features, but will keep strong step-like changes in brightness
// this is the first pass, will only compute values at even cols/rows
void removeSpeckles1(uint8_t *restrict im, uint8_t *restrict tmp, uint8_t *restrict im2x, int valid[6], int w, int h)
{
#ifdef HAVE_NEON
    removeSpeckles1_neon(im, tmp, im2x, valid, w, h);
    return;
#endif

#ifdef HAVE_AVX2
    removeSpeckles1_avx2(im, tmp, im2x, valid, w, h);
    return;
#endif

    int x,y;

    const int d = 4;    // speckle removal distance

    // vertical pass
    for (y=0; y<d; y+=2)
        for (x=0; x<w; ++x)
            PX(tmp,x,y) = PX(im, x, y);

    for (y=d; y<h-d; y+=2)
        for (x=0; x<w; ++x)
        {
            int up = ( (PX(im, x, y-1) + PX(im, x, y-2) +1)/2 + (PX(im, x, y-3) + PX(im, x, y-4) +1)/2 +1 )/2;
            int dn = ( (PX(im, x, y+1) + PX(im, x, y+2) +1)/2 + (PX(im, x, y+3) + PX(im, x, y+4) +1)/2 +1 )/2;
            int pmax = max(up, dn);
            int pmin = min(up, dn);
            PX(tmp,x,y) = min(pmax, max(pmin, PX(im, x, y)));
        }

    for (y=h-d; y<h; y+=2)
        for (x=0; x<w; ++x)
            PX(tmp,x,y) = PX(im, x, y);

    int r2 = valid[2]*valid[2]+valid[3]*valid[3];
    int brc = (valid[5]<<16) / r2;    // brightening coeff

    // horizontal pass
    for (y=0; y<h; y+=2)
    {
        for (x=0; x<d/2; ++x)
        {
            int lt = PX(tmp, x*2, y);
            int rt = PX(tmp, x*2, y);
            HORZ_BODY(x*2)
            PX2(im2x,x,y/2) = pxi;
        }

        for (x=d/2; x<(w-d)/2; ++x)
        {
            int lt = ( (PX(tmp, 2*x-1, y) + PX(tmp, 2*x-2, y) +1)/2 + (PX(tmp, 2*x-3, y) + PX(tmp, 2*x-4, y) +1)/2  +1 )/2;
            int rt = ( (PX(tmp, 2*x+1, y) + PX(tmp, 2*x+2, y) +1)/2 + (PX(tmp, 2*x+3, y) + PX(tmp, 2*x+4, y) +1)/2  +1 )/2;
            HORZ_BODY(x*2)
            PX2(im2x,x,y/2) = pxi;
        }

        for (x=(w-d)/2; x<w/2; ++x)
        {
            int lt = PX(tmp, x*2, y);
            int rt = PX(tmp, x*2, y);
            HORZ_BODY(x*2)
            PX2(im2x,x,y/2) = pxi;
        }
    }
}


// this is the second pass, compute values at odd cols/rows within prescribed area only
void removeSpeckles2(uint8_t *restrict im, uint8_t *restrict tmp, int valid[6], int w, int h, int xs, int ys, int max_r)
{
    int x,y;

    const int d = 4;    // speckle removal distance

    int stx = max(0,xs-max_r);
    int enx = min(w,xs+max_r);
    int sty = max(0,ys-max_r);
    int eny = min(h,ys+max_r);

    // vertical pass
    for (y=sty|1; y<d; y+=2)
        for (x=stx; x<enx; ++x)
            PX(tmp,x,y) = PX(im, x, y);

    for (y=max(d+1,sty|1); y<min(h-d,eny); y+=2)
        for (x=stx; x<enx; ++x)
        {
            int up = ( (PX(im, x, y-1) + PX(im, x, y-2) +1)/2 + (PX(im, x, y-3) + PX(im, x, y-4) +1)/2 +1 )/2;
            int dn = ( (PX(im, x, y+1) + PX(im, x, y+2) +1)/2 + (PX(im, x, y+3) + PX(im, x, y+4) +1)/2 +1 )/2;
            int pmax = max(up, dn);
            int pmin = min(up, dn);
            PX(tmp,x,y) = min(pmax, max(pmin, PX(im, x, y)));
        }

    for (y=h-d+1; y<min(h,eny); y+=2)
        for (x=stx; x<enx; ++x)
            PX(tmp,x,y) = PX(im, x, y);

    int r2 = valid[2]*valid[2]+valid[3]*valid[3];
    int brc = (valid[5]<<16) / r2;    // brightening coeff

    // horizontal pass
    for (y=sty; y<eny; ++y)
    {
        for (x=stx; x<d; ++x)
        {
            int lt = PX(tmp, x, y);
            int rt = PX(tmp, x, y);
            HORZ_BODY(x)
            PX(im,x,y) = pxi;
        }

        for (x=max(d,stx); x<min(w-d,enx); ++x)
        {
            int lt = ( (PX(tmp, x-1, y) + PX(tmp, x-2, y) +1)/2 + (PX(tmp, x-3, y) + PX(tmp, x-4, y) +1)/2  +1 )/2;
            int rt = ( (PX(tmp, x+1, y) + PX(tmp, x+2, y) +1)/2 + (PX(tmp, x+3, y) + PX(tmp, x+4, y) +1)/2  +1 )/2;
            HORZ_BODY(x)
            PX(im,x,y) = pxi;
        }

        for (x=w-d; x<min(w,enx); ++x)
        {
            int lt = PX(tmp, x, y);
            int rt = PX(tmp, x, y);
            HORZ_BODY(x)
            PX(im,x,y) = pxi;
        }
    }
}
#undef HORZ_BODY


int thresholdForPupil(uint8_t *restrict im2x, int w, int h, float min_pup_r, int pupthr, int valid[6])
{
    int x,y,i;

    if (pupthr<0)
    {
        uint32_t hist[256];
        int sty = (valid[1]-valid[3])/2;
        int eny = (valid[1]+valid[3])/2;

        // find pupil pixel brightness threshold (histogram image, threshold at some cumulative percentage of pixels area)
        memset(hist, 0, 256*sizeof(uint32_t));
        for (y=sty; y<eny; ++y)
        {
            int yc = y-valid[1]/2;
            // exclude corners
            int stx = max( valid[0]/2 - (valid[4]/2-abs(yc)), (valid[0]-valid[2])/2 );
            int enx = min( valid[0]/2 + (valid[4]/2-abs(yc)), (valid[0]+valid[2])/2 );
            for (x=stx; x<enx; ++x)
                ++hist[PX2(im2x,x,y)];
        }

        uint32_t ch=0;
        uint32_t min_pup_area = (uint32_t)(F_PI*min_pup_r*min_pup_r/8);
        for (i=0; i<256; ++i)
        {
            ch += hist[i];
            if  (ch>min_pup_area) break;
        }

        // 75/64 = 1.17
        // *1.17 and +6 levels here to account for noise / nonuniformity
        pupthr = 75*(i+1)/64 + 6;
    }

    return pupthr;
}


// Note: either w, h should be divisible by 4 for this function or valid area should not touch frame edges
int getBigBlackArea(uint8_t *restrict im2x, uint8_t pupthr, int16_t *restrict digest, int16_t *restrict weight, int w, int h, int valid[6], int *xs, int *ys)
{
    int x,y,j;
    int ws, w1, w2;

    // operate on 4x sub-sampled data
    int stx = (valid[0]-valid[2]) / 4;
    int enx = (valid[0]+valid[2]) / 4;
    int sty = (valid[1]-valid[3]) / 4;
    int eny = (valid[1]+valid[3]) / 4;
    int wd = w/4;
    int hd = h/4;

    // horizontal digest of thresholded image: fill pixels with the length of continuous vertical strand they belong to
    memset(digest, 0, wd*hd*sizeof(int16_t));
    for (x=stx; x<enx; ++x)
    {
        int strand = 0;
        for (y=sty; y<eny; ++y)
        {
            int xc = x-valid[0]/4;
            int yc = y-valid[1]/4;
            // exclude corners
            if (abs(xc)+abs(yc) > valid[4]/4) continue;

            if ( PX2(im2x,x*2,y*2)<=pupthr )
            {
                if (strand==0) w1 = PX2(im2x,x*2,max(0,y-1)*2) - PX2(im2x,x*2,y*2);
                strand=strand+1;
            }
            else if (strand)
            {
                // apply weight, depending on how high is the minimum gradient at the ends of the strand
                // assumption is that pupil edge is well-defined,
                // erroneously thresholded non-pupil areas are unlikely to have such property
                w2 = PX2(im2x,x*2,y*2) - PX2(im2x,x*2,(y-1)*2);
                ws = min(4, max(0, (min(w1,w2)+2)/4-2 ));

                for (j=y-strand; j<y; ++j)
                {
                    digest[j*wd+x] = strand;
                    weight[j*wd+x] = ws;
                }
                strand = 0;
            }
        }
        if (strand)
        {
            w2 = PX2(im2x,x*2,min(hd,eny)*2) - PX2(im2x,x*2,(eny-1)*2);
            ws = min(4, max(0, (min(w1,w2)+2)/4-2 ));

            for (j=eny-strand; j<eny; ++j)
            {
                digest[j*wd+x] = strand;
                weight[j*wd+x] = ws;
            }
        }
    }

    // vertical scan: find horizontal strand that crosses vertical strands with the largest total number of pixels in these strands
    int maxstrand = 0;
    int maxarea = 0;
    *xs=-1; *ys=-1;
    for (y=sty; y<eny; ++y)
    {
        int strand = 0;
        int area = 0;
        int startx = 0;
        for (x=stx; x<enx; ++x)
        {
            if ( digest[y*wd+x] )
            {
                if (strand == 0) startx = x;
                strand += digest[y*wd+x] * weight[y*wd+x];
                area   += digest[y*wd+x];
            }
            else
            {
                if (strand > maxstrand)
                {
                    *xs = (x-1+startx)/2;
                    *ys = y;
                    maxstrand = strand;
                    maxarea = area;
                }
                strand = 0;
                area   = 0;
            }
        }
        // Note: this is optional as we do not expect pupil to be at the very edge
        if (strand > maxstrand)
        {
            *xs = (wd+startx)/2;
            *ys = y;
            maxstrand = strand;
            maxarea = area;
        }
    }

    // here: xs, ys - roughly the center of pupil
    *xs = *xs * 4 + 2;
    *ys = *ys * 4 + 2;

    // roughly the radius of a pupil
    return (int)(4*sqrtf(maxarea/F_PI));
}


void reThreshold(uint8_t *restrict im, uint8_t *restrict im2x, uint8_t *restrict imthr, int w, int h, int xs, int ys, int psr, uint8_t pupthr, float max_pup_r)
{
    int x,y,j;

    // use additional area around pupil
    int psr2 = max(32, (psr+1)*4);
    int stx = max(2,xs-psr2);
    int enx = min(w-1,xs+psr2);
    int sty = max(2,ys-psr2);
    int eny = min(h-1,ys+psr2);

    // we're using double the initially-detected pupil radius,
    // therefore there will always be 'transition' pixels (at the pupil edge)
    // and somewhat more non-pupil pixels in a tile
    uint32_t hist[256], histf[256];
    memset(hist, 0, 256*sizeof(uint32_t));

    int parea = 0;
    if (psr2>=max_pup_r/2)
    {
        // for large pupil - only use even coordinates since im is fully filtered/brightened at even rows/cols only
        stx &= ~1;
        sty &= ~1;

        for (y=sty/2; y<eny/2; ++y)
            for (x=stx/2; x<enx/2; ++x)
            {
                ++hist[PX2(im2x,x,y)];
                // pupil area under initial thresholding
                if (PX2(im2x,x,y) <= pupthr) ++parea;
            }
    }
    else
    {
        for (y=sty; y<eny; ++y)
            for (x=stx; x<enx; ++x)
            {
                ++hist[PX(im,x,y)];
                // pupil area under initial thresholding
                if (PX2(im,x,y) <= pupthr) ++parea;
            }
    }

    // smooth-out histogram
    for (j=0; j<256; ++j)
        histf[j] = ( hist[max(0,j-3)] + hist[max(0,j-2)] + hist[max(0,j-1)] + hist[j] + hist[min(255,j+1)] + hist[min(255,j+2)] + hist[min(255,j+3)] + hist[min(255,j+4)] + 4) / 8;

    // cumulative sum
    for (j=1; j<256;++j)
        histf[j] += histf[j-1];

    // once enough pixels under histogram cumulative sum to cover twice the pupil area
    // - we are definitely above the edge brightness
    for (j=pupthr; j<256; ++j)
        if (histf[j] > parea + histf[pupthr-1]) break;
    uint8_t pupthr2 = pupthr;
    if (j<256) pupthr2 = (3*pupthr + j+1) / 4;

    // re-thresholding the pupil area
    int max_r = (int)(max_pup_r+1);
    stx = max(0,xs-max_r)&(~1);
    enx = min(w,xs+max_r);
    sty = max(0,ys-max_r)&(~1);
    eny = min(h,ys+max_r);
    for (y=sty/2; y<eny/2; ++y)
        for (x=stx/2; x<enx/2; ++x)
            PX(imthr,x*2,y*2) = PX2(im2x,x,y) <= pupthr2;

    max_r = max_r/2+2;
    stx = max(0,xs-max_r);
    enx = min(w,xs+max_r);
    sty = max(0,ys-max_r);
    eny = min(h,ys+max_r);
    for (y=sty; y<eny; ++y)
        for (x=stx; x<enx; ++x)
            PX(imthr,x,y) = PX(im,x,y) <= pupthr2;
}


int extractEdgePixels(uint8_t *restrict imthr, int w, int h, int *xsp, int *ysp, int psr, int max_pup_r, int16_t *xy, float *pup_area)
{
    int x,y;
    int stx, enx, sty, eny;
    int stxpup, enxpup, stypup, enypup;
    int stxf, enxf, styf, enyf;
    int stxmax, enxmax, stymax, enymax;

    int xylen = 0;
    int parea = 0;

    int xs = *xsp;
    int ys = *ysp;

    // we can keep no more than that many coordinates
    // 4800 in case of 640x480 input frame
    int maxedges = (w/4)*(h/4)/4;

    // assuming that pieces of pupil are not torn apart by more than that, even value for optimized 2x decimation
    int d = max(8, 2*(5*psr/64));
    psr = max(16, psr);

    stxmax = max(d, xs-max_pup_r+1) & ~1;
    enxmax = min(w-d-2, xs+max_pup_r); enxmax -= (enxmax-stxmax) % d;
    stymax = max(d, ys-max_pup_r+1) & ~1;
    enymax = min(h-d-2, ys+max_pup_r); enymax -= (enymax-stymax) % d;

    // a simplistic version of flood-fill
    // not using a proper visit-all list-based approach since we have
    // 15x15 sites at most, and to cover all of them less than 10 scans are required

    // at least 4 pixels out of d x d should be flagged to assume square occupied (1 pixel for 2x decimation)
    int minpx = 4;

    stxpup = max(2, xs-psr); enxpup = min(w-2, xs+psr);
    stypup = max(2, ys-psr); enypup = min(h-2, ys+psr);

    // map of occupied squares (2nd bit)
    // and mark initial pupil area (3rd bit)
    // mark only one pixel in a site for speed
    stx = enxmax; enx = stxmax;
    sty = enymax; eny = stymax;
    for (y=stymax; y<=enymax; y+=d)
        for (x=stxmax; x<=enxmax; x+=d)
        {
            int noccupied = 0;
            for(int yy=0; yy<d; yy+=2)
                for(int xx=0; xx<d; xx+=2)
                    noccupied += PX(imthr,x+xx,y+yy);

            if (noccupied >= minpx/4) PX(imthr,x,y) |= 2;

            if ((x>=stxpup) && (x<=enxpup) && (y>=stypup) && (y<=enypup))
            {
                PX(imthr,x,y) |= 4;

                // initial bounding box
                if (stx>x) stx=x;
                if (enx<x) enx=x;
                if (sty>y) sty=y;
                if (eny<y) eny=y;
            }
        }

    // expand initially marked area in connecting sites
    while(1)
    {
        int nnew = 0;
        int stxn=stx, styn=sty, enxn=enx, enyn=eny;

        for (y=max(stymax,sty-d); y<=min(enymax,eny+d); y+=d)
            for (x=max(stxmax,stx-d); x<=min(enxmax,enx+d); x+=d)
            {
                if ((PX(imthr,x,y)&2) && ((PX(imthr,x,y)&4)==0))  // candidate to connect
                {
                    if ( (PX(imthr,x,max(stymax,y-d))&4) || (PX(imthr,max(stxmax,x-d),y)&4) || (PX(imthr,min(enxmax,x+d),y)&4) || (PX(imthr,x,min(enymax,y+d))&4) )
                    {
                        // if connecting at the very edge - pupil detection is incorrect
                        if ((x==stxmax) || (x==enxmax) || (y==stymax) || (y==enymax))
                            return 0;

                        PX(imthr,x,y) |= 4;
                        ++nnew;

                        // expanding bounding box
                        if (stxn>x) stxn=x;
                        if (enxn<x) enxn=x;
                        if (styn>y) styn=y;
                        if (enyn<y) enyn=y;
                    }
                }
            }

        stx=stxn; sty=styn; enx=enxn; eny=enyn;

        if (!nnew) break;
    }

    // reset unconnected pixels within bounding box
    for (y=sty; y<=eny; y+=d)
        for (x=stx; x<=enx; x+=d)
        {
            if (PX(imthr,x,y)&4)
                PX(imthr,x,y) -= PX(imthr,x,y)&6; // remove marks
            else
            {
                // ToDo: slight optimization possible: reset only every 2nd location outside fullres area
                for(int yy=y; yy<y+d; ++yy)
                    for(int xx=x; xx<x+d; ++xx)
                        PX(imthr,xx,yy) = 0;
            }
        }

    // increasing area slightly, as there might be small number of missed pixels at the very edge
    stx = stx-2; enx = enx+d+2;
    sty = sty-2; eny = eny+d+2;

    // boundaries of full-res pupil
    int fullres = 1;
    stxf = max(2, xs-max_pup_r/2); enxf = min(w-1, xs+max_pup_r/2);
    styf = max(2, ys-max_pup_r/2); enyf = min(h-1, ys+max_pup_r/2);
    if ( (stx<stxf) || (enx>enxf) || (sty<styf) || (eny>enyf) )
        fullres = 0;

    // center of mass
    int xw = 0;
    int yw = 0;
    if (fullres)
    {
        for (y=sty; y<eny; ++y)
            for (x=stx; x<enx; ++x)
            {
                parea += PX(imthr,x,y);
                xw += x * PX(imthr,x,y);
                yw += y * PX(imthr,x,y);
            }
    }
    else
    {
        for (y=sty; y<eny; y+=2)
            for (x=stx; x<enx; x+=2)
            {
                parea += PX(imthr,x,y);
                xw += x * PX(imthr,x,y);
                yw += y * PX(imthr,x,y);
            }
        parea *= 4;
        xw *= 4;
        yw *= 4;
    }
    xs = ( xw + parea/2) / parea;
    ys = ( yw + parea/2) / parea;

    if (fullres)
    {
        for (y=sty; y<eny && (xylen<maxedges); ++y)
            for (x=stx; (x<enx) && (xylen<maxedges); ++x)
                if(PX(imthr,x,y))
                {
                    int edge = min(2, (!PX(imthr,x-1,y)) + (!PX(imthr,x+1,y)) + (!PX(imthr,x,y-1)) + (!PX(imthr,x,y+1)) ) -
                            ((!PX(imthr,x,y+1) | !PX(imthr,x-1,y)) & ((x>=xs) & (y<=ys))) -
                            ((!PX(imthr,x,y-1) | !PX(imthr,x-1,y)) & ((x>=xs) & (y>=ys))) -
                            ((!PX(imthr,x,y+1) | !PX(imthr,x+1,y)) & ((x<=xs) & (y<=ys))) -
                            ((!PX(imthr,x,y-1) | !PX(imthr,x+1,y)) & ((x<=xs) & (y>=ys)));

                    if (edge>1)
                    {
                        xy[xylen*4]   = x-xs;
                        xy[xylen*4+1] = y-ys;
                        ++xylen;
                    }
                }
    }
    else
    {
        for (y=sty; y<eny && (xylen<maxedges); y+=2)
            for (x=stx; (x<enx) && (xylen<maxedges); x+=2)
                if(PX(imthr,x,y))
                {
                    int edge = min(2, (!PX(imthr,x-2,y)) + (!PX(imthr,x+2,y)) + (!PX(imthr,x,y-2)) + (!PX(imthr,x,y+2)) ) -
                            ((!PX(imthr,x,y+2) | !PX(imthr,x-2,y)) & ((x>=xs) & (y<=ys))) -
                            ((!PX(imthr,x,y-2) | !PX(imthr,x-2,y)) & ((x>=xs) & (y>=ys))) -
                            ((!PX(imthr,x,y+2) | !PX(imthr,x+2,y)) & ((x<=xs) & (y<=ys))) -
                            ((!PX(imthr,x,y-2) | !PX(imthr,x+2,y)) & ((x<=xs) & (y>=ys)));

                    if (edge>1)
                    {
                        xy[xylen*4]   = x-xs;
                        xy[xylen*4+1] = y-ys;
                        ++xylen;
                    }
                }
    }

    *pup_area = (float)parea;
    *xsp = xs;
    *ysp = ys;

    return xylen;
}


int angle_cmp(const void * elem1, const void * elem2)
{
    float ang1 = *((float*)((int16_t*)elem1+2));
    float ang2 = *((float*)((int16_t*)elem2+2));

    return ang1>ang2 ? 1:-1;
}


int cleanInternals(int16_t *xy, int xylen, float parea)
{
    int x,y,i,j;

    // approximation of pupil diameter
    float pup_dia = 2*sqrtf(parea/F_PI);

    // pre-computing atan2f will speed-up sorting
    for (i=0; i<xylen; ++i)
        *(float*)&xy[i*4+2] = atan2f(xy[i*4+1], xy[i*4]);

    // sort in ascending angular order
    qsort(xy, xylen, 4*sizeof(int16_t), angle_cmp);

    // remove inward-placed pixels
    int repeat;
    do
    {
        repeat = 0;

        // coordinate of previous point
        int xp0 = xy[(xylen-1)*4];
        int yp0 = xy[(xylen-1)*4+1];

        for (i=j=0; i<xylen; ++i)
        {
            int inxt = i+1 >= xylen ? 0 : i+1;
            int xn = xy[inxt*4];
            int yn = xy[inxt*4+1];

            x = xy[i*4];
            y = xy[i*4+1];

            // vectors to compute angles between
            int xp = xp0-x;
            int yp = yp0-y;
            xn -= x; yn -= y;

            // vector norms
            float normp = sqrtf( (float)(xp*xp+yp*yp) );
            float normn = sqrtf( (float)(xn*xn+yn*yn) );
            float normc = sqrtf( (float)(x*x+y*y) );

            // using an intermediate [-x -y] vector
            // it is not possible to decide if pixel is inward or outward placed with angles close to 180 degree (acos will return the same value)
            float a1 = acosf( (-x*xp-y*yp) / (normp * normc) );
            float a2 = acosf( (-x*xn-y*yn) / (normn * normc) );

            // the closer the distance between points to the approximate
            // pupil diameter, the smaller the largest acceptable angle
            float max_ang = 177 - 25 * fmaxf(0, fminf(1, fmaxf(normp, normn)/pup_dia - 0.25f ));

            // keep pixels if angle between them above 175 degree
            if ((a1+a2)*180 <= max_ang*F_PI)
            {
                xy[j*4] = x;
                xy[j*4+1] = y;
                ++j;

                xp0 = x;
                yp0 = y;
            }
            else   // if current is invalid - the previous point may become invalid as well
                repeat = 1;
        }
        xylen = j;
        
    } while (repeat);

    return xylen;
}


int EllipseDistance(float A[6], float *xy, int xylen, int min_r, int max_r, float *d_ell)
{
    // find distance norm to switch to pixel coordinate space
    // by computing how far is the center of the ellipse from its edges
    float cx, cy, semx, semy, tmp;
    int goodEllipse = EllipseGeomFromAlgebraic(A, &cx, &cy, &semx, &semy, &tmp);

    // if A does not correspond to an ellipse or ellipse is too small or too large
    if ((!goodEllipse) || (semx<min_r && semy<min_r) || (semx>max_r) || (semy>max_r)) return 0;

    float nrm = (semx+semy)/2 / fabsf( A[0]*cx*cx + A[1]*cx*cy + A[2]*cy*cy + A[3]*cx + A[4]*cy + A[5] );

    // fitting parameters are: ax^2 + bxy + cy^2 +dx + ey + f
    for (int i=0; i<xylen; ++i)
        d_ell[i] = nrm * ( A[0]*xy[2*i]*xy[2*i] + A[1]*xy[2*i]*xy[2*i+1] + A[2]*xy[2*i+1]*xy[2*i+1] + A[3]*xy[2*i] + A[4]*xy[2*i+1] + A[5] );

    return 1;
}


float ellipseFitAndDistance(float *xy, int xylen, int st1, int en1, int st2, int en2, float *d_ell, int maxdist, int min_r, int max_r, float A[6])
{
    // fit ellipse to the pixels at the pupil edge
    int fitfound = EllipseDirectFit(&xy[st1*2], en2-st1, en1-st1, st2-st1, A);
    if (!fitfound) return 0;

    if ( !EllipseDistance(A, xy, xylen, min_r, max_r, d_ell) ) return 0;

    float npx = 0;
    for (int i=0; i<xylen; ++i)
    {
        float dist = fabsf(d_ell[i]);

        // if pixel closer than maxdist - pixel is weighted higher
        // if perfectly on an ellipse - the point will count as 2
        npx += 1/(0.5f+dist/(2*maxdist));
    }

    return npx;
}

// try multiple 'sectoral' fits
// assuming either a part of pupil edge is occulded or 2 parts of pupil edge occluded (not more)
int pupilEllipse(float A[6], float *xy, int xylen, float *d_ell, int min_pup_r, int max_pup_r)
{
    int i,j;

    // dividing pixels into 8 circular sectors
    #define N_SEC 8
    int bins[N_SEC+1];

    bins[0] = 0;
    for (i=j=0; i<xylen; ++i)
    {
        float ang = atan2f(xy[i*2+1], xy[i*2]) + F_PI;
        // pixel list is already sorted in ascending order of angles
        if (ang>=(j+1)*2*F_PI/N_SEC) bins[++j] = i;
    }
    for(++j;j<=N_SEC; ++j) bins[j] = xylen;

    float maxgoodpx = 0;
    int maxnused = 0;

    // two sectors
    // Note: not quite covering all the possibilities as the circle start boundary never crossed
    int bs1, bs2, be1, be2;
    for (int s1=0; s1<N_SEC-1; ++s1)
    {
        for (int e1=s1; e1<N_SEC; ++e1)
        {
            for (int s2=e1+2; s2<N_SEC-1; ++s2)
            {
                for (int e2=s2; e2<N_SEC; ++e2)
                {
                    // use at least half of the circle
                    //if ((s2>e1+2) && (e2-s2+e1-s1<N_SEC/2)) continue;
                    if (e2-s2+e1-s1<N_SEC/2-1) continue;

                    int nused;
                    nused = bins[e1+1]-bins[s1] + bins[e2+1]-bins[s2];
                    if (nused<8) continue;     // not enough points

                    float ellipse[6];
                    float npx = ellipseFitAndDistance(xy, xylen, bins[s1], bins[e1+1], bins[s2], bins[e2+1], d_ell, 2, min_pup_r, max_pup_r, ellipse);

                    if (npx>maxgoodpx)
                    {
                        memcpy(A, ellipse, 6*sizeof(float));
                        maxgoodpx = npx;
                        maxnused = nused;
                        bs1=bins[s1];
                        bs2=bins[s2];
                        be1=bins[e1+1];
                        be2=bins[e2+1];
                    }
                }
            }
        }
    }

    // re-pack used xy coordinates (for calibration)
    for (i=bs1,j=0; i<be2; ++i)
    {
        if ((i>=be1) && (i<bs2)) continue;
        xy[j*2] = xy[i*2];
        xy[j*2+1] = xy[i*2+1];
        ++j;
    }

    return maxnused;
}


// returns how much the coordinate vector should be scaled
// depending on the angle, to change unit circle into given ellipse
float ellipseScale(float iris_dirs, float phi, float el_r)
{
    float x, y;
    sincosf(iris_dirs-phi, &y, &x);

    float sin_dirs_ell_2 = (y*y) / (x*x*el_r*el_r+y*y);
    // need limit to 1 here due to precision issues
    return fminf(1, sqrtf( 1 + (el_r*el_r-1)*sin_dirs_ell_2));
}


float camRayProjToMinor(float xyc[2], float xye[2], float iris_z)
{
    // angle between center of pupil/iris ellipse towards camera center and towards eyeball center
    float alpha = atan2f(xyc[1], xyc[0]) - atan2f(xye[1], xye[0]);

    // angle of camera ray hitting imaging plane at (x,y)
    // cra = atan2(norm(xyc), iris_z);

    // projected onto plane coinciding with minor ellipse axis and perpendicular to the imaging plane
    // ang = atan(tan(cra) * cos(alpha));
    return atanf(sqrtf(xyc[0]*xyc[0]+xyc[1]*xyc[1]) / iris_z * cosf(alpha));
}


void compensatePupilElongation(float *xy, int xylen, float el_r, float pup_ang, float pup_stretch)
{
    // reduce stretching when looking towards pupil at an angle
    // no stretching for angles above ~40 degree
    float stretch_adj = (pup_stretch-1) * fmaxf(0, 1-(1-el_r)*3);

    // apply stretch
    for (int i=0; i<xylen; ++i)
    {
        // distance from the stretch axis
        float d = stretch_adj * (xy[i*2] * -sinf(pup_ang) + xy[i*2+1] * cosf(pup_ang));

        xy[i*2]   -= d * -sinf(pup_ang);
        xy[i*2+1] -= d *  cosf(pup_ang);
    }
}


void  compensatePCA(float vec_to_eyeb[2], float iris_ratio, float pca_angle[2], float *el_r, float *phi_minor)
{
    // pupil 'gaze' vector from el_r, phi_minor
    float xy = sqrtf(1.0f - *el_r * *el_r);
    float pgv[3] = {sinf(*phi_minor)*xy, -cosf(*phi_minor)*xy, *el_r};

    // rotate gaze vector with pupillary circular axis angles
    // ToDo: move pca_mat construction into initialization
    float pca_mat[3][3];
    float pca_angle3[3] = {pca_angle[1], pca_angle[0], 0};
    rotationMat(pca_angle3, pca_mat, 0);

    float pgvr[4][3];
    int n = 2;
    pgvr[0][0] =   pca_mat[0][0]*pgv[0] + pca_mat[0][1]*pgv[1] + pca_mat[0][2]*pgv[2];
    pgvr[0][1] =   pca_mat[1][0]*pgv[0] + pca_mat[1][1]*pgv[1] + pca_mat[1][2]*pgv[2];
    pgvr[0][2] =   pca_mat[2][0]*pgv[0] + pca_mat[2][1]*pgv[1] + pca_mat[2][2]*pgv[2];
    pgvr[1][0] = - pca_mat[0][0]*pgv[0] - pca_mat[0][1]*pgv[1] + pca_mat[0][2]*pgv[2];
    pgvr[1][1] = - pca_mat[1][0]*pgv[0] - pca_mat[1][1]*pgv[1] + pca_mat[1][2]*pgv[2];
    pgvr[1][2] = - pca_mat[2][0]*pgv[0] - pca_mat[2][1]*pgv[1] + pca_mat[2][2]*pgv[2];

    // margin above which phi_minor direction may become confused
    float el_thr = 0.96f;
    if ( (iris_ratio > el_thr) || (*el_r > el_thr) )
    {
        n = 4;
        // pgv2 = possibility of el_minor/el_major swap
        float xy2 = 0.031619f; // sqrt(1-MAX_EL_RATIO^2);
        float pgv2[3] = {sinf(*phi_minor+F_PI/2)*xy2, -cosf(*phi_minor+F_PI/2)*xy2, MAX_EL_RATIO};
        pgvr[2][0] =   pca_mat[0][0]*pgv2[0] + pca_mat[0][1]*pgv2[1] + pca_mat[0][2]*pgv2[2];
        pgvr[2][1] =   pca_mat[1][0]*pgv2[0] + pca_mat[1][1]*pgv2[1] + pca_mat[1][2]*pgv2[2];
        pgvr[2][2] =   pca_mat[2][0]*pgv2[0] + pca_mat[2][1]*pgv2[1] + pca_mat[2][2]*pgv2[2];
        pgvr[3][0] = - pca_mat[0][0]*pgv2[0] - pca_mat[0][1]*pgv2[1] + pca_mat[0][2]*pgv2[2];
        pgvr[3][1] = - pca_mat[1][0]*pgv2[0] - pca_mat[1][1]*pgv2[1] + pca_mat[1][2]*pgv2[2];
        pgvr[3][2] = - pca_mat[2][0]*pgv2[0] - pca_mat[2][1]*pgv2[1] + pca_mat[2][2]*pgv2[2];
    }

    // normalize vector
    float xy_ref = sqrtf(1-iris_ratio*iris_ratio);
    float vn = xy_ref / sqrtf(vec_to_eyeb[0]*vec_to_eyeb[0] + vec_to_eyeb[1]*vec_to_eyeb[1]);
    float vec_to_eyeb3[3] = {vec_to_eyeb[0]*vn, vec_to_eyeb[1]*vn, iris_ratio};

    // choose phi_minor direction (out of 4) that best points towards eyeball center
    float best_dp = -1;
    int el_best = 0;
    for (int i=0; i<n; ++i)
    {
        // dot product
        float dp = vec_to_eyeb3[0]*pgvr[i][0] + vec_to_eyeb3[1]*pgvr[i][1] + vec_to_eyeb3[2]*pgvr[i][2];
        if (dp>best_dp) { best_dp = dp; el_best = i; }
    }

    // el_r, phi_minor from pupil 'gaze' vector
    *el_r = pgvr[el_best][2];
    *phi_minor = atan2f(pgvr[el_best][1], pgvr[el_best][0]) + F_PI/2;
}


// Refraction-induced correction of pupil ellipse parameters
//
// Values are derived from Zemax model assuming 7.75mm cornea curvature,
// pupil/iris 3.3mm behind, 4mm pupil diameter, refraction index of 1.3375
// Accounting for non-flat iris in real eyes, the offset values are lower
//
// Eq 3.5 formula for observed pupil ellipse ratio is given in
// D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 52:
// 0.99*cos((angle+5.3)/1.121), for this table - excluding 5.3 term: cos((0:5:90)*pi/180/1.121)'
//
// angle	ctr offset	--- Pupil ellipse ratio ---
//       Zemax       Zemax   cos(angle)  Eq 3.5a
//  0	0.000	    1.000   1.000       1.0000
//  5	0.048	    0.997   0.996       0.9970
// 10	0.094	    0.988   0.985       0.9879
// 15	0.143	    0.972   0.966       0.9729
// 20	0.193	    0.951   0.940       0.9519
// 25	0.245	    0.924   0.906       0.9252
// 30	0.299	    0.892   0.866       0.8929
// 35	0.357	    0.855   0.819       0.8552
// 40	0.419	    0.813   0.766       0.8123
// 45	0.484       0.768   0.707       0.7644  
// 50	0.554	    0.719   0.643       0.7120
// 55	0.63	    0.667   0.574       0.6552
// 60	0.711       0.612   0.500       0.5945
// 65	0.798	    0.556   0.423       0.5302
// 70	0.891	    0.497   0.342       0.4626
// 75	0.992	    0.437   0.259       0.3923
// 80	1.097	    0.375   0.174       0.3195
// 85                        0.087       0.2449                          
// 90                        0           0.1687
//
// true cos(angle) ratio from observed pupil ratio:             cos( acos(pu_r) * 1.121 )
// pupil-iris mm shift on imaging plane from true pupil ratio:
//  - due to pupil closer to eyeball center:                    -(iris2eye-pupil2eye) * sin( acos(pu_r) )
//  - due to refraction:                                        acos(pu_r) * 0.32 + ( acos(pu_r) * 0.58 )^2

void TruePupilFromRefraction(float semx, float semy, float phi, float pxl2mm, float eyeball2d[2], float pupil_refr, float *pca_angle, float iris_ratio, float iris_z, int ictrx, int ictry, float iris2pupil,
    float cxy[2], float *el_major, float *el_minor, float *phi_minor)
{
    *el_major = max(semx, semy);
    float el_r = min(semx, semy) / *el_major;

    // phi_minor = minor axis angle from vertical
    *phi_minor = phi + (semx < semy)*F_PI/2;

    // ----- Correct pupil ellipse widening due to refraction and camera projection

    // camera ray projected onto plane perpendicular to imaging plane and coinciding with minor ellipse axis
    float xyc[2] = {ictrx-cxy[0], ictry-cxy[1]};
    float xye[2] = {eyeball2d[0]-cxy[0], eyeball2d[1]-cxy[1]};
    float cra_proj = camRayProjToMinor(xyc, xye, iris_z);

    // perspective- and refraction- corrected ellipse ratio
    el_r = cosf(acosf(el_r*cosf(cra_proj))*pupil_refr - cra_proj);

    // minimize bias in ellipse ratio from pixelization-induced error (affects eye positions looking straight into the camera)
    // el_r = min(1, el_r+0.0025f);

    // correction for pupillary circular axis via 3d re-projection
    compensatePCA(xye, iris_ratio, pca_angle, &el_r, phi_minor);
    // phi_minor points toward 2d eyeball center, range: [-pi/2 .. 3*pi/2]

    *el_minor = *el_major * el_r;

    // ----- Center offset due to refraction

    // using stabilized ellipse ratio predicted from previous frame (large jitter in gaze angles otherwise)
    float mm_shift = -iris2pupil*sinf(acosf(iris_ratio)) + acosf(iris_ratio)*0.32f + (acosf(iris_ratio)*0.58f)*(acosf(iris_ratio)*0.58f);
    // float mm_shift = -iris2pupil*sinf(acosf(el_r)) + acosf(el_r)*0.32f + (acosf(el_r)*0.58f)*(acosf(el_r)*0.58f);
    float d[2];
    d[0] =  mm_shift/pxl2mm*sinf(*phi_minor);
    d[1] = -mm_shift/pxl2mm*cosf(*phi_minor);
    cxy[0] += d[0]; cxy[1] += d[1];
}


// Center offset due to pupil decenter
void correctPupilDecenter(float cxy[2], float el_r, float phi, float pxl2mm, float decenter_mm[2], float camera_rot, float ixy[2])
{
    float decenter[2];
    decenter[0] = decenter_mm[0] / pxl2mm;
    decenter[1] = decenter_mm[1] / pxl2mm;
    // rotate decenter to camera space
    rotcam(decenter, -camera_rot, decenter);
    // reduce decentering according to ellipse scale in horizontal/vertical direction
    float el_xy[2];
    el_xy[0] = ellipseScale(0-camera_rot*F_PI/180, phi, el_r);
    el_xy[1] = ellipseScale(F_PI/2-camera_rot*F_PI/180, phi, el_r);
    ixy[0] = cxy[0] - decenter[0]*el_xy[0];
    ixy[1] = cxy[1] - decenter[1]*el_xy[1];
}

#define FAST_H_Z

// rasterize iris sector in polar coordinates
void rasterizeIris(uint8_t * restrict irim, int irim_w, float *iris_d, float *iris_scl, uint8_t *im2x, int w, int h, float cxy[2],
    float min_r, float iris_z, float phi_minor, float camera_rot)
{
    float h_min = min_r;
    float h_max = min_r+irim_w;
    
    for (int i=0; i<IRIM_H; ++i)
    {
        float dxy0[2] = {iris_scl[i] * cosf(iris_d[i]), iris_scl[i] * sinf(iris_d[i]) };
        float dxy0rot[2], dxys[2];

        // correcting for human iris being 11.5x10.5 mm
        // The cornea is transparent and it's size is about 11.5mm h x 10.5mm v (= 0.913 ratio)
        // D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 3
        rotcam(dxy0, camera_rot, dxy0rot);
        dxy0rot[1] *= 0.95f;
        rotcam( dxy0rot, -camera_rot, dxys);

        // phi_minor always points away from camera (towards eyeball center)
        // range of iris_d: [-pi/6 .. pi-pi/6]
        // rrange of phi_minor: [-pi/2 .. 3*pi/2]
        // phi_minor origin is at 12 o'clock; iris_d origin is at 3 o'clock
        // +pi/2 to match origins; +2*pi to prevent negative values in mod
        float ang_diff = fmodf(iris_d[i]-phi_minor+F_PI/2+2*F_PI, 2*F_PI);
        // negative height - this point on iris is closer to the camera than the imaging plane
        int neg = 1 - ((ang_diff < F_PI/2) || (ang_diff > 2*F_PI-F_PI/2))*2;

        float h_scl = sqrtf(1.0f - iris_scl[i]*iris_scl[i]);

        // h_z calculation that avoids divide by value changing in a loop at a cost of 1% error
        // (which for most part disappear in rounding due to correction being about 10 pixels or less)
        float h_z_min = 1 - iris_z / (h_min*neg*h_scl + iris_z);
        float h_z_max = 1 - iris_z / (h_max*neg*h_scl + iris_z);
        float h_z = h_z_min;
        float h_z_step = (h_z_max-h_z_min)/irim_w;
        for (int j=0; j<irim_w; ++j, h_z+=h_z_step)
        {
            float x = cxy[0] + (min_r+j) * dxys[0];
            float y = cxy[1] + (min_r+j) * dxys[1];

            // camera projection correction
#ifndef FAST_H_Z
            float height = (float)((min_r+j)*neg) * h_scl;
            float h_z = height / (height + iris_z);
#endif
            x -= h_z * (x-w/2);
            y -= h_z * (y-h/2);

            // only use even coordinates since im is fully filtered/brightened at even rows/cols only
            x = fmaxf(0, fminf(x, w-1.0f));
            y = fmaxf(0, fminf(y, h-1.0f));
            irim[i*irim_w+j] = PX2(im2x, (int)((x+1)/2), (int)((y+1)/2));
        }
    }
}


float quad_fit(float ym1, float y0, float y1)
{
    // find parabola extremum (where derivative of y=ax^2+bx+c goes to zero)
    if (y0>(ym1+y1)/2)
        return (ym1 - y1) / (2*ym1 + 2*y1 - 4*y0);
    else
    {
        // if parabola seem to be inverted (or just straight line) - return location of maximum
        if (y1>ym1) return 1;
            else return -1;
    }
}


// Simple method of iris radius detection: find maximum contrast break on polar image of iris
//
// The line can be slanted up to ~1/15th of the raster width (with -pi/6..pi/6 iris raster scan)
// due to up to 0.9mm decenter of pupil relative to the iris
// (assuming raster width about equals iris radius =0.9*sin(pi/6)/iris_r_mm)
// This smears the peak a bit but still provides good enough peak amplitude
// 
void detectIrisRadius(uint8_t *restrict irim, uint8_t *restrict tmp, int w, int h, float min_iris_r, float el_major, int *ir_r, int *irpeak, float *decent_v)
{
    int pupil_edge = max(0, (int)(el_major - min_iris_r + 0.5f));

    // contrast window length
    int tlen = (w+9)/10;

    // slant step (to detect vertical decentration)
    int slstep = (w+14)/15;

    // partition temporary storage to hold needed arrays
    memset(tmp, 0, 5*w*sizeof(int));
    int *restrict brtl  = (int *)tmp;
    int *restrict brtr  = brtl+w;
    int *restrict devr  = brtr+w;
    int *restrict meanv = devr+w;
    int *restrict devv  = meanv+w;

    // eliminate glints
    for (int y=0; y<h; ++y)
        for (int x=0; x<w; ++x)
            meanv[x] += PX(irim, x, y); // mean over iris scan verticals
    for (int x=0; x<w; ++x)
    {
        devv[x] = meanv[x]*5/(4*h);     // re-using devv stroage for thresholds
        meanv[x] /= h;
    }
    for (int y=0; y<h; ++y)
        for (int x=0; x<w; ++x)
            if (PX(irim, x, y) >= devv[x]) PX(irim, x, y) = meanv[x];

    // mean and deviation from mean over iris scan verticals
    // multiplied by h to keep precision
    for (int x=0; x<w; ++x)
    {
        meanv[x] = 0;
        devv[x] = 0;
    }

    int sum_meanv=0;
    for (int y=0; y<h; ++y)
        for (int x=0; x<w; ++x)
        {
            meanv[x] += PX(irim, x, y);
            sum_meanv += PX(irim, x, y);
        }

    for (int y=0; y<h; ++y)
        for (int x=0; x<w; ++x)
            devv[x] += abs(PXI(irim, x, y)*h - meanv[x]);
    for (int x=0; x<w; ++x)
        devv[x] /= h;

    
    // sliding window
    for (int i=0; i<tlen; ++i)
    {
        brtl[tlen-1] += meanv[i];
        brtr[tlen-1] += meanv[i+tlen];
        devr[tlen-1] += devv[i+tlen];
    }
    for (int i=tlen+pupil_edge; i<w-tlen; ++i)
    {
        brtl[i] = brtl[i-1] + meanv[i] - meanv[i-tlen];
        brtr[i] = brtr[i-1] + meanv[i+tlen] - meanv[i];
        devr[i] = devr[i-1] + devv[i+tlen] - devv[i];
    }

    int peak=0;
    *ir_r=0;
    for (int i=tlen+pupil_edge; i<w-tlen; ++i)
    {
        int cpeak = 2*brtr[i]-devr[i] - 2*brtl[i];
        if (cpeak > peak)
        {
            *ir_r = i;
            peak = cpeak;
        }
    }
    // looking for at least 10% diff in brightness (additional divider of 10 is from tlen=w/10)
    *irpeak = max(0, peak-(sum_meanv+50)/100);
    
    if ( (*ir_r<=tlen+pupil_edge) || (*ir_r>=w-tlen-1) || (*irpeak==0) )
    {
        *ir_r=0;
        *irpeak=0;
        *decent_v=0;
    }
    else
    {
        int peak_slant[2];
        for (int s=0; s<2; ++s)
        {
            float slant = (float)(s*2-1)*slstep;

            // mean and deviation from mean over iris scan slanted verticals
            // multiplied by h to keep precision
            int dx;
            for (int x=*ir_r-tlen; x<*ir_r+tlen; ++x)
            {
                meanv[x] = 0;
                devv[x] = 0;
            }

            for (int y=0; y<h; ++y)
                for (int x=*ir_r-tlen; x<*ir_r+tlen; ++x)
                {
                    dx = (int)(x + slant*(y-h/2)/(h/2) + 0.5f);
                    dx = max(0, min(w-1, dx));
                    meanv[x] += PX(irim, dx, y);
                }

            for (int y=0; y<h; ++y)
                for (int x=*ir_r-tlen; x<*ir_r+tlen; ++x)
                {
                    dx = (int)(x + slant*(y-h/2)/(h/2) + 0.5f);
                    dx = max(0, min(w-1, dx));
                    devv[x] += abs(PXI(irim, dx, y)*h - meanv[x]);
                }

            for (int x=*ir_r-tlen; x<*ir_r+tlen; ++x)
                devv[x] /= h;

            int brtls=0, brtrs=0, devrs=0;
            for (int i=*ir_r-tlen; i<*ir_r; ++i)
            {
                brtls += meanv[i];
                brtrs += meanv[i+tlen];
                devrs += devv[i+tlen];
            }

            peak_slant[s] = 2*brtrs-devrs - 2*brtls;
        }

        // decenter here is measured in horizontal iris raster pixels,
        // which equals camera pixels along major ellipse axis
        // vertical decenter computed with formula below is giving somewhat reduced values for some reason
        float peak_fit = quad_fit((float)peak_slant[0], (float)peak, (float)peak_slant[1]);
        *decent_v = slstep * max(-1, min(1, peak_fit )) / sinf(IRIS_SCTR/2);
        *ir_r += (int)(min_iris_r+0.5f);
    }
}


void getIris(eye_cfg *ei, float ir_r1, float ir_r2, float peak1, float peak2, float iris2d_p[2], int iris_recenter, float decent_v1, float decent_v2, float el_ratio,
    float *ir_r, float *iris_r1, float *iris_r2, float *peak, float iris2d[2], float *phi_minor)
{
    // initialize from previous frame values
    *ir_r = ei->iris_r;
    *peak = ei->iris_peak;
    iris2d[0] = iris2d_p[0];
    iris2d[1] = iris2d_p[1];

    // if iris detected - update iris radius
    float ir_r_new, peak_new;
    if (ir_r1>0 || ir_r2>0)
    {
        // reject much weaker edge altogether
        if (peak1>4*peak2) { peak2=0; ir_r2=0; }
        if (peak2>4*peak1) { peak1=0; ir_r1=0; }
        // reject high difference in radii (>40% = 0.95mm pupil decenter)
        if ( fminf(ir_r1, ir_r2)*1.4f < fmaxf(ir_r1, ir_r2) )
        {
            if (peak1<peak2) ir_r1=0;
            else ir_r2=0;
        }

        if (ir_r1>0 && ir_r2>0)
        {
            ir_r_new = (ir_r1*peak1+ir_r2*peak2)/(peak1+peak2);
            peak_new = (peak1+peak2)/2;
        }
        else
        {
            if (ir_r1>0)
            {
                ir_r_new = ir_r1;
                peak_new = peak1;
            }
            else
            {
                ir_r_new = ir_r2;
                peak_new = peak2;
            }
        }

        // update iris radius, weighting in the current peaks vs average observed peaks
        if (*peak) peak_new = fminf(peak_new, *peak*4);
        *ir_r = (((1<<ei->smooth_iris)-1)* *ir_r * *peak + ir_r_new*peak_new) / (((1<<ei->smooth_iris)-1) * *peak + peak_new);
        *peak = (((1<<ei->smooth_iris)-1)* *peak + peak_new) / (1<<ei->smooth_iris);
    }
    else
    {
        // if no iris detection in current frame - compute expected iris radius from new distance from eyeball center
        float iris2d_mm[2]      = {iris2d[0] * ei->pxl2mm, iris2d[1] * ei->pxl2mm};
        float iris2d_prev_mm[2] = {ei->iris2d_prev[0] * ei->pxl2mm, ei->iris2d_prev[1] * ei->pxl2mm};
        float d_prv_mm2 = norm_euc2(ei->eyeball3d_cam, iris2d_prev_mm);
        float d_new_mm2 = norm_euc2(ei->eyeball3d_cam, iris2d_mm);
        // using larger iris-to-eyeball-center distance to avoid error magnification
        // in z due to noisy x/y measurements when angle from camera is high
        float iris2eye = fmaxf(ei->iris2hrot_mm, ei->iris2vrot_mm);
        float iris2eye2 = iris2eye*iris2eye;
        float z_prv = ei->iris_z + ( iris2eye - sqrtf(fmaxf(0, iris2eye2 - d_prv_mm2)) ) / ei->pxl2mm;
        float z_new = ei->iris_z + ( iris2eye - sqrtf(fmaxf(0, iris2eye2 - d_new_mm2)) ) / ei->pxl2mm;
        *ir_r *= z_prv/z_new;
    }

    // re-position horizontal iris center from detected iris curves
    // Note: effectively negates coordinate correction in TruePupilFromRefraction from pupil refraction and decenter
    if (ir_r1 && ir_r2)
    {
		if (iris_recenter >= 1)
            iris2d[0] += ellipseScale(0, *phi_minor, el_ratio) * (ir_r2-ir_r1)/2;
		if (iris_recenter >= 2)
            iris2d[1] += ellipseScale(0, *phi_minor+F_PI/2, el_ratio) * (decent_v1+decent_v2)/2;
        
        // rotate phi_minor according to inclination of iris edges
		if (iris_recenter >= 3)
            *phi_minor += asinf((decent_v2-decent_v1) / *ir_r);
    }

    // used in calibration
    *iris_r1 = ir_r1;
    *iris_r2 = ir_r2;
}

int time_cmp(const void * elem1, const void * elem2)
{
    int t1 = ((timedbins*)elem1)->time;
    int t2 = ((timedbins*)elem2)->time;

    return t1>t2 ? 1:-1;
}

void Estimate2dEyeCenter(eye_cfg *ei, float iris_ratio, float phi_minor, float eyeball2d[2], float iris2ctr_pxl, float exy[2], float *eyeball_confid)
{
    float quals = 1;
    exy[0] = eyeball2d[0]; exy[1] = eyeball2d[1];
    *eyeball_confid = 0.1f;

    // reduce weight of old observations: halving after 3 seconds
    float decay = 1-0.2f/ei->fps_ec2d;
    for (int i=0; i<N_EYE_OBSERVATIONS; ++i)
    {
        ei->ctr_weight[i][0] *= decay;
        ei->ctr_weight[i][1] *= decay;
        ++ei->obs_time[i];
    }

    // iris not detected or pupil not fully visible or too low angle to camera - do not compute / update eyeball center
    // 0.175=10, 0.2=11.5 0.25=14.5 degree; 0.99=8  0.98=11.5  0.965=15.2 degree
    float iris2ctr[2] = { (ei->ctr_tail[0]-ei->iris2d[0])*ei->iris2vrot_mm/ei->iris2hrot_mm, ei->ctr_tail[1]-ei->iris2d[1] };
    if ( ei->iris_r && (ei->blink==0) && (iris_ratio < 0.98f) && ( (iris2ctr[0]*iris2ctr[0] + iris2ctr[1]*iris2ctr[1]) / (iris2ctr_pxl*iris2ctr_pxl) > 0.22f*0.22f) )
    {
        // derive iris center to eyeball center distance in 2d space by computing trigonometry from el_minor/el_major
        // Note: this method has precision issues if single frame only is used:
        //  - low precision for low angles as angle precision from acos is very low if el_minor is close to el_major
        //  - low precision for very high angles due to refraction effects causing elongation of el_minor (correction done in TruePupilFromRefraction)
        float below15deg = fmaxf(0.1f, fminf(1, (0.98f-iris_ratio)*20 ));
        float above30deg = fmaxf(0.4f, fminf(1, (iris_ratio-0.76f)*10 ));
        float occluded = fmaxf(0.1f, (1-3*ei->blink));

        // origin and direction of 2d line from current observation
        float dv[2];     // direction vector
        dv[0] =  sinf(phi_minor);
        dv[1] = -cosf(phi_minor);

        float proj_ang = acosf(iris_ratio);
        float len2dy = iris2ctr_pxl*sinf(proj_ang);
        // horizontal rotation has axis farther then eyeball center (about 12mm from iris center)
        float len2dx = ei->iris2hrot_mm / ei->iris2vrot_mm * len2dy;

        // eyeball center in 2d is at the len2d distance from the iris center
        // in the direction of minor ellipse axis (one of two possible locations)
        // e1v/e2v - possible locations of eyeball center (coincides with point on rotation axis for the vertical gaze direction)
        float iris2hrot[2], iris2vrot[2];
        float e1[2], e2[2], e1v[2], e2v[2];
        iris2hrot[0] = len2dx*dv[0];  iris2hrot[1] = len2dx*dv[1];
        iris2vrot[0] = len2dy*dv[0];  iris2vrot[1] = len2dy*dv[1];
        e1v[0] = ei->iris2d[0] + iris2vrot[0];  e1v[1] = ei->iris2d[1] + iris2vrot[1];
        e2v[0] = ei->iris2d[0] - iris2vrot[0];  e2v[1] = ei->iris2d[1] - iris2vrot[1];
        e1[0]  = ei->iris2d[0] + iris2hrot[0];  e1[1]  = ei->iris2d[1] + iris2vrot[1];
        e2[0]  = ei->iris2d[0] - iris2hrot[0];  e2[1]  = ei->iris2d[1] - iris2vrot[1];

        // - distance from designed 'perfect' center
        //   considering the eyeball center farther than twice the iris to-eyeball distance (about 24mm) from designed position very unlikely
        float w1d = fmaxf(0.1f, fminf(1, 1 - fmaxf(0, sqrtf(norm_euc2(e1, ei->eyeball_ctr))/iris2ctr_pxl-1) ));
        float w2d = fmaxf(0.1f, fminf(1, 1 - fmaxf(0, sqrtf(norm_euc2(e2, ei->eyeball_ctr))/iris2ctr_pxl-1) ));
        // - distance from previous center detection
        float w1p = fmaxf(0.2f, 1 - sqrtf(norm_euc2(e1v, eyeball2d))/iris2ctr_pxl );
        float w2p = fmaxf(0.2f, 1 - sqrtf(norm_euc2(e2v, eyeball2d))/iris2ctr_pxl );

        // if observation is of good quality - select which of the previous observations to replace
        // quality reduction from:
        // - low weights above
        // - highly occluded pupil (blink)
        // - high or low pupil angles
        float w1 = w1p*w1d * occluded * below15deg * above30deg;
        float w2 = w2p*w2d * occluded * below15deg * above30deg;
        float qual = fmaxf(w1, w2);                                            // quality to save to the bin
        quals = occluded * below15deg * fmaxf(w1p, w2p) * ei->slippage_rate;   // quality metric for smoothing

        // prevent same direction from filling all the bins by limiting range of bins to occupy according to radial direction
        // max(0, ...) here just to make MSVC happy acessing ctr_weight[k]
        int binstart = max( 0, ( (int)(phi_minor/(2*F_PI)*N_EYE_OBSERVATIONS-N_EYE_OBSERVATIONS/8 +2*N_EYE_OBSERVATIONS+0.5) ) & (N_EYE_OBSERVATIONS-1) );

        // find empty bin
        int obin = -1;
        int j,k;
        for (j=k=binstart; j<binstart+N_EYE_OBSERVATIONS/4; ++j, k=(k+1)&(N_EYE_OBSERVATIONS-1))
            if (!ei->valid[k]) { obin = k; break; }

        if (obin < 0)
        {
            // no empty bins - replace bin with the lowest quality and similar radial direction
            float minqual = qual;
            for (j=k=binstart; j<binstart+N_EYE_OBSERVATIONS/4; ++j, k=(k+1)&(N_EYE_OBSERVATIONS-1))
            {
                float qualbin = fmaxf(ei->ctr_weight[k][0], ei->ctr_weight[k][1]);
                float qualnew = qual * fmaxf( 0, cosf(phi_minor - ei->phi_obs[k]) * 2 );
                if (qualnew-qualbin > minqual)
                {
                    minqual = qualnew-qualbin;
                    obin = k;
                }
            }
        }

        if (obin >= 0)
        {
            ei->obs_time[obin] = 0;
            ei->valid[obin] = 1;
            ei->phi_obs[obin] = phi_minor;
            ei->ctr_obs[obin][0][0] = e1[0];  ei->ctr_obs[obin][0][1] = e1[1];  ei->ctr_obs[obin][1][0] = e2[0];  ei->ctr_obs[obin][1][1] = e2[1];
            ei->ctr_weight[obin][0] = w1; ei->ctr_weight[obin][1] = w2;
            ei->iris2d_obs[obin][0] = ei->iris2d[0]; ei->iris2d_obs[obin][1] = ei->iris2d[1];
        }

        // construct averaged eyeball position
        // from a set of prior observations of possible eyeball center locations and current one

        int n_valid = 0;
        for (j=0; j<N_EYE_OBSERVATIONS; ++j) n_valid += ei->valid[j];

        // select cluster centers based on both current detection and previous smoothed one
        float c1[2] = {e1[0], e1[1]};
        float c2[2] = {e2[0], e2[1]};

        if (n_valid == N_EYE_OBSERVATIONS)
        {
            if (norm_euc2(e1, ei->ctr_tail) < norm_euc2(e2, ei->ctr_tail))
            {
                c1[0] = (e1[0] + ei->ctr_tail[0] *3) / 4;
                c1[1] = (e1[1] + ei->ctr_tail[1] *3) / 4;
            }
            else
            {
                c2[0] = (e2[0] + ei->ctr_tail[0] *3) / 4;
                c2[1] = (e2[1] + ei->ctr_tail[1] *3) / 4;
            }
        }

        // limit to fresh observations only if enough high-quality observations with wide enough angular spread
        // sort from fresh to old
        for (int i=0; i<N_EYE_OBSERVATIONS; ++i)
        {
            ei->timesorted_bins[i].time = ei->obs_time[i];
            ei->timesorted_bins[i].bin = i;
        }
        qsort(ei->timesorted_bins, N_EYE_OBSERVATIONS, sizeof(timedbins), time_cmp);
        memset(ei->valid_bins, 0, N_EYE_OBSERVATIONS*sizeof(int));
        int bin = ei->last_bin;
        float min_ang = 0;
        if (bin >= 0) min_ang = ei->ang_obs[bin];
        float max_ang = min_ang;
        int good_ang_spread = 0;
        int n = 0;
        for (int i=0; i<N_EYE_OBSERVATIONS; ++i)
        {
            bin = ei->timesorted_bins[i].bin;
            if (ei->valid[bin])
            {
                ei->valid_bins[bin] = 1;
                if (fmaxf(ei->ctr_weight[bin][0], ei->ctr_weight[bin][1]) > 0.7f)       // limit to high-quality observations only
                {
                    // check if new observation extends angular spread
                    if (fmodf(ei->ang_obs[bin]-max_ang, 2*F_PI) < F_PI) max_ang = ei->ang_obs[bin];
                    if (fmodf(min_ang-ei->ang_obs[bin], 2*F_PI) < F_PI) min_ang = ei->ang_obs[bin];
                    if (n>=8 && fmodf(max_ang-min_ang, 2*F_PI) > 45*F_PI/180)
                    {
                        good_ang_spread = 1;
                        break;
                    }
                }
                ++n;
            }
        }

        if (obin>=0) ei->last_bin = obin;

        // divide into clusters by comparing distance to all previous observations
        // and find which cluster is more compact
        for (j=0,k=0; j<N_EYE_OBSERVATIONS; ++j)
            if (ei->valid_bins[j])
            {
                ei->obs[k][0] = ei->ctr_obs[j][0][0];
                ei->obs[k][1] = ei->ctr_obs[j][0][1];
                ei->weight[k] = ei->ctr_weight[j][0];
                ++k;
                ei->obs[k][0] = ei->ctr_obs[j][1][0];
                ei->obs[k][1] = ei->ctr_obs[j][1][1];
                ei->weight[k] = ei->ctr_weight[j][1];
                ++k;
            }
        
        n_valid = k/2;
        int n_valid2 = k;
        float meand1=0;
        float meand2=0;
        for (j=0; j<n_valid2; ++j)
        {
            ei->d1[j] = norm_euc2(c1, ei->obs[j]);
            ei->d2[j] = norm_euc2(c2, ei->obs[j]);
            meand1 += ei->d1[j];
            meand2 += ei->d2[j];
        }
        meand1 /= n_valid2;
        meand2 /= n_valid2;

        for (j=0; j<n_valid2; ++j)
        {
            ei->near1[j] = ei->d1[j] <= meand1;
            ei->near2[j] = ei->d2[j] <= meand2;
        }

        // re-center clusters and measure cluster size as mean distance from it's center
        int k1=0, k2=0;
        c1[0]=c1[1]=c2[0]=c2[1]=0;
        for (j=0; j<n_valid2; ++j)
        {
            if (ei->near1[j])
            {
                c1[0] += ei->obs[j][0];
                c1[1] += ei->obs[j][1];
                ++k1;
            }
            if (ei->near2[j])
            {
                c2[0] += ei->obs[j][0];
                c2[1] += ei->obs[j][1];
                ++k2;
            }
        }
        if (k1) { c1[0] /= k1; c1[1] /= k1; }
        if (k2) { c2[0] /= k2; c2[1] /= k2; }

        for (j=0; j<n_valid2; ++j)
        {
            ei->d1[j] = norm_euc2(c1, ei->obs[j]);
            ei->d2[j] = norm_euc2(c2, ei->obs[j]);
        }

        // discard outliers in a cluster
        float meanD1near1=0;
        float meanD2near2=0;
        k1=k2=0;
        for (j=0; j<n_valid2; ++j)
        {
            if (ei->near1[j])
            {
                meanD1near1 += ei->d1[j];
                ++k1;
            }
            if (ei->near2[j])
            {
                meanD2near2 += ei->d2[j];
                ++k2;
            }
        }
        if (k1) meanD1near1 /= k1;
        if (k2) meanD2near2 /= k2;

        if (k1 >=8)
        {
            k1 = 0;
            for (j=0; j<n_valid2; ++j)
                ei->near1[j] = ei->d1[j] <= 2*meanD1near1;
            meanD1near1 = 0;
            for (j=0; j<n_valid2; ++j)
                if (ei->near1[j]) { meanD1near1 += ei->d1[j]; ++k1; }
            if (k1) meanD1near1 /= k1;
        }
        if (k2 >=8)
        {
            k2 = 0;
            for (j=0; j<n_valid2; ++j)
                ei->near2[j] = ei->d2[j] <= 2*meanD2near2;
            meanD2near2 = 0;
            for (j=0; j<n_valid2; ++j)
                if (ei->near2[j]) { meanD2near2 += ei->d2[j]; ++k2; }
            if (k2) meanD2near2 /= k2;
        }

        // Selection from two possible eyeball center locations is a weighted product of:
        // - distance from designed 'perfect' center
        // - distance from previous center detection
        // - mean distance to nearby centers (equivalent of cluster radius) from previous frames
        //   considering avg cluster of 2mm or above very unlikely
        float w1n = fmaxf(0.1f, 1 - 4*sqrtf(meanD1near1)/iris2ctr_pxl);
        float w2n = fmaxf(0.1f, 1 - 4*sqrtf(meanD2near2)/iris2ctr_pxl);

        int *near;
        if (w2d*w2p*w2n > w1d*w1p*w1n)
        {
            near = ei->near2;
            exy[0] = e2v[0]; exy[1] = e2v[1];
        }
        else
        {
            near = ei->near1;
            exy[0] = e1v[0]; exy[1] = e1v[1];
        }

        // initialize tail to the first frame if no history
        if (n_valid<=1)
        {
            ei->ctr_tail[0] = ei->iris2d[0] + (exy[0]-ei->iris2d[0]) * ei->iris2hrot_mm/ei->iris2vrot_mm;
            ei->ctr_tail[1] = exy[1];
        }

        // stabilization of 2d eyeball center by weighted-averaging with last readings
        // Note: this stabilization works best for on-axis camera
        // For off-axis camera the separate x/y stabilization is flawed
        // because eye h/v rotations projected on camera 2d plane are not straight lines anymore
        float meanx=0, meany=0, sumw=0;
        k=0;
        for (j=0; j<n_valid2; ++j)
        {
            if (near[j])
            {
                meanx += ei->obs[j][0] * ei->weight[j];
                meany += ei->obs[j][1] * ei->weight[j];
                sumw += ei->weight[j];
                ++k;
            }
        }
        if (sumw)
        {
            meanx /= sumw;
            meany /= sumw;

            *eyeball_confid = sumw / k;
        }

        // Update long-term observation tail (keeping about 2 sec)
        // readings of higher confidence provide stronger updates to the tail
        float tail_upd = fmaxf(0.1f, 3*(*eyeball_confid * *eyeball_confid));
        ei->ctr_tail[0] = (ei->ctr_tail[0]*(1<<ei->smooth_2dtail) + meanx*tail_upd) / ((1<<ei->smooth_2dtail) + tail_upd);
        ei->ctr_tail[1] = (ei->ctr_tail[1]*(1<<ei->smooth_2dtail) + meany*tail_upd) / ((1<<ei->smooth_2dtail) + tail_upd);
    }

    // derive horizontal center position from position of horizontal rotation axis
    float tail[2];
    tail[0] = ei->iris2d[0] + (ei->ctr_tail[0]-ei->iris2d[0])*ei->iris2vrot_mm / ei->iris2hrot_mm;
    tail[1] = ei->ctr_tail[1];

    // smooth, but no smoothing for high-quality detections
    exy[0] = (tail[0]*((1<<ei->smooth_2deye)-1)*(1-quals) + exy[0]*quals) / (((1<<ei->smooth_2deye)-1)*(1-quals) + quals);
    exy[1] = (tail[1]*((1<<ei->smooth_2deye)-1)*(1-quals) + exy[1]*quals) / (((1<<ei->smooth_2deye)-1)*(1-quals) + quals);
}


void EyeballCenter3d(float iris_z, float old_z, float iris2d[2], float eyeball2d[2], float iris2ctr_pxl, float pxl2mm, float eyeball3d[3])
{
    // Note: here we are working in parallel-projection mode (camera perspective projection is unrelated)

    // squared distance between iris and eyeball center projection on an image plane
    float iris2eyeball2d_sq = norm_euc2(iris2d, eyeball2d);

    // eyeball z distance
    if (iris2eyeball2d_sq > iris2ctr_pxl*iris2ctr_pxl)
        // eyeball center too far may happen due to position smoothing in Estimate2dEyeCenter
        eyeball3d[2] = old_z;
    else
        eyeball3d[2] = pxl2mm * ( iris_z + sqrtf(iris2ctr_pxl * iris2ctr_pxl - iris2eyeball2d_sq) );

    eyeball3d[0] = pxl2mm * eyeball2d[0];
    eyeball3d[1] = pxl2mm * eyeball2d[1];
}


void FilterEyeballPosition(float eyeball3d[3], float eyeball3d0x[3], float iris3d[3], float iris2eye_mm, float iris2hrot_mm, int strength, float old_confid, float eyeball_confid, float eyeball3dstab[3])
{
    // axis for vertical eye rotation is closer to cornea than horizontal rotation axis (which is coinciding with eyeball center)
    // adjust coordinates of stabilized point for that
    float gaze_norm[3];
    gaze_norm[0] = (eyeball3d[0]-iris3d[0]) / iris2eye_mm;
    gaze_norm[1] = (eyeball3d[1]-iris3d[1]) / iris2eye_mm;
    gaze_norm[2] = (eyeball3d[2]-iris3d[2]) / iris2eye_mm;
    float dx = -(iris2hrot_mm-iris2eye_mm) * gaze_norm[0];

    if (old_confid)   // no filtering for the first reading
    {
        float upd_weight = fminf(1, eyeball_confid / old_confid);

        // reduce correction of 3d eyeball center stabilized Z location from high-angles and XY for very low-angles
        // x,y: full correction > 30 degree
        float wx = fmaxf(0.1f, fminf(1, (1-gaze_norm[2])*8));
        float wy = wx;
        // full correction < 30 degree, almost no correction > 45
        float wz = fmaxf(0.1f, fminf(1, (gaze_norm[2]-0.5f)*5));

        eyeball3d0x[0] = ( eyeball3d0x[0]*((1<<strength)-1) + (eyeball3d[0] - dx) * upd_weight*wx ) / ( ((1<<strength)-1) + upd_weight*wx );
        eyeball3d0x[1] = ( eyeball3d0x[1]*((1<<strength)-1) + (eyeball3d[1] - 0)  * upd_weight*wy ) / ( ((1<<strength)-1) + upd_weight*wy );
        eyeball3d0x[2] = ( eyeball3d0x[2]*((1<<strength)-1) + (eyeball3d[2] - 0)  * upd_weight*wz ) / ( ((1<<strength)-1) + upd_weight*wz );

        eyeball3dstab[0] = eyeball3d0x[0] + dx;
        eyeball3dstab[1] = eyeball3d0x[1] + 0;
        eyeball3dstab[2] = eyeball3d0x[2] + 0;
    }
    else
    {
        eyeball3d0x[0] = eyeball3d[0] - dx;
        eyeball3d0x[1] = eyeball3d[1] - 0;
        eyeball3d0x[2] = eyeball3d[2] - 0;

        eyeball3dstab[0] = eyeball3d[0];
        eyeball3dstab[1] = eyeball3d[1];
        eyeball3dstab[2] = eyeball3d[2];
    }
}


// estimate pixel to millimeters scaling factor from
// 3d location of iris center in pixels and 3d location of eyeball center in mm
float Pxl2mmFromEye2Iris(float eye_mm[3], float iris_pxl[3], float iris2eye_mm)
{
    // solve:
    // norm( pxl2mm * iris_pxl - eye_mm ) = iris2eye_mm
    // iris2eye_mm^2 = sum((pxl2mm * iris_pxl - eye_mm).^2)

    //  px=iris_pxl(1); py=iris_pxl(2); pz=iris_pxl(3);
    //  ex=eye_mm(1);    ey=eye_mm(2);    ez=eye_mm(3);
    //
    // iris2eye_mm^2 = (pxl2mm*px-ex)^2 + (pxl2mm*py-ey)^2 + (pxl2mm*pz-ez)^2;
    // iris2eye_mm^2 = (pxl2mm*px)^2 - 2*pxl2mm*px*ex + ex^2  + ..y.. + ..z..;
    // solve quadratic equation:
    // (pxl2mm^2)*(px^2+py^2+pz^2) - 2*pxl2mm*(px*ex+py*ey+pz*ez) + ex^2+ey^2+ez^2-iris2eye_mm^2 = 0

    float a = iris_pxl[0]*iris_pxl[0] + iris_pxl[1]*iris_pxl[1] + iris_pxl[2]*iris_pxl[2];
    float b = -2*(eye_mm[0]*iris_pxl[0] + eye_mm[1]*iris_pxl[1] + eye_mm[2]*iris_pxl[2]);
    float c = eye_mm[0]*eye_mm[0] + eye_mm[1]*eye_mm[1] + eye_mm[2]*eye_mm[2] - iris2eye_mm*iris2eye_mm;

    float D = b*b-4*a*c;

    // if eyeball center too far
    if (D<0) return 0;

    // minus sign below results in selecting smaller pxl2mm of two possible solutions
    // this reflects that the pupil is closer to the camera than the eyeball center
    return (-b-sqrtf(D)) / (2*a);
}


void externalPupilData(int *param, int16_t *edge, int16_t *digest, int *xylen, int *xs, int *ys, float *parea)
{
    int xmin,xmax,ymin,ymax;

    xmin=xmax=edge[0];
    ymin=ymax=edge[1];
    *xylen = param[0];
    for (int i=0; i<*xylen; ++i)
    {
        int16_t x = edge[i*2];
        int16_t y = edge[i*2+1];

        digest[i*4] = x; digest[i*4+1] = y;
        if (x < xmin) xmin=x;
        if (x > xmax) xmax=x;
        if (y < ymin) ymin=y;
        if (y > ymax) ymax=y;
    }

    if (param[1] || param[2])
    {
        *xs = param[1];
        *ys = param[2];
    }
    else
    {
        *xs = (xmax+xmin)/2;
        *ys = (ymax+ymin)/2;
    }

    if (param[3])
        *parea = (float)param[3];
    else
        *parea = F_PI * (xmax-xmin)/2 * (ymax-ymin)/2;
}


// simple bilateral-like filter to reduce gaze direction jitter
void smoothGazeVector(float gaze_vector[3], float gv_smooth[3], float gv_smooth_diff, int fps)
{
    for (int i=0; i<3; ++i)
    {
        float d = fabsf(gv_smooth[i]-gaze_vector[i]);
        float w = fmaxf(0, fminf(fps/4, (1-d/gv_smooth_diff)*fps/4 ));

        gv_smooth[i] = (gv_smooth[i]*w + gaze_vector[i]) / (w+1.0f);
    }

    float gv_z2 = 1 - (gv_smooth[0]*gv_smooth[0] + gv_smooth[1]*gv_smooth[1]);

    if (gv_z2>0)
        gv_smooth[2] = sqrtf(gv_z2);
    else
    {
        gv_smooth[2] = gaze_vector[2];
        // re-normalize smoothed gaze vector
        float norm_egaze_vec = sqrtf(gv_smooth[0]*gv_smooth[0] + gv_smooth[1]*gv_smooth[1] + gv_smooth[2]*gv_smooth[2]);
        gv_smooth[0] /= norm_egaze_vec;
        gv_smooth[1] /= norm_egaze_vec;
        gv_smooth[2] /= norm_egaze_vec;
    }
}


// ToDo:
//   - Need a condition to set ei->iris_r to 0 when no eye is visible for a long time
//   - need to change initEyeConfig api and add ei=setCalibrated(ei, ..) call to separate calibration-related parameters
//
// Eye pose estimation is from a single frame, except for the following pieces where information from multiple combined:
//   - 2d eyeball center location averaging from ei.N recent frames in Estimate2dEyeCenter: controlled with  ei.smooth_2deye, ei.smooth_2dtail
//   - 3d eyeball location filtering in FilterEyeballPosition: ei.smooth_3deye
//   - smoothing of detected iris radius in pixels: ei.smooth_iris
//   - adjustment of user iris radius in mm (adjusted only if smooth_3deye is non-zero): ei.adj_iris_mm
// Set above control coefficients in eye config to 0 to disable correspondent smoothing or adjustment
//
void getEyePose(eye_cfg *ei, uint8_t *im, int *pup_param, int16_t *pup_edge)
{
    float min_iris_r  = (float)ei->min_iris_r;
    float max_iris_r  = (float)ei->max_iris_r;
    float hmn_iris_mm = ei->hmn_iris_mm;
    float iris2eye_mm = ei->iris2eye_mm;
    float iris_z      = ei->iris_z;
    float camera_rot  = ei->camera_rot;

    int w = ei->w;
    int h = ei->h;

    // ------------ Find pupil center and edge location on the camera image

    // minimum pupil radius in pixels: 1.5mm pupil diameter, farthest distance from the camera
    float min_pup_r = min_iris_r*1.5f/hmn_iris_mm;
    
    // maximum pupil radius in pixels: 8mm pupil diameter, closest distance to the camera
    float max_pup_r = max_iris_r*8.0f/hmn_iris_mm;

    // Filtering goals:
    // 1st filter (over full frame):
    // - reduce likelyhood of pupil threshold level being affected by small/thin dark features
    // - brighten outward frame area to avoid false pupil detection in eye corners
    //
    // 2nd filter (area covering pupil):
    // - ideally - remove small bright objects (lashes, glints) while keeping edges intact
    //
    // No filtering needed for iris radius detection stage,
    // except may be for removal of bright glints and lashes aligned with raster vertically
    // (but likely they have no influence due to large integration area)

    // lashes / glints / eyewear reflections removal
    // using tmp as a temporary storages here
    removeSpeckles1(im, ei->tmp, ei->im2x, ei->valid_area, w, h);

    int xylen;
    int xs, ys;
    float parea;
    if (pup_edge && pup_param) // if external pupil detector used
        externalPupilData(pup_param, pup_edge, ei->digest, &xylen, &xs, &ys, &parea);
    else
    {
        uint8_t pupthr = thresholdForPupil(ei->im2x, w, h, min_pup_r, ei->pup_thr, ei->valid_area);

        // find the biggest black area via horizontal and vertical scans
        int psr = getBigBlackArea(ei->im2x, pupthr, ei->digest, ei->dweight, w, h, ei->valid_area, &xs, &ys);

        if ((xs<0) || (ys<0)) { ei->blink=1; return; }   // no pupil found (either eyelids closed or out of valid area)

        // compute filtered/brightened frame at full resolution over small pupils area
        // +4 here to cover neighbor pixels required during extractEdgePixels
        removeSpeckles2(im, ei->tmp, ei->valid_area, w, h, xs, ys, (int)(max_pup_r/2+4));

        // optional - re-threshold pupil using local levels around initial detection
        // useful in contrast-reduced situations, eg reflections from glasses at oblique angle
        reThreshold(im, ei->im2x, ei->tmp, w, h, xs, ys, psr, pupthr, max_pup_r);

        // repurposing digest storage to keep list of [x y] coordinates of pupil edge pixels
        xylen = extractEdgePixels(ei->tmp, w, h, &xs, &ys, psr, (int)max_pup_r, ei->digest, &parea);

        if (xylen < 8) { ei->blink=1; return; }         // require at least 8 pixels to assume the pupil ellipse found
    }


    // ------------ Fit pupil to ellipse

    // clean the pixel coord data by eliminating pixels that are 'internal' to the radially-distant pixels
    xylen = cleanInternals(ei->digest, xylen, parea);

    // stop here and return previously estimated value if no pupil detection or very low confidence in detection
    if (xylen < 8) { ei->blink=1; return; }         // require at least 8 pixels to assume the pupil ellipse found

    // converting xy from int to float here
    for (int i=0; i<xylen; ++i)
    {
        float x=ei->digest[i*4];
        float y=ei->digest[i*4+1];
        ei->xy[i*2] = x;
        ei->xy[i*2+1] = y;
    }

    // smart fit of pupil to ellipse
    float A[6];
    xylen = pupilEllipse(A, ei->xy, xylen, (float*)(ei->tmp), (int)min_pup_r, (int)max_pup_r);

    // stop here and return previously estimated value if no pupil detection
    if (xylen < 8) { ei->blink=1; return; }         // require at least 8 fitted pixels to assume the pupil ellipse found

    // negate pupil non-circularity, if any
    float cxy[2], semx, semy, phi;
    EllipseGeomFromAlgebraic(A, &cxy[0], &cxy[1], &semx, &semy, &phi);
    ei->el_ratio0 = fminf(semx,semy)/fmaxf(semx,semy);       // pupil ellipse ratio before any corrections (used in calibration)
    for (int i=0; i<xylen; ++i)
    {
        ei->xy[i*2] -= cxy[0];
        ei->xy[i*2+1] -= cxy[1];
    }

    compensatePupilElongation(ei->xy, xylen, ei->el_ratio0, ei->pupil_angle*F_PI/180, ei->pupil_stretch);
    int fitfound = EllipseDirectFit(ei->xy, xylen, 0, 0, A);

    // stop here and return previously estimated value if no pupil detection
    if (!fitfound) { ei->blink=1; return; }

    float ctmp[2];
    EllipseGeomFromAlgebraic(A, &ctmp[0], &ctmp[1], &semx, &semy, &phi);

    // final list of pupil ellipse edge pixels, used in calibration
    ei->xylen = xylen;

    cxy[0] += xs;
    cxy[1] += ys;
    // edge detection with ellipse fitting cause undershoot of major axis for small pupils
    // subtracting 1 pixel from each side makes result closer to ground truth without affecting large or high-ratio pupils
    semx = fmaxf(1, semx-2);
    semy = fmaxf(1, semy-2);

    // figure out how much of pupil is visible and set blink ratio from that
    float peli_area = F_PI*semx*semy;
    float blink = fmaxf(0, 1 - 1.1f * parea / peli_area);

    if (ei->smooth_2deye>=0)  // only react to sudden changes if 2d smoothing used
    {
        // detecting sudden change in pupil area, and if such - discard this frame
        // sudden enlargement
        if ((peli_area > 2*ei->peli_area) && (blink > 0.5f))
        {
            ei->peli_area = (15*ei->peli_area + fmaxf(ei->peli_area/2, fminf(ei->peli_area*2, peli_area)))/16;
            return;
        }

        // sudden contraction
        if ((ei->peli_area > 2*peli_area) && (ei->blink < 0.5f))
        {
            ei->peli_area = (ei->peli_area + peli_area)/2;
            return;
        }
    }
    
    // maintain smoothed pupil area in eye instance
    if (ei->blink > 0.8f)
        ei->peli_area = (ei->peli_area + peli_area)/2;
    else
        ei->peli_area = (15*ei->peli_area + fmaxf(ei->peli_area/2, fminf(ei->peli_area*2, peli_area)))/16;

    ei->blink = blink;


    // ------------ Correct pupil position and minor axis due to refraction

    // eyeball center needed to disambiguate position correction direction, but center only dicovered after iris detection
    // so, use eyeball center, iris radius and iris ratio prediction from the previous frame
    float eyeball2d_unrot[2];
    rotcam(ei->eyeball2d, -camera_rot, eyeball2d_unrot);
    eyeball2d_unrot[0] += w/2; eyeball2d_unrot[1] += h/2;

    float el_major, el_minor, phi_minor;
    float pc[2] = {cxy[0], cxy[1]};
    float ixy[2] = {cxy[0], cxy[1]};
    TruePupilFromRefraction(semx, semy, phi, ei->pxl2mm, eyeball2d_unrot, ei->pupil_refr, ei->pca_angle, ei->iris_ratio, iris_z, w/2, h/2, iris2eye_mm-ei->pupil2eye_mm,
        pc, &el_major, &el_minor, &phi_minor);

    if (ei->el_known) // if calibrating - use supplied ground truth phi_minor and el_ratio
    {
        phi_minor = ei->phi_minor;
        el_minor = el_major * ei->el_ratio;
    }
    else
    {
        ei->phi_minor = phi_minor; // used in calibration
        // ellipse main axis is the radius of the pupil/iris (it is perpindicular to the camera ray pointing to the iris center)
        ei->el_ratio  = el_minor/el_major;
    }

    // initial estimation on iris center
    correctPupilDecenter(pc, ei->el_ratio, phi_minor, ei->pxl2mm, ei->pupil_decenter_mm, camera_rot, ixy);

    pc[0] -= w/2; pc[1] -= h/2;
    rotcam(pc, camera_rot, ei->pupil2d);    // used in calibration


    // ------------ Rasterize iris

    ei->iris2d_prev[0] = ei->iris2d[0];
    ei->iris2d_prev[1] = ei->iris2d[1];
    float ixyctr[2];
    ixyctr[0] = ixy[0]-w/2;
    ixyctr[1] = ixy[1]-h/2;
    rotcam(ixyctr, camera_rot, ei->iris2d_p);
    float rotrad = camera_rot * F_PI/180;

    // fill the list of radial 'scalers' computed from ellipse shape which will turn iris shape into circle
    for (int i=0; i<IRIM_H; ++i)
    {
        ei->iris_scl1[i] = ellipseScale(ei->iris_d1[i], phi_minor, ei->el_ratio);
        ei->iris_scl2[i] = ellipseScale(ei->iris_d2[i], phi_minor, ei->el_ratio);
    }

    // rasterize iris sectors in polar coordinates
    int ir_r1=0, ir_r2=0;
    int peak1=0, peak2=0;
    float decent_v1=0, decent_v2=0;

    // at horizontal iris angles above roughly 40 degree (iris_scl=0.78) the iris edge farthest from the camera is occluded by corneal bulge and should not be used
    // at such angles the farther edge is also significantly farther in Z distance and iris radius readings would require correction according to camera ray angle
    if ((ei->iris2d_p[0] >= ei->eyeball_ctr[0]) || (ei->iris_scl1[IRIM_MIDH] > 0.79f))
    {
        rasterizeIris(ei->iris_raster, ei->irim_w, ei->iris_d1, ei->iris_scl1, ei->im2x, w, h, ixy, min_iris_r, iris_z, phi_minor, camera_rot);
        detectIrisRadius(ei->iris_raster, ei->tmp, ei->irim_w, IRIM_H, min_iris_r, el_major, &ir_r1, &peak1, &decent_v1);
    }
    if ((ei->iris2d_p[0] <= ei->eyeball_ctr[0]) || (ei->iris_scl2[IRIM_MIDH] > 0.79f))
    {
        rasterizeIris(ei->iris_raster, ei->irim_w, ei->iris_d2, ei->iris_scl2, ei->im2x, w, h, ixy, min_iris_r, iris_z, phi_minor, camera_rot);
        detectIrisRadius(ei->iris_raster, ei->tmp, ei->irim_w, IRIM_H, min_iris_r, el_major, &ir_r2, &peak2, &decent_v2);
    }

    // ToDo: discard detection (set eyeball_confidence=0 and return data from previous frame) if:
    // - new iris center is far from previous and detected peaks are significantly below ones for previous frame
    // - new iris radius is far from previous frame and eyeball_confidence was good for previous frame


    // ------------ Find the iris radius and center, in camera pixels

    // starting from here, switching from frame pixel coordinates to coordinate system with matched horizontal/vertical eye directions
    phi_minor += rotrad;
    
    // get combined iris radius from two measurements
    // adjusting refraction-/decenter-corrected iris2d and phi_minor direction from iris edege curves
    getIris(ei, (float)ir_r1, (float)ir_r2, (float)peak1, (float)peak2, ei->iris2d_p, ei->iris_recenter, decent_v1, decent_v2, ei->el_ratio,
        &ei->iris_r, &ei->iris_r1, &ei->iris_r2, &ei->iris_peak, ei->iris2d, &phi_minor);

    // used in calibration
    rotcam(ei->iris2d, -camera_rot, ei->cxy);
    ei->cxy[0] += w/2; ei->cxy[1] += h/2;


    // ------------ From above - 3d iris center location in camera space and pixels to millimeters ratio

    // if iris not detected - do not compute iris 3d position, use previous result as a stub
    if (ei->iris_r)
        ei->pxl2mm = hmn_iris_mm/2/ei->iris_r;    // pixels to millimeters ratio

    float iris2eye_pxl = iris2eye_mm/ei->pxl2mm;
    ei->iris3d_cam[0] = ei->pxl2mm * ei->iris2d[0];
    ei->iris3d_cam[1] = ei->pxl2mm * ei->iris2d[1];
    ei->iris3d_cam[2] = ei->pxl2mm * iris_z;


    // ------------ From (mulitple) ellipses - estimate center around which they are moving (= 2d eyeball center)

    float eyeball_confid = 1;
    if (ei->smooth_2deye < 0)     // negative smoothing means calibration pass, and eyeball center is fixed and given
    {
        // derive center position from position of rotation axes
        ei->eyeball2d[0] = ei->iris2d[0]+(ei->ctr_tail[0]-ei->iris2d[0])*ei->iris2vrot_mm/ei->iris2hrot_mm;
        ei->eyeball2d[1] = ei->ctr_tail[1];
    }
    else
    {
        if (ei->frame % ei->fps30)  // no need to estimate eyeball center more frequently than 30 times a second
            eyeball_confid = ei->eyeball_confidence;
        else
        {
            float exy[2];
            Estimate2dEyeCenter(ei, ei->el_ratio, phi_minor, ei->eyeball2d, iris2eye_pxl, exy, &eyeball_confid);
            ei->eyeball2d[0] = exy[0]; ei->eyeball2d[1] = exy[1];
        }
    }

    // ------------ From 2d eyeball center - 3d eyeball center location in camera 3d space.

    // use previous location if it will be impossible to compute new one below
    float stabilized_eyeball3d_cam[3];
    stabilized_eyeball3d_cam[0] = ei->eyeball3d_cam[0];
    stabilized_eyeball3d_cam[1] = ei->eyeball3d_cam[1];
    stabilized_eyeball3d_cam[2] = ei->eyeball3d_cam[2];

    if (ei->iris_r && eyeball_confid>0)
    {
        float eyeball3d_cam[3];
        EyeballCenter3d(iris_z, ei->eyeball3d_cam[2], ei->iris2d, ei->eyeball2d, iris2eye_pxl, ei->pxl2mm, eyeball3d_cam);
		// debug: disable center stabilization
        // stabilized_eyeball3d_cam[0] = eyeball3d_cam[0]; stabilized_eyeball3d_cam[1] = eyeball3d_cam[1]; stabilized_eyeball3d_cam[2] = eyeball3d_cam[2];
        FilterEyeballPosition(eyeball3d_cam, ei->eyeball3d0x, ei->iris3d_cam, iris2eye_mm, ei->iris2hrot_mm, ei->smooth_3deye, ei->eyeball_confidence, eyeball_confid, stabilized_eyeball3d_cam);
    }


    // ------------ Recompute 3d iris Z from stabilized 3d eyeball center and 2d iris location

    float iris2d_mm[2] = {ei->pxl2mm*ei->iris2d[0], ei->pxl2mm*ei->iris2d[1]};
    float iris2eye_2d = sqrtf(norm_euc2(iris2d_mm, stabilized_eyeball3d_cam));

    if (iris2eye_mm <= iris2eye_2d)
    {
        // if impossible to compute 3d iris location from stabilized eyeball center (too far) -
        // scale direction towards stabilized center to position eyeball center at correct distance from iris
        float eyeball_vec[3];
        eyeball_vec[0] = stabilized_eyeball3d_cam[0] - ei->iris3d_cam[0];
        eyeball_vec[1] = stabilized_eyeball3d_cam[1] - ei->iris3d_cam[1];
        eyeball_vec[2] = stabilized_eyeball3d_cam[2] - ei->iris3d_cam[2];
        float norm_eyeball_vec = sqrtf(eyeball_vec[0]*eyeball_vec[0] + eyeball_vec[1]*eyeball_vec[1] + eyeball_vec[2]*eyeball_vec[2]);
        stabilized_eyeball3d_cam[0] = ei->iris3d_cam[0] + eyeball_vec[0] * (iris2eye_mm/norm_eyeball_vec);
        stabilized_eyeball3d_cam[1] = ei->iris3d_cam[1] + eyeball_vec[1] * (iris2eye_mm/norm_eyeball_vec);
        stabilized_eyeball3d_cam[2] = ei->iris3d_cam[2] + eyeball_vec[2] * (iris2eye_mm/norm_eyeball_vec);
        iris2eye_2d = fminf( iris2eye_mm, sqrtf(norm_euc2(iris2d_mm, stabilized_eyeball3d_cam)) );
    }
    float iris_dz = sqrtf(iris2eye_mm*iris2eye_mm - iris2eye_2d*iris2eye_2d);
    ei->iris3d_cam[2] = stabilized_eyeball3d_cam[2]-iris_dz;

    ei->eyeball3d_cam[0] = stabilized_eyeball3d_cam[0];  ei->eyeball3d_cam[1] = stabilized_eyeball3d_cam[1];  ei->eyeball3d_cam[2] = stabilized_eyeball3d_cam[2];
    ei->eyeball_confidence = eyeball_confid;

    // estimate stabilized iris ellipse ratio (used in TruePupilFromRefraction() for the next frame)
    float iris_ratio = (ei->eyeball3d_cam[2]-ei->iris3d_cam[2]) / iris2eye_mm;
    cxy[0] = cxy[0]-w/2; cxy[1] = cxy[1]-h/2;
    if (ei->pxl2mm*ei->pxl2mm * norm_euc2(cxy, ei->iris2d_s) > 2*2)   // less smoothing if iris has moved far away (2mm displacement, possibly accumulated over multiple frames)
    {
        ei->iris_ratio = (ei->iris_ratio + iris_ratio) / 2;
        ei->iris2d_s[0] = (ei->iris2d_s[0] + cxy[0]) / 2;
        ei->iris2d_s[1] = (ei->iris2d_s[1] + cxy[1]) / 2;
    }
    else
    {
        ei->iris_ratio = ( 15*ei->iris_ratio + iris_ratio ) / 16;
        ei->iris2d_s[0] = (15*ei->iris2d_s[0] + cxy[0]) / 16;
        ei->iris2d_s[1] = (15*ei->iris2d_s[1] + cxy[1]) / 16;
    }


    // ------------ Linear transformation of pupil/eyeball position into HMD coordinate space.

    ei->iris3d[0] = ei->cam2hmd[0][0]*ei->iris3d_cam[0] + ei->cam2hmd[0][1]*ei->iris3d_cam[1] + ei->cam2hmd[0][2]*ei->iris3d_cam[2] + ei->cam2hmd[0][3];
    ei->iris3d[1] = ei->cam2hmd[1][0]*ei->iris3d_cam[0] + ei->cam2hmd[1][1]*ei->iris3d_cam[1] + ei->cam2hmd[1][2]*ei->iris3d_cam[2] + ei->cam2hmd[1][3];
    ei->iris3d[2] = ei->cam2hmd[2][0]*ei->iris3d_cam[0] + ei->cam2hmd[2][1]*ei->iris3d_cam[1] + ei->cam2hmd[2][2]*ei->iris3d_cam[2] + ei->cam2hmd[2][3];

    ei->eyeball3d[0] = ei->cam2hmd[0][0]*ei->eyeball3d_cam[0] + ei->cam2hmd[0][1]*ei->eyeball3d_cam[1] + ei->cam2hmd[0][2]*ei->eyeball3d_cam[2] + ei->cam2hmd[0][3];
    ei->eyeball3d[1] = ei->cam2hmd[1][0]*ei->eyeball3d_cam[0] + ei->cam2hmd[1][1]*ei->eyeball3d_cam[1] + ei->cam2hmd[1][2]*ei->eyeball3d_cam[2] + ei->cam2hmd[1][3];
    ei->eyeball3d[2] = ei->cam2hmd[2][0]*ei->eyeball3d_cam[0] + ei->cam2hmd[2][1]*ei->eyeball3d_cam[1] + ei->cam2hmd[2][2]*ei->eyeball3d_cam[2] + ei->cam2hmd[2][3];


    // ------------ Gaze direction calculated from eyeball center vs iris center 3d locations.

    float optaxis_vector[3];
    optaxis_vector[0] = (ei->eyeball3d[0]-ei->iris3d[0]) / iris2eye_mm;
    optaxis_vector[1] = (ei->eyeball3d[1]-ei->iris3d[1]) / iris2eye_mm;
    optaxis_vector[2] = (ei->eyeball3d[2]-ei->iris3d[2]) / iris2eye_mm;
    // rotate and possibly stretch gaze direction by alpha angle matrix
    ei->gaze_vector[0] = ei->alpha_mat[0][0]*optaxis_vector[0] + ei->alpha_mat[0][1]*optaxis_vector[1] + ei->alpha_mat[0][2]*optaxis_vector[2];
    ei->gaze_vector[1] = ei->alpha_mat[1][0]*optaxis_vector[0] + ei->alpha_mat[1][1]*optaxis_vector[1] + ei->alpha_mat[1][2]*optaxis_vector[2];
    ei->gaze_vector[2] = ei->alpha_mat[2][0]*optaxis_vector[0] + ei->alpha_mat[2][1]*optaxis_vector[1] + ei->alpha_mat[2][2]*optaxis_vector[2];
    // simple filter to reduce gaze direction jitter (less filtering if angle difference is above 6 degree)
    smoothGazeVector(ei->gaze_vector, ei->gaze_vector_smooth, ei->gv_smooth_diff, ei->fps);

    // record angle from HMD center in observations
    if (ei->last_bin>=0)
        ei->ang_obs[ei->last_bin] = atan2f(ei->gaze_vector[1], ei->gaze_vector[0]);
    ++ei->frame;
}
