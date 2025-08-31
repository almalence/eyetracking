/* ------------------------------------------------------------------------- *\

    Structure holding eye config and state; API functions declarations
 
    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
    SPDX-License-Identifier: GPL-3.0-only
    All rights reserved.

    Author: Dmitry Shmunk

\* ------------------------------------------------------------------------- */

#pragma once

#include <stdint.h>

#ifdef  __cplusplus
extern "C" {
#endif

// -------- Eye calibration errors
#define EYECAL_ALL_OK        0
#define EYECAL_NO_MEMORY     1
#define EYECAL_NO_PUPIL      2      // not enough frames with visible pupil
#define EYECAL_AXES_RANGE    3      // eye rotational axes are out of range
#define EYECAL_MAX_ERR       4      // total count of error codes

// keep 64 recent positions to re-estimate eyeball center position from
#define N_EYE_OBS_LOG2       6
#define N_EYE_OBSERVATIONS   (1<<N_EYE_OBS_LOG2)

#define IRIM_H               80     // iris raster height

typedef struct {
    int time;
    int bin;
} timedbins;


// -------- Eye instance struct

typedef struct {

    // -------- various arrays

    // temporary buffers for operations with camera frame
    uint8_t * tmp;
    uint8_t * im2x;
    int16_t * digest;
    int16_t * dweight;
    uint8_t * iris_raster;

    // direction angles and scale for iris ellipse scan (left and right side)
    float iris_scl1[IRIM_H], iris_scl2[IRIM_H];
    float iris_d1[IRIM_H], iris_d2[IRIM_H];

    // -------- Camera-related constants

    int fps;                                    // camera framerate (frames per second)
    int w, h;                                   // camera frame dimension (width, height)
    int irim_w;                                 // iris scan width
    int min_iris_r;                             // minimum and maximum iris radius, in pixels
    int max_iris_r;
    float iris_z;                               // distance to image plane, pixels

    int fps30;                                  // how many frames per 1/30th of a second
    int fps_ec2d;                               // fps at which 2d eyeball center estimations are performed
    uint32_t frame;

    // valid_area: Array describing valid area to constrain search for pupil in a form of:
    //   [ctr_x ctr_y horz_dist vert_dist diag_dist corner_brighten] (distances from ctr)
    //   specify diag_dist=horz_dist+vert_dist to have rectangular crop area
    //   specify horz_dist=vert_dist=diag_dist to have diamond-shaped crop area
    //   specify diag_dist = 0.7 * (horz_dist+vert_dist) to have nearly round crop area
    //   set corner_brighten to 0 to assume equally-bright image from the camera
    //   set corner_brighten to 128 to bright-up image corners by factor 1.5
    //   set corner_brighten to 256 to bright-up image corners by factor 2.0, etc (factor = 1+corner_brighten/256).
    int   valid_area[6];

    int   pup_thr;

    float camera_rot;                           // how much to rotate (degree clockwise) camera frame to make eye nearly-horizontal
    float cam2hmd[3][4];                        // Camera (after camera_rot applied) to HMD rotation/translation matrix (constructed from cam2hmd[6] vector)
    float hmd2cam[3][4];                        // Inverse (HMD to Camera) rotation matrix, used in calibration

    // --------  Default temporal filtering strength

    // - zero means no adjustment or filtering
    // - the higher the number - the stronger the smoothing (proportional to 2^n)
    int smooth_2deye;    // 2d eyeball coordinate
    int smooth_3deye;    // 3d eyeball coordinate
    int smooth_2dtail;   // 2d eyeball coordinate tail
    int smooth_iris;     // iris radius in pixels

    float gv_smooth_diff;   // apply smoothing to gaze vector if readings are within ~2 degree

    // -------- Anatomic constants

    // expected eyeball center location in camera space (relative to frame center)
    float eyeball_ctr[2];

    // Iris to rotation axes
    // 15mm and 12mm from cornea apex; minus roughly 3mm from cornea apex to plane where iris edge is
    // A.Ohlendorf etc. 2022. Positions of the horizontal and vertical centre of rotation in eyes with different refractive errors    
    // Horizontal rotation has axis farther from the iris than the true eyeball center
    // Vertical rotation has axis at about the true eyeball center
    // D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 10
    float iris2hrot_mm;
    float iris2vrot_mm;

    // Iris to eyeball center
    // 7.92 mm 
    // Model of the human eye developed for MCNPX (model A) based on the dimensions provided in
    // NCRP Report no 130 (NCRP 1999) and the eye model from Charles and Brown (1975)
    // https://www.researchgate.net/figure/Model-of-the-human-eye-developed-for-MCNPX-model-A-based-on-the-dimensions-provided-in_fig2_51703220
    // Estimate2dEyeCenter and FilterEyeballPosition assume that eyeball center coincide with vertical rotation center
    float iris2eye_mm;

    // Pupil to eyeball center
    // more or less the same as iris-to-eye
    float pupil2eye_mm;

    // Matrix for angle alpha (between optical and visual axes)
    // typical values: 4 degree nasal, 3 degree upwards
    // D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 64
    float alpha_mat[3][3];

    // -------- nearly-constants (may be adjusted from detection data)

    // Iris diameter
    float hmn_iris_mm;

    float pupil_decenter_mm[2];                 // horizontal decenter in mm (>0 - to the right)

    float pupil_refr;                           // pupil widening due to refraction
    float pupil_angle;                          // pupil elongation (non-circularity)
    float pupil_stretch;
    float pca_angle[2];                         // pupillary circular axis horizontal/vertical angles

    // rate at which HMD slippage (eyeball center shift) may occur [0..1]
    // may need to be set higher for synthetic datasets where center shifts abruptly
    float slippage_rate;

    // -------- state data

    float pxl2mm;                               // pixels-to-mm ratio

    float pupil2d[2];                           // pupil location (refraction-corrected) on image plane
    float iris2d_s[2];                          // smoothed iris center location on image plane
    float iris2d_p[2];                          // rough iris center location on image plane, estimated from pupil center
    float iris2d[2];                            // final estimated iris center location on image plane (from iris ellipse)
    float iris2d_prev[2];                       // same for the previous frame
    int   iris_recenter;                        // re-estimate iris center from iris edges: 0=off 1=horizontal only 2=both 3=both + ellipse angle(phi_minor)
    float eyeball2d[2];                         // estimated eye center location on image plane

    float iris3d_cam[3];                        // estimated 3d iris and eye center locations in camera space
    float eyeball3d_cam[3];
    float eyeball3d0x[3];                       // stabilized eye center in camera space, assuming 0-degree vertical rotation
    float iris3d[3];                            // estimated 3d iris and eye center locations in hmd space
    float eyeball3d[3];
    
    float eyeball_confidence;

    float gaze_vector[3];
    float gaze_vector_smooth[3];                // gaze vector with jitter smoothed at fixation positions

    float blink;                                // proportional to amount of pupil occlusion (0=fully visible, 1=eye fully closed)
    float peli_area;                            // area of pupil ellipse (pixels)

    float iris_r;                               // currently detected iris radius in pixels
    float iris_peak;                            // iris detection strength
    float iris_ratio;                           // stabilized ratio of iris ellipse

    // values used during calibration
    float *xy;                                  // list of pupil edge points
    int   xylen;
    float cxy[2];                               // pupil center in camera frame
    float el_ratio0;                            // pupil ellipse ratio before any corrections
    float el_ratio;                             // pupil ellipse ratio (after refraction/projection/etc corrections)
    float phi_minor;                            // minor ellipse axis direction (after refraction/projection/etc corrections)
    float iris_r1;                              // iris radius in pixels, estimation from the left edge
    float iris_r2;                              // iris radius in pixels, estimation from the right edge
    float el_known;                             // flag: do not estimate iris ellipse ratio from image, true ratio provided

    // -------- per-observation data

    int obs_time[N_EYE_OBSERVATIONS];
    int valid[N_EYE_OBSERVATIONS];
    float phi_obs[N_EYE_OBSERVATIONS];
    float ang_obs[N_EYE_OBSERVATIONS];
    float ctr_obs[N_EYE_OBSERVATIONS][2][2];
    float iris2d_obs[N_EYE_OBSERVATIONS][2];
    float ctr_weight[N_EYE_OBSERVATIONS][2];
    float ctr_tail[2];
    int   last_bin;

    // per-observation data internal to Estimate2dEyeCenter(), moved here to reduce pressure on stack
    timedbins timesorted_bins[N_EYE_OBSERVATIONS];
    int valid_bins[N_EYE_OBSERVATIONS];
    float obs[N_EYE_OBSERVATIONS*2][2];
    float weight[N_EYE_OBSERVATIONS*2];
    float d1[N_EYE_OBSERVATIONS*2], d2[N_EYE_OBSERVATIONS*2];
    int near1[N_EYE_OBSERVATIONS*2], near2[N_EYE_OBSERVATIONS*2];

} eye_cfg;


// -------- Eye calibration struct

typedef struct {
    float iris2hrot;
    float iris2vrot;
    float iris2pupil;        // how much pupil is closer to the eyeball center than the iris plane
    float pupil_decenter[2];
    float eyeball_ctr[2];
    float pca_angle[2];
    float pupil_refr;
    float pupil_angle;
    float pupil_stretch;
    float gaze_alpha[2];
} eye_calib;


// -------- API functions - eye tracking

// Initialize eye pose estimator instance (eye config)
// 
// Parameters:
//
// w, h        - camera frame width and height
// cam2hmd     - Camera-to-HMD coordinate space (X-Y-Z Euler rotations used, default in eg Blender):
//               - x/y/z rotations (degree, around horizontal axis, around vertical axis, around axis perpendicular to the sensor)
//               - followed by x/y/z displacements (mm)
// gaze_alpha  - angle alpha (between optical and visual axes)
//               horizontal angle, then vertical; for the left eye coordinate system is starting at top-left corner (looking at HMD screen)
//
// See descriptions in structure above for info about the rest of the parameters
// Returns NULL if not enough memory
eye_cfg * initEyeConfig(
    int w,
    int h,
    int min_iris_r,
    int max_iris_r,
    int valid_area[6],
    int pup_thr,
    float camera_hfov,
    float camera_rot,
    int fps,
    float cam2hmd[6]
    );


// Aaply calibrated eye parameters
//
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
void applyCalibration(eye_cfg *ei, eye_calib *ec);


// Release/deallocate all resources
void releaseEyeConfig(eye_cfg *ei);


// Estimate eye pose
//
// Eye pose estimation is from a single frame, except for the following pieces where information from multiple combined:
//   - 2d eyeball center location averaging from N recent frames in Estimate2dEyeCenter: controlled with  ei->smooth_2deye, ei->smooth_2dtail
//   - 3d eyeball location filtering in FilterEyeballPosition: ei->smooth_3deye
//   - smoothing of detected iris radius in pixels: ei->smooth_iris
// Set above control coefficients in eye config to 0 to disable correspondent smoothing or adjustment
//
// Parameters:
//  ei        - pointer to eye instance
//  im        - camera frame (single channel, usually IR)
//  pup_edge  - if non-zero - use externally detected list of pupil edge points (pairs of x,y coordinates)
//  pup_param - array with externally detected pupil parameters: [npoints ctr_x ctr_y area]:
//    npoints - number of edge points in pup_edge (at least 8 should be provided)
//    ctr_x,ctr_y - approximate coordinates of ellpise center (optional, 0 can be passed if unknown)
//    area    - pupil area (optional, 0 can be passed if unknown)
//
// Results of estimation are in:
// ei->iris3d          - iris center location in 3d space (in mm)
// ei->eyeball3d       - eye ball center location in 3d space (in mm)
// ei->gaze_vector     - normalized gaze vector
void getEyePose(eye_cfg *ei, uint8_t *im, int *pup_param, int16_t *pup_edge);


// -------- API functions - eye calibration

// Perform eye calibration
//
// M input frames should contain eye views equally distributed between anchor locations (approximately M/5 frames for each)
// Anchor locations should be in the following order: Left of center, Right of center, Above center, Below center, Center
//
// Parameters:
//  ei        - pointer to eye instance, pre-initialized with initEyeConfig()
//  frames    - array of pointers to frames captured from eye-tracking camera
//  M         - total number of frames
//  anc       - location of anchor points used for calibration, horizontal and vertical position in degree
//
//  Output:
//  err       - EYECAL_ALL_OK if no errors, or error code
//
//  Return:   structure with calibrated eye parameters, to be used with applyCalibration()
eye_calib * calibrateEye(
    eye_cfg   *ei,
    uint8_t  **frames,
    int        M,
    float      anc[5][2],
    int *err
    );


// Free memory allocated for calibration parameters in calibrateEye()
void releaseCalibration(eye_calib * ec);


#ifdef  __cplusplus
}
#endif
