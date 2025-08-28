% -------------------------------------------------------------------------
%
%    Eye-tracking instance initialization
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% Notes on input parameters:
%
% valid_area: Array describing valid area to constrain search for pupil in a form of:
%   [ctr_x ctr_y horz_dist vert_dist diag_dist corner_brighten] (distances from ctr)
%   specify diag_dist=horz_dist+vert_dist to have rectangular crop area
%   specify horz_dist=vert_dist=diag_dist to have diamond-shaped crop area
%   specify diag_dist = 0.7 * (horz_dist+vert_dist) to have nearly round crop area
%   set corner_brighten to 0 to assume equally-bright image from the camera
%   set corner_brighten to 128 to bright-up image corners by factor 1.5
%   set corner_brighten to 256 to bright-up image corners by factor 2.0, etc (factor = 1+corner_brighten/256).
%

function ei = initEyeConfig(w, h, min_iris_r, max_iris_r, valid_area, pup_thr, camera_hfov, camera_rot, fps, cam2hmd)

    % keep recent positions to re-estimate eyeball center position from
    N = 64;
    ei.N = N;

    % -------- Camera-related constants

    ei.fps         = fps;
    ei.min_iris_r  = min_iris_r;                      % minimum and maximum iris radius, in pixels
    ei.max_iris_r  = max_iris_r;
    ei.valid_area  = valid_area;
    ei.pup_thr     = pup_thr;
    ei.camera_rot  = camera_rot;                      % how much to rotate (clockwise) camera frame to make eye nearly-horizontal
    ei.iris_z = (w/2)/tan(camera_hfov/2*pi/180);      % z-axis distance between camera and pupil (or image plane) in pixels

    ei.fps30       = ceil(fps/30);                    % how many frames per 1/30th of a second
    ei.fps_ec2d    = round(fps/ei.fps30);             % fps at which 2d eyeball center estimations are performed
    ei.frame       = uint32(0);

    % Camera pose matrix from XYZ Euler rotations and offset
    cam_xyzrot = rotationMat(cam2hmd(1:3), 0);
    ei.cam2hmd     = [cam_xyzrot cam2hmd(4:6)'];     % Camera (after camera_rot applied) to HMD rotation matrix
    cam_xyzrot = rotationMat(cam2hmd(1:3), 1);
    ei.hmd2cam     = [cam_xyzrot -cam2hmd(4:6)'];    % Inverse (HMD to Camera) rotation matrix, used in calibration

    
    % --------  Default temporal filtering strength

    % - zero means no adjustment or filtering
    % - the higher the number - the stronger the smoothing (proportional to 2^n)
    f = floor(log2(ei.fps_ec2d/30));
    ei.smooth_2deye  = max(1, 3+f);   % 2d eyeball coordinate
    ei.smooth_2dtail = max(2, 4+f);   % 2d eyeball coordinate tail (about 1-2 secs)
    f = floor(log2(fps/30));
    ei.smooth_3deye  = max(2, 3+f);   % 3d eyeball coordinate
    ei.smooth_iris   = max(1, 5+f);   % iris radius in pixels

    ei.gv_smooth_diff = 0.035;          % apply smoothing to gaze vector if readings are within ~2 degree (asin(0.035))

    % -------- Anatomic constants (may be adjusted in calibration)

    % expected eyeball center location in camera space (relative to frame center)
    ei.eyeball_ctr = [0 0];

    % Iris to rotation axes
    % 15mm and 12mm from cornea apex; minus roughly 3mm from cornea apex to plane where iris edge is
    % A.Ohlendorf etc. 2022. Positions of the horizontal and vertical centre of rotation in eyes with different refractive errors    
    % Horizontal rotation has axis farther from the iris than the true eyeball center
    % Vertical rotation has axis at about the true eyeball center
    % D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 10
    ei.iris2hrot_mm = 11.5;
    ei.iris2vrot_mm = 9.4;

    % Iris to eyeball center
    % 7.92 mm 
    % Model of the human eye developed for MCNPX (model A) based on the dimensions provided in
    % NCRP Report no 130 (NCRP 1999) and the eye model from Charles and Brown (1975)
    % https://www.researchgate.net/figure/Model-of-the-human-eye-developed-for-MCNPX-model-A-based-on-the-dimensions-provided-in_fig2_51703220
    % Estimate2dEyeCenter() and FilterEyeballPosition() assume that eyeball center coincide with vertical rotation center
    ei.iris2eye_mm = ei.iris2vrot_mm;

    % Pupil to eyeball center, more or less the same as iris-to-eye, but:
    % - moves back and forth with accommodation by up to 0.4mm (accommodated lens push entrance pupil forward)
    % - can be different in synthetic eyes
    ei.pupil2eye_mm = ei.iris2eye_mm-0.3;

    % Matrix for angle alpha (between optical and visual axes)
    % typical values: 4 degree nasal, 3 degree upwards
    % D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 64
    % horizontal angle, then vertical; for the left eye coordinate system is starting at top-left corner (looking at HMD screen)
    ei.alpha_mat = rotationMat([4 2 0], 0);

    % Iris diameter
    % 12 mm
    % H. Gross. 2008. Handbook of Optical Systems: Vol. 4 Survey of Optical Instruments.
    % Wiley-VCH Verlag GmbH+Co. KGaA, Weinheim.
    %
    % 11.6mm
    % Guillaume Francois, Pascal Gautron, Gaspard Breton, and Kadi Bouatouch
    % Image-Based Modeling of the Human Eye
    % IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 15, NO. 5, SEPTEMBER/OCTOBER 2009
    ei.hmn_iris_mm = 11.6;

    ei.pupil_decenter_mm = [0 0];               % pupil decenter in mm (>0 - to the right/down, looking from the camera)

    ei.pupil_refr = 1.121;                      % pupil widening due to refraction
    ei.pupil_angle = 0;                         % pupil elongation (non-circularity)
    ei.pupil_stretch = 1;
    ei.pca_angle = [0 0];                       % pupillary circular axis horizontal/vertical angles

    % rate at which HMD slippage (eyeball center shift) may occur [0..1]
    % may need to be set higher for synthetic datasets where center shifts abruptly
    ei.slippage_rate = 0.3;

    % -------- state data

    ei.max_el_ratio = 0.9995;                   % <1 to help with phi direction in 1st frame

    ei.pxl2mm = ei.hmn_iris_mm/(min_iris_r+max_iris_r);  % pixels-to-mm ratio

    ei.pupil2d = [0 0];                         % pupil location (refraction-corrected) on image plane
    ei.iris2d_s = [0 0];                        % smoothed iris center location on image plane
    ei.iris2d_p = [0 0];                        % rough iris center location on image plane, estimated from pupil center
    ei.iris2d = [0 0];                          % final estimated iris center location on image plane (from iris ellipse)
    ei.iris2d_prev = [0 0];                     % same for the previous frame
    ei.iris_recenter = 1;                       % re-estimate iris center from iris edges: 0=off 1=horizontal only 2=both 3=both + ellipse angle(phi_minor)
    % estimated eye center location on image plane, for the first frame assume eyeball center is in designed position
    ei.eyeball2d = ei.eyeball_ctr;

    ei.iris3d_cam = [0 0 ei.iris_z*ei.pxl2mm];  % estimated 3d iris and eye center locations in camera space
    ei.eyeball3d_cam = ei.iris3d_cam + [0 0 ei.iris2eye_mm];
    ei.eyeball3d0x = ei.eyeball3d_cam;          % stabilized eye center in camera space, assuming 0-degree vertical rotation
    ei.iris3d = [0; 0; ei.iris2eye_mm];         % estimated 3d iris and eye center locations in hmd space
    ei.eyeball3d = [0 0 0];
    
    ei.eyeball_confidence = 0;

    ei.gaze_vector = [0; 0; 1];
    ei.gaze_vector_smooth = [0; 0; 1];          % gaze vector with jitter smoothed at fixation positions

    ei.blink = 1;                               % assume fully-closed eye initially
    ei.peli_area = 0;                           % area of pupil ellipse (pixels)

    ei.iris_r = 0;                              % currently detected iris radius in pixels
    ei.iris_peak = 0;                           % iris detection strength
    ei.iris_ratio = ei.max_el_ratio;            % stabilized ratio of iris ellipse

    % values used during calibration
    ei.pxy = [];                                % list of pupil edge points
    ei.cxy = [w/2 h/2];                         % pupil center in camera frame
    ei.el_ratio0 = ei.max_el_ratio;             % pupil ellipse ratio before any corrections
    ei.el_ratio = ei.max_el_ratio;              % pupil ellipse ratio (after refraction/projection/etc corrections)
    ei.phi_minor = 0;                           % minor ellipse axis direction (after refraction/projection/etc corrections)
    ei.iris_r1 = 0;                             % iris radius in pixels, estimation from the left edge
    ei.iris_r2 = 0;                             % iris radius in pixels, estimation from the right edge
    ei.el_known = false;                        % flag: do not estimate iris ellipse ratio from image, true ratio provided

    % -------- per-observation data

    ei.obs_time       = zeros(N,1,'single');
    ei.valid          = zeros(N,1,'logical');
    ei.phi_obs        = zeros(N,1,'single');
    ei.ang_obs        = zeros(N,1,'single');
    ei.ctr_obs        = zeros(N,2,2,'single');
    ei.ctr_weight     = zeros(N,2,'single');
    ei.last_bin       = -1;

    ei.ctr_tail = ei.eyeball_ctr;

    % -------- debug
    ei.verbose = 0;
end

