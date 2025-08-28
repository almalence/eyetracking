% -------------------------------------------------------------------------
%
%    Calibration of pupil non-circularity (elongation) and pupil circular axis (PCA) angles
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ei, ctrn, status] = calibratePupilElongationPCA(ei, cam, gv_cam, T2anc, Tuse, ctr, pxl2mm, N, verbose)

M = length(cam);    % total frames
status = false;
ctrn = [];

ei.smooth_2deye  = -1;      % use eyeball center found in the first pass
ei.eyeball2d = ctr; ei.ctr_tail = ctr; ei.eyeball_ctr = ctr;

pxy = cell(N,1);
cxy = cell(N,1);
el_ratio0 = cell(N,1);
iris_ratio = cell(N,1);
numpos = zeros(N,1);
for i=1:M
    ianc = T2anc(i);
    ei = getEyePose(cam{i}, ei);
    if Tuse(i) && ei.blink < 0.5
        pxy{ianc} = [pxy{ianc} {ei.pxy}];
        cxy{ianc} = [cxy{ianc} {ei.cxy}];
        el_ratio0{ianc} = [el_ratio0{ianc} {ei.el_ratio0}];
        iris_ratio{ianc} = [iris_ratio{ianc} {ei.iris_ratio}];
        
        numpos(ianc) = numpos(ianc) + 1;
    end
end

if any(numpos([1:4 N])<3)
    warning('Pupil calibration failed, not enough good frames...');
    return;
end

% compute el_ratio/phi_minor corrections (pca, elongation) from reference given by gaze vectors in camera space vs measured

% there might be slight differences between eyeball centers found in this pass and the next one
% since iris centers are used here for camRayProjToMinor() and compensatePCA()
% and pupil centers are used in the next pass (and in normal tracking)

% find the best pup_ang and pup_stretch via iterative minimization
best_pup_refr = ei.pupil_refr;
best_pup_ang = 0;
best_pup_stretch = 1;
best_ax = 0;
best_ay = 0;

% first run - find pupilaary circular axis (pca)
[h,w] = size(cam{1});
img_ctr = [w/2 h/2]+1; % +1 to correct for Matlab 1-based indexes
ctr_unrot = rotcam(ctr, -ei.camera_rot) + img_ctr;

% only calibrate for pca if positions above 9 degree from both sides relative to the eye center are present
pcax_range = 0;
pcay_range = 0;
if min(gv_cam(:,1)) < -0.15 && max(gv_cam(:,1)) > 0.15
    pcax_range = -15:0.25:15;
end
if min(gv_cam(:,2)) < -0.15 && max(gv_cam(:,2)) > 0.15
    pcay_range = -15:0.25:15;
end

% assuming angle order: L R T B ... C

best_diff = 1e10;
for ax=pcax_range
    diff = diffFromPCA(pxy([1 2 N]), cxy([1 2 N]), gv_cam([1 2 N],:), ctr_unrot, img_ctr, ei.iris_z, ax, 0, 0);
    if diff<best_diff
        best_ax = ax;
        best_diff = diff;
    end
end

best_diff = 1e10;
for ay=pcay_range
    diff = diffFromPCA(pxy([3 4 N]), cxy([3 4 N]), gv_cam([3 4 N],:), ctr_unrot, img_ctr, ei.iris_z, best_ax, ay, 0);
    if diff<best_diff
        best_ay = ay;
        best_diff = diff;
    end
end

% re-run horizontal pass with best vertical pca estimation
best_diff = 1e10;
for ax=pcax_range
    diff = diffFromPCA(pxy([1 2 N]), cxy([1 2 N]), gv_cam([1 2 N],:), ctr_unrot, img_ctr, ei.iris_z, ax, best_ay, 0);
    if diff<best_diff
        best_ax = ax;
        best_diff = diff;
    end
end

ei.pca_angle = [best_ax best_ay];

% second run - find optimal parameters to best match direction towards center (but not scale) via pupil elongation:
% - ignore camera perspective projection  ... cra-proj has heavy influence on angles ??? ... perform same correction as to gv_cam in first run ... ???
% - ignore pupil refraction
% - parameters affecting: pup_ang, pup_stretch

% only calibrate for pupil elongation if sufficienctly front-looking position is present
if max(gv_cam(:,3)) < 0.98
    pup_ang_range1 = 0;
    pup_ang_range2 = 0;
    pup_stretch_range1 = 1;
    pup_stretch_range2 = 0;
else
    pup_ang_range1 = 0:pi/8:pi;
    pup_ang_range2 = -pi/8:pi/32:pi/8;
    pup_stretch_range1 = 0.8:0.05:1;
    pup_stretch_range2 = -0.05:0.005:0.05;
end

best_diff = 1e10;
for pup_ang = pup_ang_range1
    for pup_stretch = pup_stretch_range1
        v = [pup_ang; pup_stretch; best_ax; best_ay; best_pup_refr];
        diff = ctrDiffFromPupilRatio(pxy(1:N-1), cxy(1:N-1), el_ratio0(1:N-1), iris_ratio(1:N-1), ctr_unrot, ...
            img_ctr, pxl2mm(1:N-1), ei.max_el_ratio, ei.iris_z, ei.iris2eye_mm, ei.iris2hrot_mm, ei.iris2vrot_mm, v, 0);
        if diff<best_diff
            best1_pup_ang = pup_ang;
            best1_pup_stretch = pup_stretch;
            best_diff = diff;
        end
    end
end

best_diff = 1e10;
for pup_ang = best1_pup_ang+pup_ang_range2
    for pup_stretch = best1_pup_stretch+pup_stretch_range2
        v = [pup_ang; pup_stretch; best_ax; best_ay; best_pup_refr];
        diff = ctrDiffFromPupilRatio(pxy(1:N-1), cxy(1:N-1), el_ratio0(1:N-1), iris_ratio(1:N-1), ctr_unrot, ...
            img_ctr, pxl2mm(1:N-1), ei.max_el_ratio, ei.iris_z, ei.iris2eye_mm, ei.iris2hrot_mm, ei.iris2vrot_mm, v, 0);
        if diff<best_diff
            best_pup_ang = pup_ang;
            best_pup_stretch = pup_stretch;
            best_diff = diff;
        end
    end
end

if verbose, fprintf("Preliminary pupil stretching: %6.4f\n", best_pup_stretch); end

[~, ctrn] = ctrDiffFromPupilRatio(pxy(1:N-1), cxy(1:N-1), el_ratio0(1:N-1), iris_ratio(1:N-1), ctr_unrot, ...
        img_ctr, pxl2mm(1:N-1), ei.max_el_ratio, ei.iris_z, ei.iris2eye_mm, ei.iris2hrot_mm, ei.iris2vrot_mm, ...
        [best_pup_ang; best_pup_stretch; best_ax; best_ay; best_pup_refr] , 0);
ctrn = rotcam(ctrn-img_ctr, ei.camera_rot);

ei.pupil_refr = best_pup_refr;
ei.pupil_angle = best_pup_ang*180/pi;
ei.pupil_stretch = best_pup_stretch;

status = true;
