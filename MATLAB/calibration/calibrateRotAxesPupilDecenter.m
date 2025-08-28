% -------------------------------------------------------------------------
%
%    Calibration of eye rotational axes positions and pupil decentration
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ei, ctr, iris2pupil, pxl2mm, slippage, fig, status] = calibrateRotAxesPupilDecenter(ei, cam, gv_cam, T2anc, Tuse, primary_anc, N, verbose)

M = length(cam);    % total frames
status = false;
fig = [];
ctr = [0 0];
iris2pupil = 0.3;
slippage = zeros(N,3); slippage(:,3) = 1;

ei.smooth_2deye  = 1; % eyeball center detection jitter may cause center 'flipping' to erroneous cluster without smoothing
ei.smooth_2dtail = 0;
ei.smooth_3deye  = 0;
ei.smooth_iris   = 3; % erroneous iris radius detecitons with low confidence will get in the way without iris smoothing

% Note: in case of poor iris detection - pupil decenter calibraiton will fail but pupil decentration can be ignored
% since even typically high decenter of 0.9mm will only result in
% 0.9*sin(15) = 0.22mm systematic error in rotation axis distance in the very worst case

pxl2mm   = zeros(N,1);
pupil2d  = zeros(N,2);
iris2dx  = zeros(N,1);
iris2dy  = zeros(N,1);
nirisx   = zeros(N,1);
nirisy   = zeros(N,1);
npupil   = zeros(N,1);

iris_r1  = cell(N,1);
iris_r2  = cell(N,1);

ei.el_known = true;
ei.iris_recenter = 2; % re-estimate iris both horizontally and vertically (for pupil decenter estimation)
for i=1:M
    ianc = T2anc(i);

    % pre-fill iris ratio of previous frame with ground-truth value
    ei.iris_ratio = gv_cam(ianc,3);
    % pre-fill ground truth ratio and phi (will overwrite values detected from pupil)
    ei.el_ratio =  gv_cam(ianc,3);
    ei.phi_minor = atan2(gv_cam(ianc,2), gv_cam(ianc,1))+pi/2;
    ei = getEyePose(cam{i}, ei);

    if Tuse(i) && ei.blink < 0.5
        pxl2mm(ianc)   = pxl2mm(ianc) + ei.pxl2mm;
        pupil2d(ianc,:) = pupil2d(ianc,:) + ei.iris2d_p;
        npupil(ianc)   = npupil(ianc) + 1;

        % for decenter - use locations where iris2d was estimated from 2 radii
        if ei.iris_r1 && ei.iris_r2
            iris2dx(ianc)   = iris2dx(ianc) + ei.iris2d(1);
            nirisx(ianc)   = nirisx(ianc) + 1;

            iris2dy(ianc)   = iris2dy(ianc) + ei.iris2d(2);
            nirisy(ianc)   = nirisy(ianc) + 1;
        end

        % record separate estimations of iris radius from left and right edge
        % to correct later on with pupil decenter and adjust resulting pxl2mm values
        if ei.iris_r1>0, iris_r1{ianc} = [iris_r1{ianc} ei.iris_r1]; end
        if ei.iris_r2>0, iris_r2{ianc} = [iris_r2{ianc} ei.iris_r2]; end
    end
end
% set back default values
ei.el_known = false;
ei.iris_recenter = 1;

if any(npupil<3)
    warning('Rotation axes/slippage calibration fail, pupil heavily occluded...');
    return;
end

% averaged per-position pxl2mm, iris2d
pxl2mm   = pxl2mm ./ npupil;
pupil2d = pupil2d ./ npupil;
iris2d   = [iris2dx./nirisx iris2dy./nirisy];
iris2d(isnan(iris2d))=pupil2d(isnan(iris2d));

% averaged separate estimations of iris radius from left and right edge
% excluding outliers
ir_r1 = zeros(N,1);
ir_r2 = zeros(N,1);
n_r1 = zeros(N,1);
n_r2 = zeros(N,1);
for i=1:N
    if ~isempty(iris_r1{i})
        ir_r1(i) = mean(iris_r1{i});
        valid = abs(iris_r1{i}-ir_r1(i)) < 0.2*ir_r1(i);
        ir_r1(i) = mean(iris_r1{i}(valid));
        n_r1(i) = sum(valid);
    end

    if ~isempty(iris_r2{i})
        ir_r2(i) = mean(iris_r2{i});
        valid = abs(iris_r2{i}-ir_r2(i)) < 0.2*ir_r2(i);
        ir_r2(i) = mean(iris_r2{i}(valid));
        n_r2(i) = sum(valid);
    end
end

% estimate pxl2mm for all locations from location closest to the camera
% using single base value of pxl2mm for better consistency
% pxl2mm(:) = pxl2mm(primary_anc); % flat model

z_pri = ei.iris_z + ei.iris2eye_mm * (1-gv_cam(primary_anc,3)) / pxl2mm(primary_anc);
for i=1:N
    if i ~= primary_anc
        z = ei.iris_z + ei.iris2eye_mm * (1-gv_cam(i,3)) / pxl2mm(primary_anc);
        pxl2mm(i) = pxl2mm(primary_anc) * z/z_pri;
    end
end

% compute pupil decenter
pup_dec_h = 0; pup_dec_v = 0;
nh = 0; nv = 0;
for i=1:N
    if nirisx(i)>2
        sh = sqrt(1-gv_cam(i,1)^2);
        if sh>0.5
            pup_dec_h = pup_dec_h + (pupil2d(i,1)-iris2d(i,1)) / sh * pxl2mm(i);
            nh = nh+1;
        end
    end

    if nirisy(i)>2
        sv = sqrt(1-gv_cam(i,2)^2);
        if sv>0.5
            pup_dec_v = pup_dec_v + (pupil2d(i,2)-iris2d(i,2)) / sv * pxl2mm(i);
            nv = nv+1;
        end
    end
end

% pupil decenter is of low importance for off-axis camera, so it can be skipped
if nh
    pup_dec_h = pup_dec_h / nh;
else
    if verbose, disp('Horizontal pupil decenter calibration skipped...'); end
end
if nv
    pup_dec_v = pup_dec_v / nv;
else
    if verbose, disp('Vertical pupil decenter calibration skipped...'); end
end

% re-estimate pxl2mm accounting for pupil decenter
ir_r1 = ir_r1 - pup_dec_h./pxl2mm;
ir_r2 = ir_r2 + pup_dec_h./pxl2mm;
iris_r = (ir_r1.*n_r1 + ir_r2.*n_r2) ./ (n_r1+n_r2);
pxl2mm_obs = pxl2mm;
valid = n_r1+n_r2 >= 3;
pxl2mm_obs(valid) = ei.hmn_iris_mm/2./iris_r(valid);    % observed pxl2mm, used for detection of 3d slippage during calibration
pxl2mm = pxl2mm .* pxl2mm_obs(primary_anc)/pxl2mm(primary_anc);

% compute distance to axis (as gaze vector length) from the length between
% reference vector projections and distances between iris2d's
% assuming angle order: L R T B ...
d_ref_h = abs(gv_cam(1,1) - gv_cam(2,1));
dist_mm_p_h = abs( pupil2d(1,1)*pxl2mm(1) -  pupil2d(2,1)*pxl2mm(2)) / d_ref_h;
dist_mm_h = dist_mm_p_h;
d_ref_v = abs(gv_cam(3,2) - gv_cam(4,2));
dist_mm_p_v = abs( pupil2d(3,2)*pxl2mm(3) -  pupil2d(4,2)*pxl2mm(4)) / d_ref_v;
dist_mm_v = dist_mm_p_v;

if min(nirisx(1), nirisx(2)) >= 2
    dist_mm_h = abs( iris2d(1,1)*pxl2mm(1) -  iris2d(2,1)*pxl2mm(2)) / d_ref_h;
else
    if verbose, fprintf('Iris-to-pupil calibration skipped ...\n'); end
end

if min(nirisy(3), nirisy(4)) >= 2
    dist_mm_v = abs( iris2d(3,2)*pxl2mm(3) -  iris2d(4,2)*pxl2mm(4)) / d_ref_v;
end

iris2pupil_h = dist_mm_h-dist_mm_p_h+0.3;
iris2pupil_v = dist_mm_v-dist_mm_p_v+0.3;
% horizontal iris ellipse estimation is more reliable
%iris2pupil = ( iris2pupil_h + iris2pupil_v ) /2 ;
% not expecting pupil-iris distance to be above 1mm
iris2pupil = clip(iris2pupil_h, -1, 1);

% very small and very large vlaues are possible with prescription glasses
if min(dist_mm_h,dist_mm_v)<4 || max(dist_mm_h,dist_mm_v)>25
    warning('Rotation axes calibration fail, calibration out of range: H: %5.2fmm  V: %5.2fmm ...', dist_mm_h, dist_mm_v);
    return;
else
    ei.iris2hrot_mm = dist_mm_h;
    ei.iris2vrot_mm = dist_mm_v;
    ei.iris2eye_mm = ei.iris2vrot_mm;
    ei.pupil2eye_mm = ei.iris2vrot_mm - iris2pupil;
    if verbose, fprintf('Rotation axes calibration (from iris position) complete: H axis: %5.2fmm  V axis: %5.2fmm\n', ei.iris2hrot_mm, ei.iris2vrot_mm); end

    ei.pupil_decenter_mm = [pup_dec_h pup_dec_v];
    if verbose
        fprintf('Pupil decenter calibration complete: [%5.2f %5.2f]\n', ei.pupil_decenter_mm(1), ei.pupil_decenter_mm(2));
        fprintf('Iris2pupil: %5.2f  (H: %5.2fmm V: %5.2fmm)\n', iris2pupil, iris2pupil_h, iris2pupil_v);
    end
end

% estimate eyeball center on image plane and slippage of each location
clr = {"r","g","b","c","m","y","k","#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"};
[ctr, slippage, fig] = getCtrAndSlip(iris2d, -gv_cam(:,1:2), N, verbose, clr);

slippage(:,3) = pxl2mm_obs./pxl2mm;

status = true;
