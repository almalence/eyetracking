% -------------------------------------------------------------------------
%
%    Calibration of optical vs visal axis angle alpha
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ei, alpha, gaze_deg, ir_ratio, el_ratio, el_pos, ctr_pos] = calibrateAlphaAngle(ei, cam, anc, T2anc, Tuse, ctr, pxl2mm, N)

M = length(cam);    % total frames


% Note: alpha angle is due to the angle between optical and visual eye axis
%       alpha stretching (scale) - either due to the pupil swim effects of HMD or eyewear

ei.valid(:) = 0; % using this to start center estimations afresh (removes trash from the first pass)
ei.smooth_2deye  = 1;
ei.smooth_2dtail = 4;
ei.smooth_3deye  = 0;
ei.smooth_iris   = 3;
ei.eyeball2d = ctr; ei.ctr_tail = ctr; ei.eyeball_ctr = ctr;

rotrad      = ei.camera_rot * pi/180;

gaze_deg  = zeros(M,2);
ir_ratio  = zeros(N,1);
el_ratio  = zeros(N,1);
el_pos    = zeros(N,2);
ctr_pos   = zeros(N,2);
numratpos = zeros(N,1);
use = zeros(size(Tuse), 'logical');

for i=1:M
    ianc = T2anc(i);
    ei = getEyePose(cam{i}, ei);
    gaze_deg(i,:) = 180/pi * asin(-ei.gaze_vector(1:2));

    if Tuse(i) && ei.blink < 0.5
        use(i) = true;

        dv = [sin(ei.phi_minor+rotrad) -cos(ei.phi_minor+rotrad)];
        proj_ang = acos(ei.el_ratio);
        len2dy = ei.iris2eye_mm/pxl2mm(ianc)*sin(proj_ang);
        len2dx = ei.iris2hrot_mm/ei.iris2vrot_mm * len2dy;
        e1h = ei.iris2d + len2dx.*dv;
        e1v = ei.iris2d + len2dy.*dv;

        ctr_pos(ianc,:) = ctr_pos(ianc,:) + [e1h(1) e1v(2)];
        el_pos(ianc,:) = el_pos(ianc,:) + ei.iris2d;
        el_ratio(ianc) = el_ratio(ianc) + ei.el_ratio;
        ir_ratio(ianc) = ir_ratio(ianc) + ei.iris_ratio;
        numratpos(ianc) = numratpos(ianc) + 1;
    end
end

el_ratio = el_ratio./numratpos;
el_pos   = el_pos./numratpos;
ctr_pos  = ctr_pos./numratpos;
ir_ratio = ir_ratio./numratpos;

% Alpha angle
%alpha = mean(anc(primary_anc,:)-gaze_deg(Tuse_primary,:));
alpha = mean(anc(T2anc(use),:)-gaze_deg(use,:));
ei.alpha_mat = rotationMat([alpha(2) -alpha(1) 0], 0);
