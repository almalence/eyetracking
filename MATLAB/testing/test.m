% -------------------------------------------------------------------------
%
%    Test eye pose estimation
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

clearvars; close all;

addpath([pwd '/dataload']);
addpath([pwd '/..'], [pwd '/../ellipse'], [pwd '/../draw']);

dataset = 'dataset/SomniumVR1/01';
[eye_cfg, D, datazip] = parseSomniumVR1Dataset(dataset, 50, 1440);

eye_calib = struct( 'iris2hrot',12.19, 'iris2vrot',10.15, 'iris2pupil',0.30, 'pupil_decenter',[0.72 -0.21], 'eyeball_ctr',[-102 44], 'pca_angle',[0.00 4.75], 'pupil_refr',1.121, 'pupil_angle',157.50, 'pupil_stretch',0.945, 'gaze_alpha',[-0.53 -2.26] );

eye_cfg = applyCalibration(eye_cfg, eye_calib);

zf = io.ZipFile([dataset datazip]);

M=length(D);
gaze_deg  = zeros(M,2);
no_pupil = 0;

% display intermediary detection algorithm images
eye_cfg.verbose = 2;

for i=1:M
    cam = im2gray(zf.readImage(D{i,1}));           % image from IR-camera
    cam = somniumCamDewarp(cam, 0);

    % dump data for cmdline C version
    % imwrite(cam, sprintf('%05d.bmp', i-1));

    eye_cfg = getEyePose(cam, eye_cfg);

    gaze_deg(i,:) = 180/pi * asin(-eye_cfg.gaze_vector(1:2));

    fprintf ("frame: %03d iris: [%6.2f %6.2f %6.2f]  eyeball: [%6.2f %6.2f %6.2f]  gazevec: [%5.2f %5.2f %5.2f]  gaze angle: [%5.1f %5.1f]\n", ...
        i, eye_cfg.iris3d(1), eye_cfg.iris3d(2), eye_cfg.iris3d(3), ...
        eye_cfg.eyeball3d(1), eye_cfg.eyeball3d(2), eye_cfg.eyeball3d(3), ...
        eye_cfg.gaze_vector(1),eye_cfg.gaze_vector(2),eye_cfg.gaze_vector(3), gaze_deg(i,1), gaze_deg(i,2) );

    if eye_cfg.blink==1, no_pupil = no_pupil + 1; end

    if eye_cfg.verbose
        figure(1); title(i); pause;
        % do the verbose run for the first 20 frames only
        if i>=20, eye_cfg.verbose = 0; end
    end
end
fprintf("\n");

% -------- Statistics

gaze_deg_true = [[D{:,2}]' [D{:,3}]'];

% measure performance according to Tobii white paper guidelines
vid = [D{:,5}];
% skip first 160ms after each scene change in statistic calculation (as stipulated by Tobii guidelines)
v=~[D{:,4}];
[prec, accuracy, ranges] = precisionAndAccuracyPerAngleRange(gaze_deg(v,:), gaze_deg_true(v,:), vid(v));

% gaze direction measurements
fprintf("Gaze overall precision: %5.2f degree, accuracy: %5.2f degree\n", mean(prec), mean(accuracy));
rangestr={" 0-10", "10-20", "20-25", "25-30"};
for i=1:4
    fprintf("Gaze range %s precision: %5.2f degree, accuracy: %5.2f degree\n", rangestr{i}, mean(prec(ranges==i)), mean(accuracy(ranges==i)));
end

% gaze direction graphs
figure;
vx = 1:M; vx = vx(v);
subplot(2,1,1); plot(vx, gaze_deg(v,1), vx, gaze_deg_true(v,1)); title('X degree'); grid on; grid minor; ylim([-30 30]);
subplot(2,1,2); plot(vx, gaze_deg(v,2), vx, gaze_deg_true(v,2)); title('Y degree'); grid on; grid minor; ylim([-30 30]); axis ij;
