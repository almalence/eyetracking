% -------------------------------------------------------------------------
%
%    Parse eye images dataset captured with Somnium VR1
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [eye_cfg, D, datazip] = parseSomniumVR1Dataset(dataset, fps, N)

    % ------------ Parse dataset CSV file

    dataparams = [dataset '.csv'];
    opts = detectImportOptions(dataparams);
    opts.SelectedVariableNames = ["imagefile", "eye", "vid", "gaze_x", "gaze_y"];
    opts.DataLines = [2 N+1];
    D = readcell(dataparams, opts);

    skip = round(25/fps);

    % fill in gaze degree gronund truth
    for i=1:length(D)
        vid    = D{i,3};
        D{i,2} = D{i,4};    % gaze X
        D{i,3} = -D{i,5};   % gaze Y
        D{i,5} = vid;       % position number count
        D{i,4} = 0;         % flag indicating the new position start
        D{i,6} = 0; D{i,7} = 0; D{i,8} = 0; % eyeball ground truth position - unknown
    end

    % scene change in past 300ms flag 
    % For a reference, Tobii whitepaper recommend: "80ms before and 80ms after eye gaze direction change"
    % but this is for 'guided' direction change, while we change location instantaneously
    v = [D{:,5}];
    v(2:end) = v(2:end) ~= v(1:end-1);
    v(1) = 1;
    n300 = ceil(25*0.3/skip);
    v = conv(v, ones(n300,1));
    v = v(1:end-(n300-1));
    D(:,4) = num2cell(v);

    % decimate to match fps requested
    D = D(1:skip:end,:);

    datazip = '.zip';

    % ------------ Initialize eye config

    % Horizontal FOV of the eye-capturing camera
    % Negative sign here indicates chief ray converging from camera entrance pupil towards eye
    % This may happen in cases where camera is behind HMD lens and image is formed
    % through some sort of concave mirror surface
    camera_hfov = -30;

    % maximum/minimum iris radius in pixels: closest/farthest distance from the camera
    max_iris_r = 120;
    min_iris_r = 65;

    % how much to rotate camera frame to make eye nearly-horizontal
    camera_zrot = 0;

    % camera-to-HMD space conversion: rotations (around horizontal axis, vertical axis, in frame plane), translations (x,y,z)
    cam2hmd = [0 22 0  25 0 48];
    
    % Array describing valid area to constrain search for pupil in a form of:
    % ctr_x, ctr_y, horz_dist, vert_dist, diag_dist (distances from ctr), corner_brighten*256
    % diag_dist = 0.7*(horz_dist+vert_dist) to cut corners that can not be reached by pupil
    valid_area = int32([280 200  230 200  280  128]);
    
    % -1 for pupil brightness threshold = per frame auto-detect
    eye_cfg = initEyeConfig(540, 400, min_iris_r, max_iris_r, valid_area, -1, camera_hfov, camera_zrot, fps, cam2hmd);

end
