% -------------------------------------------------------------------------
%
%    Load frames from dataset and calibrate eye parameters
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function eye_calib = calibrateSet(dataset, Ncalib, Nframes, parseFn, dewarpFn, verbose)

% anchor locations for desired gaze angles, should follow the order: L R T B .. [additional locations] .. C
if Ncalib==9
    anc = [-15 0; 15 0; 0 -15; 0 15; -15 -15; 15 -15; -15 15; 15 15; 0 0];
else
    anc = [-15 0; 15 0; 0 -15; 0 15; 0 0];
end

% using 25/30 fps for calibration
[ei, D, datazip] = eval([parseFn '(dataset, 30, Nframes)']);

T=1:length(D);
v=[D{T,4}]; % flags for initial 160ms of the scene
gd_true = [[D{T,2}]' [D{T,3}]'];

% select N locations closest to the desired anchor points
N = size(anc,1); % total calibration points

Tst = zeros(1,N);                   % start and
Ten = zeros(1,N);                   % end position of each recording
Tanc = cell(1,N);
T2anc = cell(1,N);
for i=1:N
    d = vecnorm((gd_true-anc(i,:))');
    [~,idx]=min(d);

    Tst(i) = idx;                           % calibration position start
    idx = idx-1 + find(~v(idx:end), 1);     % skip frames flagged as transition
    idx = idx-1 + find(v(idx:end), 1);      % skip till end of scene
    if isempty(idx), idx = length(v)+1; end % if last scene
    Ten(i) = idx-1;
    len = Ten(i)-Tst(i)+1;

    Tanc{i} = Tst(i):Ten(i);
    T2anc{i} = i*ones(1,len);
end

% check if sufficient locations for calibration
for i=1:N
    for j=1:N
        if i~=j && all(anc(i,:) == anc(j,:))
            disp('Not enough calibration points...');
            eye_calib=[];
            return;
        end
    end
end

T = [Tanc{:}];
T2anc = [T2anc{:}];

M=length(T);
zf = io.ZipFile([dataset datazip]);

% read frames into memory
cam=cell(1,M);
for i=1:M
    cam{i} = im2gray(zf.readImage(D{T(i),1}));           % image from IR-camera
    if isa(cam{i}, 'uint16'), cam{i} = uint8(cam{i}/256); end
    if ~isempty(dewarpFn), cam{i} = eval([dewarpFn '(cam{i}, 0)']); end

    % dump data for cmdline C version
    % imwrite(cam{i}, sprintf('%05d.bmp', i-1));
end

disp('Calibration start ...');
[~, eye_calib] = calibrateEye(ei, cam, anc, T2anc, verbose);
