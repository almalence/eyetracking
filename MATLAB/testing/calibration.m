% -------------------------------------------------------------------------
%
%    Eye parameters calibration script
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

clearvars; close all;

addpath([pwd '/dataload']);
addpath([pwd '/..'], [pwd '/../calibration'], [pwd '/../ellipse'], [pwd '/../draw']);

dataset = 'dataset/SomniumVR1_calib/01';
parseFn = 'parseSomniumVR1Dataset'; dewarpFn = 'somniumCamDewarp'; Nframes = 199;

eye_calib = calibrateSet(dataset, 5, Nframes, parseFn, dewarpFn, true);
