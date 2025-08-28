% -------------------------------------------------------------------------
%
%    Fill-in eye config fields from eye calibration data
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

%  ec.iris2hrot
%  ec.iris2vrot
%  ec.iris2pupil        how much pupil is closer to the eyeball center than the iris plane
%  ec.pupil_decenter
%  ec.eyeball_ctr
%  ec.pca_angle
%  ec.pupil_refr
%  ec.pupil_angle
%  ec.pupil_stretch
%  ec.gaze_alpha
%
% Default values:
% struct('iris2hrot', 11, 'iris2vrot', 9, 'iris2pupil', 0.3, 'pupil_decenter', [0 0], 'eyeball_ctr', [0 0], 'pca_angle', [0 0], 'pupil_refr', 1.121, 'pupil_angle', 0, 'pupil_stretch', 1, 'gaze_alpha', [4 -2])
%
function ei = applyCalibration(ei, ec)

    ei.iris2hrot_mm      = ec.iris2hrot;
    ei.iris2vrot_mm      = ec.iris2vrot;
    ei.iris2eye_mm       = ei.iris2vrot_mm;
    ei.pupil2eye_mm      = ei.iris2vrot_mm - ec.iris2pupil;

    ei.pupil_decenter_mm = ec.pupil_decenter;

    ei.eyeball_ctr       = ec.eyeball_ctr;
    ei.ctr_tail          = ei.eyeball_ctr;
    ei.eyeball2d         = ei.eyeball_ctr;

    ei.pca_angle         = ec.pca_angle;
    ei.pupil_refr        = ec.pupil_refr;
    ei.pupil_angle       = ec.pupil_angle;
    ei.pupil_stretch     = ec.pupil_stretch;

    ei.alpha_mat         = rotationMat([ec.gaze_alpha(2) -ec.gaze_alpha(1) 0], 0);
end
