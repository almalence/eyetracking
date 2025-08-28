% -------------------------------------------------------------------------
%
%    Smooth eyeball center location in 3d temporally
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [eyeball3d, eyeball3d0x] = FilterEyeballPosition(eyeball3d, old3d0x, iris3d, iris2eye_mm, iris2hrot_mm, strength, old_confid, eyeball_confid)

    % axis for vertical eye rotation is closer to cornea than horizontal rotation axis (which is coinciding with eyeball center)
    % adjust coordinates of stabilized point for that
    gaze_norm = (eyeball3d-iris3d) / iris2eye_mm;
    dx = -(iris2hrot_mm-iris2eye_mm) * gaze_norm(1);

    eyeball3d0x = eyeball3d - [dx 0 0];

    if old_confid   % no filtering for the first reading

        upd_weight = min(1, eyeball_confid / old_confid);
        % x,y: full correction > 30 degree            
        wx = clip((1-gaze_norm(3))*8, 0.1, 1);
        wy = wx;
        % z: full correction < 30 degree, almost no correction > 45
        wz = clip((gaze_norm(3)-0.5)*5, 0.1, 1);
        upd_weight = upd_weight * [wx wy wz];

        eyeball3d0x = (old3d0x.*(2^strength-1) + eyeball3d0x.*upd_weight) ./ ((2^strength-1) + upd_weight);

        eyeball3d = eyeball3d0x + [dx 0 0];
    end
end
