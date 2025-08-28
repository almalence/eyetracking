% -------------------------------------------------------------------------
%
%    Simple bilateral-like filter to reduce gaze direction jitter
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function gv_smooth = smoothGazeVector(gaze_vector, gv_smooth, gv_smooth_diff, fps)

    d = abs(gv_smooth-gaze_vector);
    w = clip((1-d/gv_smooth_diff)*fps/4, 0, fps/4);

    gv_smooth = (gv_smooth.*w + gaze_vector) ./ (w+1);

    gv_z2 = 1-sum(gv_smooth(1:2).^2);
    if gv_z2>0
        gv_smooth(3) = sqrt(gv_z2);
    else
        gv_smooth(3) = gaze_vector(3);
        gv_smooth = gv_smooth ./ norm(gv_smooth);
    end

end
