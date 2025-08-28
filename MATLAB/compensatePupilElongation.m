% -------------------------------------------------------------------------
%
%    Compensate pupil eged points coordinates for pupil elongation
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function xy = compensatePupilElongation(xy, el_r, pup_ang, pup_stretch)

    % reduce stretching when looking towards pupil at an angle
    % no stretching for angles above ~40 degree
    stretch_adj = (pup_stretch-1) * max(0, 1-(1-el_r)*3);

    % distance from the stretch axis
    d = -xy(:,1) * sin(pup_ang) + xy(:,2) * cos(pup_ang);

    % apply stretch
    xy = xy - d * stretch_adj .* [-sin(pup_ang) cos(pup_ang)];
    
end
