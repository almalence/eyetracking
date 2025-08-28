% -------------------------------------------------------------------------
%
%    2D coordinte rotation of camera image plane
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function out = rotcam(crd, ang)
    ang = ang * pi/180;

    out = zeros(size(crd), 'single');
    out(:,1) =  crd(:,1)*cos(ang) - crd(:,2)*sin(ang);
    out(:,2) =  crd(:,1)*sin(ang) + crd(:,2)*cos(ang);
end
