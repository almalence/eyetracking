% -------------------------------------------------------------------------
%
%    Pupil decentration correction
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function ixy = correctPupilDecenter(cxy, el_ratio, phi_minor, pxl2mm, decenter_mm, camera_rot)

    decenter = decenter_mm / pxl2mm;
    % rotate decenter to camera space
    decenter = rotcam(decenter, -camera_rot);
    % reduce decentering according to ellipse scale in horizontal/vertical direction
    el_xy(1) = ellipseScales(0-camera_rot*pi/180, phi_minor, el_ratio);
    el_xy(2) = ellipseScales(pi/2-camera_rot*pi/180, phi_minor, el_ratio);
    ixy = cxy - [decenter(1)*el_xy(1) decenter(2)*el_xy(2)];

end
