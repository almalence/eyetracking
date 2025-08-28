% -------------------------------------------------------------------------
%
%    Returns how much the coordinate vector should be scaled
%    depending on the angle, to change unit circle into given ellipse
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function scl = ellipseScales(iris_dirs, phi, el_r)
    x=cos(iris_dirs-phi);
    y=sin(iris_dirs-phi);

    dirs_ell = atan2(y, x*el_r);
    % need limit to 1 here due to precision issues
    scl = min(1, sqrt( cos(dirs_ell).^2 + (el_r*sin(dirs_ell)).^2) );
end
