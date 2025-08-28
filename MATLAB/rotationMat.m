% -------------------------------------------------------------------------
%
%    Rotation matrix from rotation angles in degree
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% X-Y-Z Euler rotation
function m = rotationMat(r, inv)
    rx = r(1)*pi/180;
    ry = r(2)*pi/180;
    rz = r(3)*pi/180;

    if inv
        rx = -rx;
        ry = -ry;
        rz = -rz;
    end

    RX = [     1       0        0
               0   cos(rx)  -sin(rx)
               0   sin(rx)   cos(rx) ];
    RY = [ cos(ry)     0     sin(ry)
               0       1        0
          -sin(ry)     0     cos(ry) ];
    RZ = [ cos(rz)  -sin(rz)    0
           sin(rz)   cos(rz)    0
               0       0        1   ];
    if inv
        m = RX * RY * RZ;
    else
        m = RZ * RY * RX;
    end
end
