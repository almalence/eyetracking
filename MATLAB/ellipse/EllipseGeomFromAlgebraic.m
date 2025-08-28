% -------------------------------------------------------------------------
%
%    Extract ellipse center, semiaxes, rotation from algebraic parameters
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [cxy, semx, semy, phi] = EllipseGeomFromAlgebraic(A)

    a = A(1); b = A(2)/2; c = A(3); d = A(4)/2; f = A(5)/2; g = A(6);

    cxy = [ (c*d - b*f)/(b^2-a*c)  (a*f - b*d)/(b^2-a*c) ];
    
    semx = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
    semy = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));

    % phi range: [-pi/4 .. 3*pi/4]
    phi = (mod(atan2(2*b, a-c)+pi+pi/2,2*pi)-pi/2)/2;
end
