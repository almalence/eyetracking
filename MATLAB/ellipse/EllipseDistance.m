% -------------------------------------------------------------------------
%
%    Calculate distance between points and ellipse edge
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ell_ok, d_ell] = EllipseDistance(A, x, y, min_r, max_r)

    d_ell = [];
    ell_ok = true;

    % find distance norm to switch to pixel coordinate space
    % by computing how far is the center of the ellipse from its edges
    [cxy, semx, semy, ~] = EllipseGeomFromAlgebraic(A);

    % if A does not correspond to an ellipse or ellipse is too small or too large
    if ~isreal(semx) || ~isreal(semy) || (semx<min_r && semy<min_r) || (semx>max_r) || (semy>max_r)
        ell_ok = false;
        return;
    end

    nrm = (semx+semy)/2 / abs( A(1)*cxy(1)^2 + A(2)*cxy(1)*cxy(2) + A(3)*cxy(2)^2 + A(4)*cxy(1) + A(5)*cxy(2) + A(6) );

    % fitting parameters are: ax^2 + bxy + cy^2 +dx + ey + f
    d_ell = nrm * ( A(1)*x.^2 + A(2)*x.*y + A(3)*y.^2 + A(4)*x + A(5)*y + A(6) );

end
