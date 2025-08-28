% -------------------------------------------------------------------------
%
%    Fit ellipse and calculate distance-weighted fit error
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [A, npx] = ellipseFitAndDistance(x, y, use, maxdist, min_r, max_r)

    % fit ellipse to the pixels at the pupil edge
    A = EllipseDirectFit([x(use) y(use)]);

    [ell_ok, d_ell] = EllipseDistance(A, x, y, min_r, max_r);
    if ~ell_ok, A = []; npx = 0; return; end

    dist = abs(d_ell);

    % if pixel closer than maxdist - pixel is weighted higher
    % if perfectly on an ellipse - the point will count as 2
    npx = sum(1./(0.5+dist/(2*maxdist)));

end
