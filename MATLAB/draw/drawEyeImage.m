% -------------------------------------------------------------------------
%
%    Draw eye image
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function drawEyeImage(im, verbose)

    if verbose==0, return; end

    figure(1); subplot(1,2,2);
    imshow(im);
end
