% -------------------------------------------------------------------------
%
%    Plot centers of eye features
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function drawEyeCenter(exy, ixy, clr, img_ctr, ang, verbose)

    if verbose==0, return; end

    exy = rotcam(exy, -ang) + img_ctr;
    if ~isempty(ixy)
        ixy = rotcam(ixy, -ang) + img_ctr;
    end

    figure(1);
    
    subplot(1,2,1);
    hold on;
    plot (exy(1), exy(2), [clr 'o']);
    if ~isempty(ixy)
        plot ([exy(1) ixy(1)], [exy(2) ixy(2)], clr);
    end
    hold off;

    subplot(1,2,2);
    hold on;
    plot (exy(1), exy(2), [clr '*']);
    hold off;

end
