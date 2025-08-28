% -------------------------------------------------------------------------
%
%    Histogram-based adaptive pupil thresholding of eye image
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [imthr, pupthr] = thresholdForPupil(im, min_pup_r, pupthr, valid)

    if pupthr<0
        stx = valid(1)-valid(3)+1;
        enx = valid(1)+valid(3);
        sty = valid(2)-valid(4)+1;
        eny = valid(2)+valid(4);
        imvalid = im(sty:eny, stx:enx);

        % exclude corners
        [y,x] = ndgrid(-valid(4)+1:valid(4), -valid(3)+1:valid(3));
        imvalid = max( imvalid, uint8(255 * (abs(x)+abs(y) > valid(5))) );

        % find pupil pixel brightness threshold (histogram image, threshold at some cumulative percentage of pixels area)
        ch = cumsum(histcounts(imvalid, 0:256));
        min_pup_area = pi*min_pup_r*min_pup_r/2;
        pupthr = find(ch>min_pup_area);
        % *1.17 and +6 levels here to account for noise / nonuniformity
        % pupthr = 1.17*pupthr(1) + 6;
        % matching to C code, 75/64=1.17
        pupthr = idivide(75*int32(pupthr(1)), 64) + 6;
    end

    imthr = im <= pupthr;
end
