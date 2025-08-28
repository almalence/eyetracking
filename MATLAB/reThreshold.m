% -------------------------------------------------------------------------
%
%    Pupil image re-thresholding using local neighborhood of the pupil
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function imthr = reThreshold(im, xs, ys, psr, pupthr, max_pup_r)

    [h, w] = size(im);

    % use additional area around pupil
    psr2 = max(32, ceil(psr)*4);
    stx = max(2,xs-psr2); enx = min(w-1,xs+psr2)-1;
    sty = max(2,ys-psr2); eny = min(h-1,ys+psr2)-1;
    % for speed reasons only use even rows/columns to compute histogram, if pupil is large
    if psr2>=max_pup_r/2
        stx = stx-mod(stx-1,2);
        sty = sty-mod(sty-1,2);
        x = stx:2:enx; y = sty:2:eny;
    else
        x = stx:enx; y = sty:eny;
    end

    tile = im(y,x);

    % pupil area under initial thresholding
    parea = sum(tile <= pupthr, 'all');

    % we're using double the initially-detected pupil radius,
    % therefore there will always be 'transition' pixels (at the pupil edge)
    % and somewhat more non-pupil pixels in a tile
    hst = histcounts(tile, 0:256);

    % smooth-out histogram
    hstf = filter(ones(1,8)/8, 1, [hst([1 1 1]) hst hst([255 255 255 255])]);
    hstf = round(hstf(7+(1:256))); % round to match C version

    sumhstf = cumsum(hstf);
    % once enough pixels under histogram cumulative sum to cover twice the pupil area
    % - we are definitely above the edge brightness
    edge = find(sumhstf > parea + sumhstf(pupthr), 1);
    if ~isempty(edge)
        pupthr2 = idivide(3*pupthr + edge, 4);
    else
        pupthr2 = pupthr;
    end

    imthr = im <= pupthr2;
end
