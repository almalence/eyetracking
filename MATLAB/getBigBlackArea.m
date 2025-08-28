% -------------------------------------------------------------------------
%
%    Find pupil location on eye image
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [xs, ys, psr] = getBigBlackArea(im, imthr, valid)

    [h, w] = size(imthr);

    % operate on 4x sub-sampled data (input is 2x subsampled already)
    h = floor(h/2);
    w = floor(w/2);
    imthr = imthr(1:2:h*2, 1:2:w*2);
    im = im(1:2:h*2, 1:2:w*2);
    valid = single(valid / 2);

    stx = valid(1)-valid(3)+1;
    enx = valid(1)+valid(3);
    sty = valid(2)-valid(4)+1;
    eny = valid(2)+valid(4);

    % for matlab it is faster to pre-compute validArea
    [yc,xc] = ndgrid(0:h-1, 0:w-1);
    xc = abs( xc-valid(1) );
    yc = abs( yc-valid(2) );
    validArea = xc <= valid(3) & yc <= valid(4) & xc+yc <= valid(5);
    
    % horizontal digest of thresholded image: fill pixels with the length of continuous vertical strand they belong to
    dighor = zeros(h,w);
    weight = zeros(h,w);
    for x=stx:enx
        strand = 0;
        for y=sty:eny
            if validArea(y,x)
                if imthr(y,x)
                    if strand==0
                        w1 = im(max(1,y-1),x)-im(y,x);
                    end
                    strand=strand+1;
                elseif strand
                    % apply weight, depending on how high is the minimum gradient at the ends of the strand
                    % assumption is that pupil edge is well-defined,
                    % erroneously thresholded non-pupil areas are unlikely to have such property
                    w2 = im(y,x)-im(y-1,x);
                    ws = clip(floor(min(w1,w2)/4)-2, 0, 4);
                    dighor(y-strand:y-1, x) = strand;
                    weight(y-strand:y-1, x) = ws;
                    strand = 0;
                end
            end
        end
        if strand
            w2 = im(min(h,eny+1),x)-im(eny,x);
            ws = clip(floor(min(w1,w2)/4)-2, 0, 4);
            dighor(eny-strand+1:eny, x) = strand;
            weight(eny-strand+1:eny, x) = ws;
        end
    end

    % vertical scan: find horizontal strand that crosses vertical strands with the largest total number of pixels in these strands
    maxstrand = 0;
    maxarea = 0;
    xs=0; ys=0;
    for y=sty:eny
        strand = 0;
        area   = 0;
        startx = 0;
        for x=stx:enx
            if dighor(y,x)
                if strand == 0, startx = x; end
                strand = strand + dighor(y,x) * weight(y,x);
                area   = area + dighor(y,x);
            else
                if strand > maxstrand
                    xs = floor((x-1+startx)/2);
                    ys = y;
                    maxstrand = strand;
                    maxarea   = area;
                end
                strand = 0;
                area   = 0;
            end
        end
        % Note: this is optional as we do not expect pupil to be at the very edge of the frame
        if strand > maxstrand
            xs = floor((w+startx)/2);
            ys = y;
            maxstrand = strand;
            maxarea   = area;
        end
    end

    % here: xs, ys - roughly the center of pupil
    xs = (xs-1)*4+3;
    ys = (ys-1)*4+3;
    
    % roughly the radius of a pupil
    psr = 4*sqrt(maxarea/pi);

end
