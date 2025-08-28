% -------------------------------------------------------------------------
%
%    Extract pupil edge points
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [xs, ys, xp, yp, parea] = extractEdgePixels(imthr, xs, ys, psr, max_pup_r)

    [h, w] = size(imthr);
    imthr = int16(imthr);
    max_pup_r = floor(max_pup_r);
    xp = []; yp = [];
    parea = 0;

    % assuming that pieces of pupil are not torn apart by more than that, even value for optimized code
    d = single(max(8, 2*idivide(5*int32(floor(psr)),64)));
    psr = max(16, floor(psr));

    stxmax = max(d+1, xs-max_pup_r);
    enxmax = min(w-d-2, xs+max_pup_r); enxmax = enxmax - mod(enxmax-stxmax, d);
    stymax = max(d+1, ys-max_pup_r);
    enymax = min(h-d-2, ys+max_pup_r); enymax = enymax - mod(enymax-stymax, d);

    % a simplistic version of flood-fill
    % not using a proper visit-all list-based approach since we have
    % 15x15 sites at most, and to cover all of them less than 10 full scans are required

    % at least 4 pixels out of d x d should be flagged to assume square occupied
    minpx = 4;

    stxpup = max(2, xs-psr); enxpup = min(w-1, xs+psr);
    stypup = max(2, ys-psr); enypup = min(h-1, ys+psr);

    % map of occupied squares (2nd bit)
    % and mark initial pupil area (3rd bit)
    % mark only one pixel in a site for speed
    stx = enxmax; enx = stxmax;
    sty = enymax; eny = stymax;
    for y=stymax:d:enymax
        for x=stxmax:d:enxmax
            % skipping every second pixel to match optimized C version
            if sum(imthr(y:2:y+d-1, x:2:x+d-1), 'all') >= minpx/4
                imthr(y,x) = imthr(y,x) + 2;
            end
            if x>=stxpup && x<=enxpup && y>=stypup && y<=enypup
                imthr(y,x) = imthr(y,x) + 4;

                % initial bounding box
                if stx>x, stx=x; end
                if enx<x, enx=x; end
                if sty>y, sty=y; end
                if eny<y, eny=y; end
            end
        end
    end

    % expand initially marked area in connecting sites
    while true
        nnew = 0;
        for y=max(stymax,sty-d):d:min(enymax,eny+d)
            for x=max(stxmax,stx-d):d:min(enxmax,enx+d)
                if bitand(imthr(y,x),2) && bitand(imthr(y,x),4)==0  % candidate to connect
                    if bitand(imthr(max(stymax,y-d),x),4) || ...
                       bitand(imthr(y,max(stxmax,x-d)),4) || bitand(imthr(y,min(enxmax,x+d)),4) || ...
                       bitand(imthr(min(enymax,y+d),x),4)
                        % if connecting at the very edge - pupil detection is incorrect
                        if x==stxmax || x==enxmax || y==stymax || y==enymax
                            return;
                        end

                        imthr(y,x) = imthr(y,x) + 4;
                        nnew = nnew+1;

                        % expanding bounding box
                        if stx>x, stx=x; end
                        if enx<x, enx=x; end
                        if sty>y, sty=y; end
                        if eny<y, eny=y; end

                    end
                end
            end
        end

        if nnew == 0, break; end
    end

    % reset unconnected pixels within bounding box
    for y=sty:d:eny
        for x=stx:d:enx
            if bitand(imthr(y,x),4)==0
                imthr(y:y+d-1, x:x+d-1) = 0;
            else
                imthr(y, x) = imthr(y, x) - bitand(imthr(y,x),6); % remove marks
            end
        end
    end

    % increasing area slightly, as there might be small number of missed pixels at the very edge
    stx = stx-2; enx = enx+d+2;
    sty = sty-2; eny = eny+d+2;

    y=sty:eny;
    x=stx:enx;
    parea = sum(imthr(y,x), 'all');

    % center of mass
    xs = round(sum( int32(x) .* int32(sum(imthr(y,x),1) ) ) / sum(imthr(y,x),'all'), TieBreaker="plusinf");
    ys = round(sum( int32(y) .* int32(sum(imthr(y,x),2)') ) / sum(imthr(y,x),'all'), TieBreaker="plusinf");

    xm = repmat(x, length(y),1)-xs;
    ym = repmat(y',1,length(x))-ys;

    % collect coordinates of the pixels at the pupil edge, checking area around pupil center
    % do not include pixels having boundaries towards center xs ys
    edges = min(2, (~imthr(y,x-1)) + (~imthr(y,x+1)) + (~imthr(y-1,x)) + (~imthr(y+1,x))) - ...
        ((~imthr(y+1,x) | ~imthr(y,x-1)) .* (xm>=0 & ym<=0)) - ...
        ((~imthr(y-1,x) | ~imthr(y,x-1)) .* (xm>=0 & ym>=0)) - ...
        ((~imthr(y+1,x) | ~imthr(y,x+1)) .* (xm<=0 & ym<=0)) - ...
        ((~imthr(y-1,x) | ~imthr(y,x+1)) .* (xm<=0 & ym>=0));
    % subtracting 1 here has an effect of selecting only outer corner pixels
    edgepixel = imthr(y,x) & max(0,edges-1);

    [xp,yp] = find(edgepixel');
    xp = xp-xs+stx-1;
    yp = yp-ys+sty-1;
end
