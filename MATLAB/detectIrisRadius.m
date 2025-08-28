% -------------------------------------------------------------------------
%
%    Simple method of iris radius detection: find maximum contrast break on pseudo-polar image of iris
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% The line can be slanted up to ~1/15th of the raster width (with -pi/6..pi/6 iris raster scan)
% due to up to 0.9mm decenter of pupil relative to the iris
% (assuming raster width about equals iris radius =0.9*sin(pi/6)/iris_r_mm)
% This smears the gradient peak a bit but still provides good enough peak amplitude
% 
function [ir_r, irpeak, decent_v] = detectIrisRadius(irim, min_iris_r, el_major, iris_sctr)

    [h,w] = size(irim);

    pupil_edge = ceil(max(0, el_major - min_iris_r));

    % eliminate glints
    meanv = sum(irim);
    thrim = irim>=floor(meanv*5/(4*h));
    meanvim = repmat(floor(meanv/h),h,1);
    irim(thrim) = meanvim(thrim);

    % contrast window length
    tlen = ceil(w/10);

    % slant step (to detect vertical decentration)
    slstep = ceil(w/15);

    brtl = zeros(w,1, 'int32');
    brtr = zeros(w,1, 'int32');
    devr = zeros(w,1, 'int32');

    % int32 here to match C version
    meanv = int32(sum(irim));
    devv  = idivide(int32(sum(abs(int32(irim)*h-meanv))), h);

    % sliding window
    for i=tlen+pupil_edge:w-tlen
        brtl(i) = sum(meanv(i-tlen+1:i));
        brtr(i) = sum(meanv(i+1:i+tlen));
        devr(i) = sum(devv(i+1:i+tlen));
    end

    [peak, ir_r] = max( 2*brtr-devr - 2*brtl );
    % looking for at least 10% diff in brightness (additional divider of 10 is from tlen=w/10)
    irpeak = max(0, peak-sum(meanv)/100);

    ir_r = ir_r-1;  % correction for Matlab 1-based indexes

    if ir_r<=tlen+pupil_edge || ir_r>=w-tlen-1 || irpeak==0
        ir_r=0;
        irpeak=0;
        decent_v=0;
    else
        % estimate vertical pupil decentration from rasterized iris edge slant
        peak_slant = [0 0];
        for s=1:2
            slant = ((s-1)*2-1)*slstep;

            irsl = irim;
            for i=1:h
                dx = round(slant*(i-1-h/2)/(h/2), TieBreaker="plusinf");
                veci = max(1,dx+1):min(w,w+dx);
                veco = max(1,-dx+1):min(w,w-dx);
                irsl(i,veco) = irim(i,veci);
            end
            meanv = int32(sum(irsl));
            devv  = idivide(int32(sum(abs(int32(irsl)*h-meanv))), h);

            brtl = sum(meanv(ir_r-tlen+1:ir_r), 'native');
            brtr = sum(meanv(ir_r+1:ir_r+tlen), 'native');
            devr = sum(devv(ir_r+1:ir_r+tlen), 'native');

            peak_slant(s) = 2*brtr-devr - 2*brtl;
        end

        % decenter here is measured in horizontal iris raster pixels,
        % which equals camera pixels along major ellipse axis
        % vertical decenter computed with formula below is giving somewhat reduced values for some reason
        decent_v = slstep * clip( quad_fit(peak_slant(1), peak, peak_slant(2)), -1, 1) / sin(iris_sctr/2);
        ir_r = ir_r + min_iris_r;
    end

end

function x = quad_fit(ym1,y0,y1)

    ym1 = single(ym1);
    y0  = single(y0);
    y1  = single(y1);
    
    % find parabola extremum (where derivative of y=ax^2+bx+c goes to zero)
    if y0>(ym1+y1)/2
        x = (ym1 - y1) / (2*ym1 + 2*y1 - 4*y0);
    else
        % if parabola seem to be inverted (or just straight line) - return location of maximum
        if y1>ym1, x=1; else x=-1; end
    end
end
