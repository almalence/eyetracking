% -------------------------------------------------------------------------
%
%    Plot pupil ellipse
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function drawPupil(x, y, xs, ys, usedpx, cx, cy, semx, semy, phi, imthr, drawleft, drawright, verbose)

    if verbose==0, return; end

    % draw the ellipse and its axes
    R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ];
    ver_line  = R * [ [0 0]; semy*[-1 1] ];
    horz_line = R * [ semx*[-1 1]; [0 0] ];
    theta_r   = linspace(0,2*pi);
    ellipse   = R * [semx*cos(theta_r); semy*sin(theta_r)];
    
    % draw
    figure(1);
    
    if drawleft
        subplot(1,2,1);
        if verbose>1
            imshow(imthr);
            hold on;
            plot(xs,ys,'*');
            plot(x,y,'.b');
            plot(x(usedpx),y(usedpx),'og');
            plot( cx+ver_line(1,:), cy+ver_line(2,:),'r' );
            plot( cx+horz_line(1,:), cy+horz_line(2,:),'r' );
            plot( cx+ellipse(1,:), cy+ellipse(2,:),'r' );
            hold off;
        else
            [h,w] = size(imthr);
            axis image ij;
            axis([1 w 1 h]);
        end
    end

    if drawright
        subplot(1,2,2);
        hold on;
        plot( cx+ellipse(1,:), cy+ellipse(2,:),'r' );
        plot(cx,cy,'.r');
        hold off;
    end

end
