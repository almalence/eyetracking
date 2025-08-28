% -------------------------------------------------------------------------
%
%    Remove inward-placed pixels
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [x, y] = cleanInternals(x,y, parea)

    % approximation of pupil diameter
    pup_dia = 2*sqrt(parea/pi);

    % sort in angular order
    ang = atan2(y, x);
    [~, idx] = sort(ang);
    xy = [x(idx) y(idx)];

    while true
        repeat = false;
        n = size(xy,1);
        valid = ones(n,1,'logical');
        p1 = n;
        for i=1:n
            if i<n, p2 = i+1;
            else
                p2 = find(valid, 1);
                if p2==n, break; end
            end

            vc = xy(i,:);
            v1 = xy(p1,:)-vc;
            v2 = xy(p2,:)-vc;
            a1 = acos( dot(-vc, v1) / (norm(v1)*norm(-vc)) );
            a2 = acos( dot(-vc, v2) / (norm(v2)*norm(-vc)) );
    
            % the closer the distance between points to the approximate
            % pupil diameter, the smaller the largest acceptable angle
            max_ang = 177 - 25 * clip(max(norm(v1),norm(v2))/pup_dia - 0.25, 0, 1);
    
            if a1+a2 > max_ang*pi/180
                valid(i) = 0;
                % if current is invalid - the previous (p1) point may become invalid as well
                repeat = true;
            else
                p1 = i;
            end
        end
        xy = xy(valid,:);
        if repeat == false, break; end
    end

    x = xy(:,1);
    y = xy(:,2);
end
