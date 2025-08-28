% -------------------------------------------------------------------------
%
%    Measure eye pose precision and accuracy according to Tobii white paper recommendations
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [prec, accuracy, ranges] = precisionAndAccuracyPerAngleRange(deg, deg_gt, vid)

    % smooth pursuit - divide into small segments to properly measure precision
    if all(vid==0)
        vid = floor((1:length(vid)) / 16);
    end

    % compute per-stimulus precision
    deg_mean = zeros(size(deg));
    prec = zeros(size(deg,1),1);
    for i=0:max(vid)
        seg = vid==i;
        if ~isempty(find(seg,1))
            deg_mean(seg,1) = mean(deg(seg,1)-deg_gt(seg,1)) + deg_gt(seg,1);
            deg_mean(seg,2) = mean(deg(seg,2)-deg_gt(seg,2)) + deg_gt(seg,2);

            devi = vecnorm( deg(seg,:)-deg_mean(seg,:) ,2,2);
            % precision = RMSD
            prec(seg) = sqrt(mean(devi.^2));
        end
    end

    % difference between measured and ground truth
    diff = abs(deg-deg_gt);
    
    % vecnorm converts horizontal/vertical angles into radial values (assumption of reasonably small angles)
    deg_gt_r = vecnorm(deg_gt, 2,2);
    accuracy = vecnorm(diff, 2,2);

    ranges = zeros(size(deg,1),1);
    ranges(deg_gt_r<10) = 1;
    ranges(deg_gt_r>=10 & deg_gt_r<20) = 2;
    ranges(deg_gt_r>=20 & deg_gt_r<25) = 3;
    ranges(deg_gt_r>=25) = 4;
end

