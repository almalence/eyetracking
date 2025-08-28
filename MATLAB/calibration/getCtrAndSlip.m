% -------------------------------------------------------------------------
%
%    Find eyeball center on image plane and slippage of each location
%    Note: camera perspective projection not accounted for here
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ctr, slippage, f] = getCtrAndSlip(p2d, gt, N, verbose, clr)

    ctr_gt = mean(gt);
    ctr_p2d = mean(p2d);
    scl_gt = sum(abs(gt-ctr_gt));
    scl_p2d = sum(abs(p2d-ctr_p2d));

    p2d_ideal = ctr_p2d + (gt-ctr_gt) .* scl_p2d ./ scl_gt;

    slippage = p2d_ideal-p2d;

    ctrp = estimateCtr(p2d_ideal, gt, N);
    ctr = mean(ctrp);

    f = 1;
    if verbose
        f = figure; hold on;
        for i=1:N
            plot(ctrp(i,1), ctrp(i,2), 'x', 'Color',clr{i});
            plot(p2d(i,1), p2d(i,2), '+', 'Color',clr{i});
            plot(p2d_ideal(i,1), p2d_ideal(i,2), 'o', 'Color',clr{i});
            plot([ctrp(i,1) p2d_ideal(i,1)], [ctrp(i,2) p2d_ideal(i,2)], ':', 'Color', clr{i});
        end
        c10 = round([ctr; p2d_ideal]/10)*10;
        hold off; axis image; axis ij; xlim([min(c10(:,1))-20 max(c10(:,1))+20]); ylim([min(c10(:,2))-20 max(c10(:,2))+20]); grid on; grid minor;
    end

end


function ctrp = estimateCtr(p2d, gt, N)
    ctrp = zeros(N,2);
    for i=1:N
        d_ref = [0 0];
        d_cam = [0 0];
        for j=1:N
            d_ref = d_ref + abs( gt(i,:) - gt(j,:) );
            d_cam = d_cam + abs( p2d(i,:) - p2d(j,:) );
        end
    
        ctrp(i,:) = p2d(i,:) - gt(i,1:2) .* d_cam ./ d_ref;
    end
end
