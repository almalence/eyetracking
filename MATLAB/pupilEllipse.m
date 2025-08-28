% -------------------------------------------------------------------------
%
%    Fit pupil edge points to ellipse via multiple 'sectoral' fits
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% assuming either a part of pupil edge is occulded or 2 parts of pupil edge occluded (not more)
function [A, usedpx] = pupilEllipse(x, y, min_pup_r, max_pup_r)

    A=[];
    
    % dividing circle into 8 sectors
    Nsec = 8;
    s_ang = atan2(y, x);
    abin = min(Nsec, 1 + floor((s_ang+pi)*Nsec/(2*pi)));

    maxgoodpx = 0;
    usedpx = [];

    % two sectors
    % Note: not quite covering all the possibilities as the circle start boundary never crossed
    for s1=1:Nsec-1
        for e1=s1:Nsec
            for s2=e1+2:Nsec-1
                for e2=s2:Nsec
                    % use at least half of the circle
                    if e2-s2+e1-s1<Nsec/2-1, continue; end

                    use1 = abin>=s1 & abin<=e1;
                    use2 = abin>=s2 & abin<=e2;
                    
                    use = use1 | use2;
                    % not enough total points
                    if sum(use)<8, continue; end
                    [ellipse, npx] = ellipseFitAndDistance(single(x), single(y), use, 2, min_pup_r, max_pup_r);

                    if npx>maxgoodpx, A = ellipse; maxgoodpx = npx; usedpx = use; end
                end
            end
        end
    end
end
