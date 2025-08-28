% -------------------------------------------------------------------------
%
%    Estimate eyeball center location on 2d image plane
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ei, exy, eyeball_confid] = Estimate2dEyeCenter(ei, iris_ratio, phi_minor, eyeball2d, iris2ctr_pxl, img_ctr, verbose)

    exy = eyeball2d;
    eyeball_confid = 0.1;
    quals = 1;

    % reduce weight of old observations: halving after 3 seconds
    ei.ctr_weight = ei.ctr_weight .* (1-0.2/ei.fps_ec2d);
    ei.obs_time = ei.obs_time + 1;
    
    % iris not detected or pupil not fully visible or too low angle to camera - do not compute / update eyeball center
    % 0.98=11.5 degree; 0.22=12.7 degree
    if ei.iris_r && ei.blink==0 && iris_ratio<0.98 && norm((ei.ctr_tail-ei.iris2d).*[ei.iris2vrot_mm/ei.iris2hrot_mm 1])/iris2ctr_pxl > 0.22
        % derive iris center to eyeball center distance in 2d space by computing trigonometry from el_minor/el_major
        % Note: this method has precision issues if single frame only is used:
        %  - low precision for low angles as angle precision from acos is very low if el_minor is close to el_major
        %  - low precision for very high angles due to refraction effects causing elongation of el_minor (correction done in TruePupilFromRefraction)
        below15deg = clip((0.98-iris_ratio)*20, 0.1, 1);
        above30deg = clip((iris_ratio-0.76)*10,  0.4, 1);
        occluded = max(0.1, (1-3*ei.blink));

        % origin and direction of 2d line from current observation
        dv = [sin(phi_minor) -cos(phi_minor)];     % direction vector

        proj_ang = acos(iris_ratio);
        len2dy = iris2ctr_pxl*sin(proj_ang);

        % horizontal rotation has axis farther then eyeball center
        len2dx = ei.iris2hrot_mm/ei.iris2vrot_mm * len2dy;

        % eyeball center in 2d is at the len2d distance from the iris center
        % in the direction of minor ellipse axis (one of two possible locations)
        % e1v/e2v - possible locations of eyeball center (coincides with point on rotation axis for the vertical gaze direction)
        iris2hrot = len2dx.*dv;
        iris2vrot = len2dy.*dv;
        e1v = ei.iris2d + iris2vrot;                    % eyeball center coordinates
        e2v = ei.iris2d - iris2vrot;
        e1  = ei.iris2d + [iris2hrot(1) iris2vrot(2)];  % rotation axes locations (x for horizontal rotation, y for vertical)
        e2  = ei.iris2d - [iris2hrot(1) iris2vrot(2)];

        drawEyeCenter(e1v, [], 'r', img_ctr, ei.camera_rot, verbose);
        drawEyeCenter(e2v, [], 'r', img_ctr, ei.camera_rot, verbose);

        % distance from designed 'perfect' center
        % considering the eyeball center farther than twice the iris to-eyeball distance (about 24mm) from designed position very unlikely
        w1d = max(0.1, min(1, 1 - max(0, norm(e1-ei.eyeball_ctr)/iris2ctr_pxl-1) ));
        w2d = max(0.1, min(1, 1 - max(0, norm(e2-ei.eyeball_ctr)/iris2ctr_pxl-1) ));
        % distance from previous center detection
        w1p = max(0.2, 1 - norm(e1v-eyeball2d)/iris2ctr_pxl );
        w2p = max(0.2, 1 - norm(e2v-eyeball2d)/iris2ctr_pxl );

        % if observation is of good quality - select which of the previous observations to replace
        % quality reduction from:
        % - low weights above
        % - highly occluded pupil (blink)
        % - high or low pupil angles
        w1 = w1p*w1d * occluded * below15deg * above30deg;
        w2 = w2p*w2d * occluded * below15deg * above30deg;
        qual = max(w1, w2);               % quality to save to the bin
        quals = occluded * below15deg * max(w1p, w2p) * ei.slippage_rate; % quality metric for smoothing

        % prevent same direction from filling all the bins by limiting range of bins to occupy according to radial direction
        N = ei.N;
        binrange = mod(round(phi_minor/(2*pi)*N + (-N/8:N/8)) + N, N) + 1;

        % find empty bin
        obin = binrange( find(~ei.valid(binrange),1) );
        if isempty(obin)
            % no empty bins - replace bin with the lowest quality and similar radial direction
            minqual = 0;
            for i=binrange
                qualbin = max(ei.ctr_weight(i,:));
                qualnew = qual * max(0, cos(phi_minor-ei.phi_obs(i))*2);
                if qualnew-qualbin > minqual
                    minqual = qualnew-qualbin;
                    obin = i;
                end
            end
        end

        if ~isempty(obin)
            ei.obs_time(obin)      = 0;
            ei.valid(obin)         = true;
            ei.phi_obs(obin)       = phi_minor;
            ei.ctr_obs(obin, :, :) = [e1; e2];
            ei.ctr_weight(obin,:)  = [w1; w2];
        end

        % construct averaged eyeball position
        % from a set of prior observations of possible eyeball center locations and current one

        % select cluster centers based on both current detection and previous smoothed one
        c1 = e1; c2 = e2;
        if sum(ei.valid) == N
            if norm(e1-ei.ctr_tail) < norm(e2-ei.ctr_tail)
                c1 = (e1+ei.ctr_tail*3)/4;
            else
                c2 = (e2+ei.ctr_tail*3)/4;
            end
        end

        % limit to fresh observations only if enough high-quality observations with wide enough angular spread
        % sort from fresh to old
        [~,bins] = sort(ei.obs_time);
        valid = zeros(N, 1, 'logical');
        bin = ei.last_bin;
        min_ang = 0;
        if bin >= 0, min_ang = ei.ang_obs(bin); end
        max_ang = min_ang;
        n = 0;
        for bin = bins'
            if ei.valid(bin)
                valid(bin) = true;
                if max(ei.ctr_weight(bin,:)) > 0.7       % limit to high-quality observations only
                    % check if new observation extends angular spread
                    if mod(ei.ang_obs(bin)-max_ang, 2*pi) < pi
                        max_ang = ei.ang_obs(bin);
                    end
                    if mod(min_ang-ei.ang_obs(bin), 2*pi) < pi
                        min_ang = ei.ang_obs(bin);
                    end
                    if n>=8 && mod(max_ang-min_ang, 2*pi) > 45*pi/180
                        break;
                    end
                end
                n = n+1;
            end
        end

        if ~isempty(obin), ei.last_bin = obin; end

        % divide into clusters by comparing distance to all previous observations
        % and find which cluster is more compact
        obs = ei.ctr_obs(valid,:,:);
        cmbnd_obs = prod(size(obs,[1 2]));
        obs = reshape(obs, [cmbnd_obs 2]);
        weight = reshape(ei.ctr_weight(valid,:), [cmbnd_obs 1]);
        d1 = vecnorm((c1-obs)').^2;
        d2 = vecnorm((c2-obs)').^2;
        % divide obervations into clusters
        near1 = d1<=mean(d1);
        near2 = d2<=mean(d2);

        % re-center clusters and measure cluster size as mean distance from it's center
        c1 = mean(obs(near1,:),1);
        c2 = mean(obs(near2,:),1);
        d1 = vecnorm((c1-obs)').^2;
        d2 = vecnorm((c2-obs)').^2;

        % discard outliers in a cluster
        meanD1near1 = mean(d1(near1));
        meanD2near2 = mean(d2(near2));
        if sum(near1)>=8
            near1 = d1 <= 2*meanD1near1;
            meanD1near1 = mean(d1(near1));
        end
        if sum(near2)>=8
            near2 = d2 <= 2*meanD2near2;
            meanD2near2 = mean(d2(near2));
        end

        % Selection from two possible eyeball center locations is a weighted product of:
        % - distance from designed 'perfect' center
        % - distance from previous center detection
        % - mean distance to nearby centers (equivalent of cluster radius) from previous frames
        %   considering avg cluster of 2mm or above very unlikely
        w1n = max(0.2, 1 - 4*sqrt(meanD1near1)/iris2ctr_pxl);
        w2n = max(0.2, 1 - 4*sqrt(meanD2near2)/iris2ctr_pxl);

        if w2d*w2p*w2n > w1d*w1p*w1n
            obs_near = obs(near2,:);
            wp = weight(near2);
            exy=e2v;
        else
            obs_near = obs(near1,:);
            wp = weight(near1);
            exy=e1v;
        end
        eyeball_confid = mean(wp);

        % initialize tail of accumulated observations to the first frame if no history
        if sum(valid)<=1
            ei.ctr_tail = [ei.iris2d(1) + (exy(1)-ei.iris2d(1))*ei.iris2hrot_mm/ei.iris2vrot_mm exy(2)];
        end

        % stabilization of 2d eyeball center by weighted-averaging with last readings
        % Note: this stabilization works best for on-axis camera
        % For off-axis camera the separate x/y stabilization is flawed
        % because eye h/v rotations projected on camera 2d plane are not straight lines anymore
        meanxy = [sum(obs_near(:,1).*wp)/sum(wp) sum(obs_near(:,2).*wp)/sum(wp)];

        % Update long-term observation tail (keeping about 2 sec)
        % readings of higher confidence provide stronger updates to the tail
        tail_upd = max(0.1, 3*(eyeball_confid^2));
        ei.ctr_tail = (ei.ctr_tail*(2^ei.smooth_2dtail) + meanxy*tail_upd) / (2^ei.smooth_2dtail + tail_upd);
    end

    % derive horizontal center position from position of horizontal rotation axis
    tail = [ei.iris2d(1)+(ei.ctr_tail(1)-ei.iris2d(1))*ei.iris2vrot_mm/ei.iris2hrot_mm  ei.ctr_tail(2)];
    % smooth, but no smoothing for high-quality detections
    exy = (tail*(2^ei.smooth_2deye-1)*(1-quals) + exy*quals) / ((2^ei.smooth_2deye-1)*(1-quals)+quals);
    
    drawEyeCenter(exy, ei.iris2d, 'm', img_ctr, ei.camera_rot, verbose);

end
