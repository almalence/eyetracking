% -------------------------------------------------------------------------
%
%    Compute errors in iris ratio given the pupil edge points PCA angle
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% ToDo: for off-axis camera, ratio at the central position will be not 1, and we may need to take direction-specific component of it
% Looks like it will be an issue for camera_rot ~= 0 situations only
function diff = diffFromPCA(pxy, cxy, gv_cam, ctr, img_ctr, iris_z, ax, ay, verbose)

    N = length(pxy);
    ratio = zeros(N,1);

    for i=1:N

        rect_true_n = 0;
        rect_obs_n = 0;

        for f=1:length(pxy{i})
            xy = pxy{i}{f};
            A = EllipseDirectFit(xy);
            if ~isreal(A), diff=1e10; return; end
            [~, semx, semy, phi] = EllipseGeomFromAlgebraic(A);

            % find enclosing rectangle for the observed ellipse
            [rw, rh] = encl_rect(semx, semy, phi);
            % observed bounding rectangle ratio at given position
            rect_obs_n = rect_obs_n + rw/rh;

            % reverse-compensate PCA
            pca_mat = rotationMat([-ay -ax 0], 0);
            gv_pca = pca_mat * gv_cam(i,:)';

            % reverse-compensate refraction and camera perspective elongation of the true ellipse ratio
            cra_proj = camRayProjToMinor(cxy{i}{f}-img_ctr, cxy{i}{f}-ctr, iris_z);
            el_r = cos((acos(gv_pca(3))+cra_proj)/1.121)/cos(cra_proj);
            phi = atan2(gv_pca(2), gv_pca(1))+pi/2;
            [rw, rh] = encl_rect(1, el_r, phi);
            % expected bounding rectangle ratio for given position
            rect_true_n = rect_true_n + rw/rh;
        end

        ratio(i) = rect_obs_n / rect_true_n;
    end

    % goal is to reach the same relation to true ratio at all positions
    diff = sum((ratio-mean(ratio)).^2);
end


function [rw, rh] = encl_rect(semx, semy, phi)
    ux = semx * cos(phi);
    uy = semx * sin(phi);
    vx = semy * cos(phi + pi/2);
    vy = semy * sin(phi + pi/2);
    rw = 2*sqrt(ux*ux + vx*vx);
    rh = 2*sqrt(uy*uy + vy*vy);
end

