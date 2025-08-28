% -------------------------------------------------------------------------
%
%    Compute errors in eyeball center estimation given the pupil parameters (pupil elongation and PCA)
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [diff, ctrn] = ctrDiffFromPupilRatio(pxy, cxy, el_ratio0, iris_ratio, ctr, img_ctr, pxl2mm, max_el_ratio, iris_z, iris2eye_mm, iris2hrot_mm, iris2vrot_mm, v, verbose)

    pup_ang = v(1);
    pup_stretch = v(2);
    pca_angle = [v(3) v(4)];
    pup_refr = v(5);
    
    N = length(pxy);

    diff = 0;
    ctrn = zeros(N,2);
    for i=1:N
        for f=1:length(pxy{i})
            xy = pxy{i}{f};
            % pupil elongation compensation
            xy = compensatePupilElongation(xy, el_ratio0{i}{f}, pup_ang, pup_stretch);
            A = EllipseDirectFit(xy);
            if ~isreal(A), diff=1e10; return; end
            [~, semx, semy, phi] = EllipseGeomFromAlgebraic(A);
            el_r = min(semx,semy) / max(semx,semy);

            phi_minor = phi + (semx<semy)*pi/2;

            % apply refraction/perspective correction
            cra_proj = camRayProjToMinor(cxy{i}{f}-img_ctr, cxy{i}{f}-ctr, iris_z);
            el_r = cos(acos(el_r*cos(cra_proj))*pup_refr - cra_proj);
            % pupillary circular axis angle compensation
            vec_to_eyeb = ctr-cxy{i}{f};
            [el_r, phi_minor] = compensatePCA(vec_to_eyeb, iris_ratio{i}{f}, el_r, phi_minor, pca_angle, max_el_ratio);
            % phi_mior points toward 2d eyeball center, range: [-pi/2 .. 3*pi/2]

            % Note: unlike in Estimate2dEyeCenter(), estimation is in unrotated camera pixel coordinates here
            dv = [sin(phi_minor) -cos(phi_minor)];
            proj_ang = acos(el_r);
            len2dy = iris2eye_mm/pxl2mm(i)*sin(proj_ang);
            len2dx = iris2hrot_mm/iris2vrot_mm * len2dy;
            e1h = cxy{i}{f} + len2dx.*dv;
            e1v = cxy{i}{f} + len2dy.*dv;
            ctr_from_el_r = [e1h(1) e1v(2)];

            ctrn(i,:) = ctrn(i,:) + ctr_from_el_r;
        end

        ctrn(i,:) = ctrn(i,:)/length(pxy{i});
        diff = diff + norm(ctrn(i,:) - ctr);

        if verbose, fprintf('%6.1f %6.1f   ', ctrn(i,:)); end
    end

    if verbose, fprintf('\n'); end
end
