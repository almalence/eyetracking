% -------------------------------------------------------------------------
%
%    Estimate iris pose on 2d image plane
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ir_r, ir_r1, ir_r2, peak, iris2d, phi_minor] = getIris(ei, ir_r1, ir_r2, peak1, peak2, iris2d, iris_recenter, decent_v1, decent_v2,  el_ratio, phi_minor)

    % initialize from previous frame values
    ir_r = ei.iris_r;
    peak = ei.iris_peak;

    % if iris detected - update iris radius
    if ir_r1>0 || ir_r2>0
    
        % reject much weaker edge altogether
        if peak1>4*peak2, peak2=0; ir_r2=0;end
        if peak2>4*peak1, peak1=0; ir_r1=0;end
        % reject high difference in radii (>40% = 0.95mm pupil decenter)
        if min(ir_r1, ir_r2)*1.4 < max(ir_r1, ir_r2)
            if peak1<peak2, ir_r1=0;
            else, ir_r2=0; end
        end

        if ir_r1>0 && ir_r2>0
            ir_r_new = (ir_r1*peak1+ir_r2*peak2)/(peak1+peak2);
            peak_new = (peak1+peak2)/2;
        elseif ir_r1>0
            ir_r_new = ir_r1;
            peak_new = peak1;
        else
            ir_r_new = ir_r2;
            peak_new = peak2;
        end

        % update iris radius, weighting current peaks vs average observed peaks
        if peak, peak_new = min(peak_new, peak*4); end
        ir_r = ((2^ei.smooth_iris-1)*ir_r*peak + ir_r_new*peak_new) / ((2^ei.smooth_iris-1)*peak+peak_new);
        peak = ((2^ei.smooth_iris-1)*peak + peak_new) / (2^ei.smooth_iris);
    else
        % if no iris detection in current frame - compute expected iris radius from new distance from eyeball center
        d_prv_mm = norm(ei.eyeball3d_cam(1:2) - ei.iris2d_prev * ei.pxl2mm);
        d_new_mm = norm(ei.eyeball3d_cam(1:2) - iris2d * ei.pxl2mm);
        % using larger iris to eyeball-center distance to avoid error magnification in z
        % due to noisy x/y measurements when angle from camera is high
        iris2eye = max(ei.iris2hrot_mm, ei.iris2vrot_mm);
        z_prv = ei.iris_z + ( iris2eye - sqrt(max(0, iris2eye^2 - d_prv_mm^2)) ) / ei.pxl2mm;
        z_new = ei.iris_z + ( iris2eye - sqrt(max(0, iris2eye^2 - d_new_mm^2)) ) / ei.pxl2mm;
        ir_r = ir_r * z_prv/z_new;
    end

    % good detection - both sides detected with good peaks and with reasonably close values
    % re-position horizontal iris center from detected iris curves
    % Note: effectively negates coordinate corrections in TruePupilFromRefraction() from pupil refraction and decenter
    if ir_r1>0 && ir_r2>0
        if iris_recenter >= 1
            iris2d(1) = iris2d(1) + ellipseScales(0, phi_minor, el_ratio) * (ir_r2-ir_r1)/2;
        end
        if iris_recenter >= 2
            iris2d(2) = iris2d(2) + ellipseScales(0, phi_minor+pi/2, el_ratio) * (decent_v1+decent_v2)/2;
        end
        
        % rotate phi_minor according to inclination of iris edges
        if iris_recenter >= 3
            phi_minor = phi_minor + asin((decent_v2-decent_v1)/ir_r);
        end
    end

end
