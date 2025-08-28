% -------------------------------------------------------------------------
%
%    Eye pose estimation
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% Eye pose estimation is from a single frame, except for the following pieces where information from multiple combined:
%   - 2d eyeball center location averaging from ei.N recent frames in Estimate2dEyeCenter: controlled with  ei.smooth_2deye, ei.smooth_2dtail
%   - 3d eyeball location filtering in FilterEyeballPosition: ei.smooth_3deye
%   - smoothing of detected iris radius in pixels: ei.smooth_iris
% Set above control coefficients in eye config to 0 to disable correspondent smoothing or adjustment
%
function ei = getEyePose(im, ei)

    verbose     = ei.verbose;
    min_iris_r  = ei.min_iris_r;
    max_iris_r  = ei.max_iris_r;
    hmn_iris_mm = ei.hmn_iris_mm;
    iris2eye_mm = ei.iris2eye_mm;
    iris_z      = ei.iris_z;
    camera_rot  = ei.camera_rot;
    rotrad      = camera_rot * pi/180;

    [h,w] = size(im);
    img_ctr = [w/2 h/2]+1; % +1 to correct for Matlab 1-based indexes

    % ------------ Find pupil center location on the camera image

    % minimum pupil radius in pixels: 1.5mm pupil diameter, farthest distance from the camera
    min_pup_r = min_iris_r*1.5/hmn_iris_mm;
    
    % maximum pupil radius in pixels: 8mm pupil diameter, closest distance to the camera
    max_pup_r = max_iris_r*8/hmn_iris_mm;

    % lashes / glints / eyewear reflections removal and brightening, with 2x subsampling to speed-up processing
    imbr2x = cornerBrighten2x(removeSpeckles2x(im), ei.valid_area);

    [imthr2x, pupthr] = thresholdForPupil(imbr2x, min_pup_r/2, ei.pup_thr, ei.valid_area/2);

    % find the biggest black area via horizontal and vertical scans
    [xs, ys, psr] = getBigBlackArea(imbr2x, imthr2x, ei.valid_area/2);

    if xs<=0 || ys<=0, ei.blink = 1; return; end      % no pupil found (either eyelids closed or out of valid area)

    if ei.verbose>2
        figure(2); imshow(imthr2x); hold on; plot(xs/2,ys/2,'*'); hold off; title ('initial thresholding');
    end

    % re-threshold pupil using local levels around initial detection, at full resolution
    % useful in contrast-reduced situations, eg reflections from glasses at oblique angle
    imbr = cornerBrighten(removeSpeckles(im), ei.valid_area);
    imthr = reThreshold(imbr, xs, ys, psr, pupthr, max_pup_r);


    % ------------ Fit pupil to ellipse

    [xs, ys, x, y, parea] = extractEdgePixels(imthr, xs, ys, psr, max_pup_r);

    if length(x) < 8, ei.blink = 1; return; end                % require at least 8 pixels to try ellipse fitting

    if ei.verbose>2
        figure(3); imshow(imthr); hold on; plot(x+xs,y+ys,'.'); hold off; title('final threshold and ellipse edge pixels');
    end

    % clean the pixel coord data by eliminating pixels that are 'internal' to the radially-distant pixels
    [x, y] = cleanInternals(x,y, parea);

    if ei.verbose>2
        figure(3); hold on; plot(xs,ys,'*g', x+xs,y+ys,'og'); hold off; title('final and cleaned ellipse edge pixels');
    end

    % stop here and return previously estimated value if no pupil detection or very low confidence in detection
    if length(x) < 8, ei.blink = 1; return; end                % require at least 8 pixels to try ellipse fitting

    % smart fit of pupil to ellipse
    [A, usedpx] = pupilEllipse(x, y, min_pup_r, max_pup_r);

    % stop here and return previously estimated value if no ellipse fit
    if isempty(A),    ei.blink = 1; return; end

    % negate pupil non-circularity
    [cxy, semx, semy] = EllipseGeomFromAlgebraic(A);
    ei.el_ratio0 = min(semx,semy)/max(semx,semy);       % pupil ellipse ratio before any corrections (used in calibration)
    xy = [x(usedpx) y(usedpx)]-cxy;
    ei.pxy = xy; % used in calibration
    xy = compensatePupilElongation(xy, ei.el_ratio0, ei.pupil_angle*pi/180, ei.pupil_stretch);
    A = EllipseDirectFit(xy);

    % stop here and return previously estimated value if no ellipse fit after non-circularity corrections
    if isempty(A),    ei.blink = 1; return; end

    [~, semx, semy, phi] = EllipseGeomFromAlgebraic(A);
    cxy = cxy + [xs ys];

    % edge detection with ellipse fitting cause undershoot of major axis for small pupils
    % subtracting 1 pixel from each side makes result closer to ground truth
    % without affecting large or high-ratio pupils
    semx = max(1, semx-2);
    semy = max(1, semy-2);

    drawEyeImage(im, verbose);
    drawEyeCenter(cxy, [], 'g', [0 0], 0, verbose);
    drawPupil(x+xs, y+ys, xs, ys, usedpx, cxy(1), cxy(2), semx, semy, phi, uint8(imthr)*64+im/4*3, 1, 0, verbose);

    % figure out how much of pupil is visible and set blink ratio from that
    peli_area = pi*semx*semy;
    blink = max(0, 1 - 1.1 * parea / peli_area);

    if ei.smooth_2deye>0    % only react to sudden changes if 2d smoothing used
        % detecting sudden change in pupil area, and if such - discard this frame
        % sudden enlargement
        if peli_area > 2*ei.peli_area && blink>0.5
            ei.peli_area = (15*ei.peli_area + clip(peli_area, ei.peli_area/2, ei.peli_area*2))/16;
            return;
        end
    
        % sudden contraction
        if ei.peli_area > 2*peli_area && ei.blink<0.5
            ei.peli_area = (ei.peli_area + peli_area)/2;
            return;
        end
    end

    % maintain smoothed pupil area in eye instance
    if ei.blink>0.8
        ei.peli_area = (ei.peli_area + peli_area)/2;
    else
        ei.peli_area = (15*ei.peli_area + clip(peli_area, ei.peli_area/2, ei.peli_area*2))/16;
    end
    ei.blink = blink;

 
    % ------------ Correct pupil position and minor axis due to refraction

    % eyeball center needed to disambiguate position correction direction, but center only discovered after iris detection
    % so, use eyeball center, iris radius and iris ratio prediction from the previous frame
    eyeball2d_unrot = rotcam(ei.eyeball2d, -camera_rot) + img_ctr;

    [pc, el_major, el_minor, phi_minor] = TruePupilFromRefraction(cxy, semx, semy, phi, ei.pxl2mm, ...
        eyeball2d_unrot, ei.pupil_refr, ei.pca_angle, ei.iris_ratio, iris_z, img_ctr, ei.iris2eye_mm - ei.pupil2eye_mm, ei.max_el_ratio);

    if ei.el_known % if calibrating - use supplied ground truth phi_minor and el_ratio
        phi_minor = ei.phi_minor;
        el_minor = el_major * ei.el_ratio;
    else
        ei.phi_minor = phi_minor; % used in calibration
        % ellipse main axis is the radius of the pupil/iris (it is perpindicular to the camera ray pointing to the iris center)
        ei.el_ratio  = el_minor/el_major;
    end
    ei.pupil2d = rotcam(pc-img_ctr, camera_rot);
    drawPupil([], [], [], [], [], pc(1), pc(2), el_major, el_minor, phi_minor, [], 0, 1, verbose);

    % initial estimation on iris center
    ixy = correctPupilDecenter(pc, ei.el_ratio, phi_minor, ei.pxl2mm, ei.pupil_decenter_mm, camera_rot);


    % ------------ Rasterize iris

    ei.iris2d_prev = ei.iris2d;     % keep iris position estimation from previous frame (may be used for iris radius estimation if no iris detection in current frame)
    ei.iris2d_p = rotcam(ixy-img_ctr, camera_rot);  % rough iris center location, estimated from pupil center

    % look to the left / right of the pupil center, at up to [-15 +30] degree from horizontal line,
    % distances to look at are limited by possible Z distances of the eye
    
    % fill the list of radial 'scalers' computed from ellipse shape which will turn iris shape into circle
    iris_sctr_up = pi/12; % angle from horizontal up to which iris is usually visible  15 degree up
    iris_sctr_dn = pi/6;  %                                                            30 degree down
    iris_sctr = iris_sctr_up + iris_sctr_dn;
    iris_d1 = linspace(pi+iris_sctr_up, pi-iris_sctr_dn, min_iris_r) - rotrad;
    iris_d2 = linspace(  -iris_sctr_up,    iris_sctr_dn, min_iris_r) - rotrad;
    iris_scl1 = ellipseScales(iris_d1, phi_minor, ei.el_ratio);
    iris_scl2 = ellipseScales(iris_d2, phi_minor, ei.el_ratio);

    % rasterize iris sectors in pseudo-polar coordinates
    % use 2x subsampled imbr for better processing speed
    % at horizontal iris angles above roughly 38 degree (iris_scl=0.79) the iris edge farthest from the camera is occluded by corneal bulge and should not be used
    % at such angles the farther edge is also significantly farther in Z distance and iris radius readings would require correction according to camera ray angle
    ir_r1=0; peak1=0; decent_v1=0;
    ir_r2=0; peak2=0; decent_v2=0;
    if ei.iris2d_p(1)>=ei.eyeball_ctr(1) || iris_scl1(floor(end*iris_sctr_up/iris_sctr))>0.79
        [ir_r1, peak1, decent_v1] = detectIrisRadius( rasterizeIris(iris_d1, iris_scl1, imbr2x, ixy, min_iris_r, max_iris_r, iris_z, phi_minor, camera_rot), min_iris_r, el_major, iris_sctr);
    end
    if ei.iris2d_p(1)<=ei.eyeball_ctr(1) || iris_scl2(floor(end*iris_sctr_up/iris_sctr))>0.79
        [ir_r2, peak2, decent_v2] = detectIrisRadius( rasterizeIris(iris_d2, iris_scl2, imbr2x, ixy, min_iris_r, max_iris_r, iris_z, phi_minor, camera_rot), min_iris_r, el_major, iris_sctr);
    end

    drawIris(min_iris_r, iris_d1, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, 3, verbose, 'w', 'w');
    drawIris(max_iris_r, iris_d1, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, 3, verbose, 'w', 'w');
    drawIris(ir_r1, iris_d1, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, (ir_r1>0)+1, verbose, 'b', 'r');
    drawIris(min_iris_r, iris_d2, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, 3, verbose, 'w', 'w');
    drawIris(max_iris_r, iris_d2, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, 3, verbose, 'w', 'w');
    drawIris(ir_r2, iris_d2, ixy, ei.el_ratio, phi_minor, iris_z, img_ctr, (ir_r2>0)+1, verbose, 'b', 'r');


    % ------------ Find iris 2d pose, in camera pixels

    % starting from here, switching from frame pixel coordinates to coordinate system with matched horizontal/vertical eye directions
    phi_minor = phi_minor+rotrad;
    % get combined iris radius from two measurements
    % adjusting refraction-/decenter-corrected iris2d and phi_minor direction from iris edege curves
    [ei.iris_r, ei.iris_r1, ei.iris_r2, ei.iris_peak, ei.iris2d, phi_minor] = ...
        getIris(ei, ir_r1, ir_r2, single(peak1), single(peak2), ei.iris2d_p, ei.iris_recenter, decent_v1, decent_v2,  ei.el_ratio, phi_minor);

    ei.cxy = rotcam(ei.iris2d, -camera_rot) + img_ctr;    % used in calibration

    
    % ------------ From above - pixels-to-millimeters ratio and 3d iris center location in camera space

    % if iris not detected - do not compute pxl2mm, use previous result as a stub
    if ei.iris_r
        ei.pxl2mm = hmn_iris_mm/2/ei.iris_r;
    end
    iris2eye_pxl = iris2eye_mm/ei.pxl2mm;
    ei.iris3d_cam = ei.pxl2mm * [ei.iris2d iris_z];


    % ------------ From (mulitple) ellipses - estimate center around which they are moving (= 2d eyeball center)

    if ei.smooth_2deye < 0     % negative smoothing means calibration pass, and eyeball center is fixed and given
        eyeball_confid = 1;
        % derive center position from position of rotation axes
        ei.eyeball2d = [ei.iris2d(1)+(ei.ctr_tail(1)-ei.iris2d(1))*ei.iris2vrot_mm/ei.iris2hrot_mm ei.ctr_tail(2)];
    else
        if mod(ei.frame, ei.fps30)  % no need to estimate eyeball center more frequently than 30 times a second
            eyeball_confid = ei.eyeball_confidence;
        else
            [ei, exy, eyeball_confid] = Estimate2dEyeCenter(ei, ei.el_ratio, phi_minor, ei.eyeball2d, iris2eye_pxl, img_ctr, verbose);
            ei.eyeball2d = exy;
        end
    end

    % ------------ From 2d eyeball center - 3d eyeball center location in camera 3d space.

    % use previous location if it will be impossible to compute new one below
    stabilized_eyeball3d_cam = ei.eyeball3d_cam;

    if ei.iris_r && eyeball_confid>0
        eyeball3d_cam = EyeballCenter3d(iris_z, ei.eyeball3d_cam(3), ei.iris2d, ei.eyeball2d, iris2eye_pxl, ei.pxl2mm);
        [stabilized_eyeball3d_cam, ei.eyeball3d0x] = FilterEyeballPosition(eyeball3d_cam, ei.eyeball3d0x, ei.iris3d_cam, iris2eye_mm, ei.iris2hrot_mm, ei.smooth_3deye, ei.eyeball_confidence, eyeball_confid);
    end


    % ------------ Recompute 3d iris Z from stabilized 3d eyeball center and 2d iris location

    iris2eye_2d = norm(ei.pxl2mm*ei.iris2d - stabilized_eyeball3d_cam(1:2));
    if iris2eye_mm <= iris2eye_2d
        % if impossible to compute 3d iris location from stabilized eyeball center (too far) -
        % scale direction towards stabilized center to position eyeball center at correct distance from iris
        eyeball_vec = stabilized_eyeball3d_cam - ei.iris3d_cam;
        stabilized_eyeball3d_cam = ei.iris3d_cam + eyeball_vec * (iris2eye_mm/norm(eyeball_vec));
        iris2eye_2d = min( iris2eye_mm, norm(ei.pxl2mm*ei.iris2d - stabilized_eyeball3d_cam(1:2)) );
    end
    iris_dz = sqrt(iris2eye_mm^2 - iris2eye_2d^2);
    ei.iris3d_cam(3) = stabilized_eyeball3d_cam(3)-iris_dz;

    ei.eyeball3d_cam = stabilized_eyeball3d_cam;
    ei.eyeball_confidence = eyeball_confid;

    % estimate stabilized iris ellipse ratio (used in TruePupilFromRefraction() for the next frame)
    iris_ratio = dot([0 0 1], (ei.eyeball3d_cam-ei.iris3d_cam) / iris2eye_mm);
    cxy = cxy-img_ctr;
    if ei.pxl2mm * norm(cxy - ei.iris2d_s) > 2   % less smoothing if iris has moved far away (2mm displacement, possibly accumulated over multiple frames)
        ei.iris_ratio = (ei.iris_ratio + iris_ratio) / 2;
        ei.iris2d_s = (ei.iris2d_s+cxy)/2;
    else
        ei.iris_ratio = ( 15*ei.iris_ratio + iris_ratio ) / 16;
        ei.iris2d_s = (15*ei.iris2d_s+cxy)/16;
    end


    % ------------ Linear transformation of iris/eyeball position into HMD coordinate space.

    ei.iris3d = ei.cam2hmd(:,1:3) * ei.iris3d_cam' + ei.cam2hmd(:,4);
    ei.eyeball3d = ei.cam2hmd(:,1:3) * ei.eyeball3d_cam' + ei.cam2hmd(:,4);


    % ------------ Gaze direction calculated from eyeball center vs iris center 3d locations.

    optaxis_vector = (ei.eyeball3d-ei.iris3d) / iris2eye_mm;
    % rotate gaze direction by alpha angle matrix
    ei.gaze_vector = ei.alpha_mat * optaxis_vector;

    % simple filter to reduce gaze direction jitter (no filtering if angle difference is above 2 degree)
    ei.gaze_vector_smooth = smoothGazeVector(ei.gaze_vector, ei.gaze_vector_smooth, ei.gv_smooth_diff, ei.fps);

    % record angle from HMD center in observations
    if ei.last_bin>=0
        ei.ang_obs(ei.last_bin) = atan2(ei.gaze_vector(2), ei.gaze_vector(1));
    end
    ei.frame = ei.frame+1;

end
