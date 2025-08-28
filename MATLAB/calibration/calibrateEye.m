% -------------------------------------------------------------------------
%
%    Calibration of eye parameters
% 
%    SPDX-FileCopyrightText: Copyright (c) 2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [ei, ec] = calibrateEye(ei, cam, anc, T2anc, verbose)

M = length(cam);    % total frames
N = size(anc,1);    % total calibration points
ec=[];

% calculate projection of reference gaze vectors onto camera image plane
gv_cam = zeros(N,3);
for i=1:N
    % reference gaze vector in HMD space
    gv = -sin(anc(i,:)*pi/180);
    gv(3) = sqrt(1-sum(gv.^2));
    gv_cam(i,:) = (ei.hmd2cam(:,1:3) * gv')';
end
% location closest to the camera (most direct look to the camera)
[~, primary_anc] = max(gv_cam(:,3));

% select which frames to use for measurements
Tancuse = cell(1,N);
for i=1:N
    % use only last half of the position (eye is better stabilized)
    len = sum(T2anc==i);
    nouse = floor(len/2);
    Tancuse{i} = [zeros(1,nouse, 'logical') ones(1,len-nouse, 'logical')];
end
Tuse = logical([Tancuse{:}]);
% Tuse_primary = Tuse&0;
% Tuse_primary(T2anc==primary_anc) = Tancuse{primary_anc};

% reset eye instance calibrations
ei.alpha_mat = eye(3);

% -------- first pass - initial estimation of position of eyeball rotation axes and pupil decentration

[ei, ctr, iris2pupil, pxl2mm, ~, fig, status] = calibrateRotAxesPupilDecenter(ei, cam, gv_cam, T2anc, Tuse, primary_anc, N, verbose);
if ~status, disp('Rotation axes calibration fail ...'); return; end

% -------- second pass - estimation of pupil elongation and pupillary circular axis angles

[ei, ctrn, status] = calibratePupilElongationPCA(ei, cam, gv_cam, T2anc, Tuse, ctr, pxl2mm, N, verbose);
if ~status, disp('Pupil calibration fail ...'); return; end

if verbose
    fprintf("Pupillary circular axis: [%5.2f, %5.2f]\n", ei.pca_angle(1), ei.pca_angle(2));
    fprintf("Pupillary refraction widening: %5.3f\n", ei.pupil_refr);
    fprintf("Pupil stretching: angle: %7.3f degree, strength: %6.4f\n", ei.pupil_angle, ei.pupil_stretch);
end

% -------- third pass - figure alpha angle and rotation axes based on ellipse ratio

[ei, alpha, gaze_deg, ir_ratio, el_ratio, el_pos, ctr_pos] = calibrateAlphaAngle(ei, cam, anc, T2anc, Tuse, ctr, pxl2mm, N);

if verbose
    fprintf('Angle alpha calibration complete: %5.2f %5.2f\n', alpha(1), alpha(2));
    
    % Secondary scale of rotation axes distances
    % secondary correction might be required due to dectenter, refraction, pupil imperfections not fully accounted for in previous passes
    gazedeg0m = gaze_deg(Tuse,:)+alpha;
    gdtrue0m = anc(T2anc(Tuse),:);
    alpha_scl = sum(gazedeg0m .* gdtrue0m) ./ sum(gazedeg0m.*gazedeg0m);
    
    % Note: alpha_scl is a scaling of X/Y components of gaze vector, not angles
    % For low angles it is nearly identical (and for high angles we have no clue how scaling behaves anyway)
    % The end result: low angles should be used during calibration to compute scaling
    fprintf('Rotation axes calibration (from ellipse ratio/rotation, for a reference only): H axis: %5.2fmm  V axis: %5.2fmm\n', ei.iris2hrot_mm/alpha_scl(1), ei.iris2vrot_mm/alpha_scl(2));
end

% -------- statistics

if verbose
    r_gt = gv_cam(:,3);
    figure; plot(1:N, el_ratio, 1:N, ir_ratio, 1:N, r_gt);
    legend('pupil ellipse ratio', 'iris ellipse ratio', 'Ground truth');
    
    figure(fig); hold on;
    clr = {"r","g","b","c","m","y","k","#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"};
    if ~isempty(ctrn)
        for i=1:N-1, plot(ctrn(i,1), ctrn(i,2), '*', 'Color', clr{i}); end
    end
    for i=1:N
        plot(ctr_pos(i,1), ctr_pos(i,2), 'x', 'Color', clr{i});
        plot(el_pos(i,1), el_pos(i,2), 'o', 'Color', clr{i});
        plot([ctr_pos(i,1) el_pos(i,1)], [ctr_pos(i,2) el_pos(i,2)], '-', 'Color', clr{i});
    end
    hold off;
    
    gaze_deg_alpha = gaze_deg+alpha;
    
    gd_true = anc(T2anc(Tuse),:);
    pv = 1:sum(Tuse);
    figure; 
    subplot(2,1,1); plot(pv, gaze_deg(Tuse,1), pv, gaze_deg_alpha(Tuse,1), pv, gd_true(:,1)); grid on; grid minor;
    legend('X', 'alpha corrected', 'Ground truth');
    subplot(2,1,2); plot(pv, gaze_deg(Tuse,2), pv, gaze_deg_alpha(Tuse,2), pv, gd_true(:,2)); grid on; grid minor;
    legend('Y', 'alpha corrected', 'Ground truth');
end

% -------- export calibraiton structure

calibration_string = [ 'struct(' ...
    ' ''iris2hrot'','  num2str(ei.iris2hrot_mm, '%5.2f') ...
    ', ''iris2vrot'',' num2str(ei.iris2vrot_mm, '%5.2f') ...
    ', ''iris2pupil'',' num2str(iris2pupil, '%4.2f') ...
    ', ''pupil_decenter'',[' num2str(ei.pupil_decenter_mm(1), '%5.2f') ' ' num2str(ei.pupil_decenter_mm(2), '%5.2f') ']' ...
    ', ''eyeball_ctr'',[' num2str(ei.eyeball_ctr(1), '%4.0f') ' ' num2str(ei.eyeball_ctr(2), '%4.0f') ']' ...
    ', ''pca_angle'',[' num2str(ei.pca_angle(1), '%7.2f') ' ' num2str(ei.pca_angle(2), '%7.2f') ']' ...
    ', ''pupil_refr'',' num2str(ei.pupil_refr, '%5.3f') ...
    ', ''pupil_angle'',' num2str(ei.pupil_angle, '%8.2f') ...
    ', ''pupil_stretch'',' num2str(ei.pupil_stretch, '%5.3f') ...
    ', ''gaze_alpha'',[' num2str(alpha(1), '%7.2f') ' ' num2str(alpha(2), '%7.2f') '] )' ];

fprintf("\nCalibration result:\n");
fprintf([calibration_string '\n']);

ec = eval(calibration_string);
