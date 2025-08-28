% -------------------------------------------------------------------------
%
%    Compensate pupil ellipse for pupillary circular axis inclination
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function [el_r, phi_minor] = compensatePCA(vec_to_eyeb, iris_ratio, el_r, phi_minor, pca_angle, max_el_ratio)

    % pupil 'gaze' vector from el_r, phi_minor
    xy = sqrt(1-el_r^2);
    pgv = [sin(phi_minor)*xy; -cos(phi_minor)*xy; el_r];

    % rotate gaze vector with pupillary circular axis angles
    pca_mat = rotationMat([pca_angle(2) pca_angle(1) 0], 0);

    % two directions possible due to ambiguity of pupil projected on image plane
    pgvr(:,1) = pca_mat * (pgv  .* [ 1;  1;  1]);
    pgvr(:,2) = pca_mat * (pgv  .* [-1; -1;  1]);

    % margin above which phi_minor direction may become confused
    % (minor ellipse axis confusion with major axis)
    el_thr = 0.96; % 0.96 = 16.2 degree
    if iris_ratio > el_thr || el_r > el_thr
        % pgv2 = possibility of el_minor/el_major swap
        xy2 = sqrt(1-max_el_ratio^2);
        pgv2 = [sin(phi_minor+pi/2)*xy2; -cos(phi_minor+pi/2)*xy2; max_el_ratio];
        pgvr(:,3) = pca_mat * (pgv2 .* [ 1;  1;  1]);
        pgvr(:,4) = pca_mat * (pgv2 .* [-1; -1;  1]);
    end
    % normalize vector
    xy_ref = sqrt(1-iris_ratio^2);
    vec_to_eyeb = [vec_to_eyeb'/norm(vec_to_eyeb)*xy_ref; iris_ratio];

    % choose phi_minor direction (out of 4) that best points towards eyeball center
    [~, el_best] = max( dot(repmat(vec_to_eyeb,1,size(pgvr,2)), pgvr) );

    % el_r, phi_minor from pupil 'gaze' vector
    el_r = pgvr(3,el_best);
    phi_minor = atan2(pgvr(2,el_best), pgvr(1,el_best))+pi/2;
end

