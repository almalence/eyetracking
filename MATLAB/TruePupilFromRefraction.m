% -------------------------------------------------------------------------
%
%    Refraction-induced correction of pupil ellipse parameters
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% Values are derived from Zemax model assuming 7.75mm cornea curvature,
% pupil/iris 3.3mm behind, 4mm pupil diameter, refraction index of 1.3375
% Accounting for non-flat iris in real eyes, the offset values are lower
%
% Eq 3.5 formula for observed pupil ellipse ratio is given in
% D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 52:
% 0.99*cos((angle+5.3)/1.121), for this table - excluding 5.3 term: cos((0:5:90)*pi/180/1.121)'
%
% angle	ctr offset	--- Pupil ellipse ratio ---
%       Zemax       Zemax   cos(angle)  Eq 3.5a
%  0	0.000	    1.000   1.000       1.0000
%  5	0.048	    0.997   0.996       0.9970
% 10	0.094	    0.988   0.985       0.9879
% 15	0.143	    0.972   0.966       0.9729
% 20	0.193	    0.951   0.940       0.9519
% 25	0.245	    0.924   0.906       0.9252
% 30	0.299	    0.892   0.866       0.8929
% 35	0.357	    0.855   0.819       0.8552
% 40	0.419	    0.813   0.766       0.8123
% 45	0.484       0.768   0.707       0.7644  
% 50	0.554	    0.719   0.643       0.7120
% 55	0.63	    0.667   0.574       0.6552
% 60	0.711       0.612   0.500       0.5945
% 65	0.798	    0.556   0.423       0.5302
% 70	0.891	    0.497   0.342       0.4626
% 75	0.992	    0.437   0.259       0.3923
% 80	1.097	    0.375   0.174       0.3195
% 85                        0.087       0.2449                          
% 90                        0           0.1687
%
% true cos(angle) ratio from observed pupil ratio:             cos( acos(pu_r) * 1.121 )
% pupil-iris mm shift on imaging plane from true pupil ratio:
%  - due to pupil closer to eyeball center:                    -(iris2eye-pupil2eye) * sin( acos(pu_r) )
%  - due to refraction:                                        acos(pu_r) * 0.32 + ( acos(pu_r) * 0.58 )^2

% Note: the output center (cxy) is in the iris ellipse plane, not pupil plane (pupillary circular axis and iris2pupil corrections applied)
function [cxy, el_major, el_minor, phi_minor] = ...
    TruePupilFromRefraction(cxy, semx, semy, phi, pxl2mm, eyeball2d, pupil_refr, pca_angle, iris_ratio, iris_z, img_ctr, iris2pupil, max_el_ratio)

    el_major = max(semx, semy);
    el_r = min(semx, semy) / el_major;

    % phi_minor = minor axis angle from vertical
    phi_minor = phi + (semx<semy)*pi/2;

    % ----- Correct pupil ellipse widening due to refraction and camera projection

    % camera ray projected onto plane perpendicular to imaging plane and coinciding with minor ellipse axis
    cra_proj = camRayProjToMinor(img_ctr-cxy, eyeball2d-cxy, iris_z);

    % perspective- and refraction- corrected ellipse ratio
    el_r = cos(acos(el_r*cos(cra_proj))*pupil_refr - cra_proj);

    % correction for pupillary circular axis via 3d re-projection
    [el_r, phi_minor] = compensatePCA(eyeball2d-cxy, iris_ratio, el_r, phi_minor, pca_angle, max_el_ratio);
    % phi_minor points toward 2d eyeball center, range: [-pi/2 .. 3*pi/2]

    el_minor = el_major * el_r;

    % ----- Center offset due to refraction

    % using stabilized ellipse ratio predicted from previous frame (large jitter in gaze angles otherwise)
    mm_shift = -iris2pupil*sin(acos(iris_ratio)) + acos(iris_ratio)*0.32 + (acos(iris_ratio)*0.58)^2;
    d =  [sin(phi_minor) -cos(phi_minor)] * mm_shift/pxl2mm;
    cxy  = cxy+d;

end
