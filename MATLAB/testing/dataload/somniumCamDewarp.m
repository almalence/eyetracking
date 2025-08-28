% -------------------------------------------------------------------------
%
%    Somnium VR1 ET camera de-warping
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% An input frame for single eye is 640x400
% An ad-hoc center is at [360,200]
% For such center selection, camera is:
% - at 22 degree horizontal angle relative to HMD
% - at roughly 82mm eye relief distance
% For dewarped output frame of 540x400, the HFOV is about 30 degree (rounded from 28.5 calculated from experiment)
% eye=0 - left eye, eye=1 - right eye
function im = somniumCamDewarp(im, eye)

    % center of the output frame in input frame
    xc = 360; yc = 200;

    % width / height of the dewarped frame
    wd = 540; hd = 400;

    [~,w] = size(im);

    % distance from the center
    [y,x] = ndgrid((1:hd)-hd/2, (1:wd)-wd/2);
    r = sqrt(x.^2 + y.^2);

    % pincussion and perspective (due to off-axis camera placement behind the lens) corrections
    xd = x+xc - 0.16*x + 0.0012*x.^2 + 0.0007*y.^2 + 0.00055*x.*r;
    yd = y+yc - 0.10*y + 0.0007*x.*y +               0.00055*y.*r;

    if eye==0, xd = w-(xd-1); end

    im = uint8(interp2(single(im), xd, yd, 'linear', 255));

    % correcting light fall-off on the side farther away from camera
    im(:,241:wd) = uint8(single(im(:,241:wd)) .* (1+(0:300-1)*0.005));

    if eye==0, im = im(:,end:-1:1); end
end
