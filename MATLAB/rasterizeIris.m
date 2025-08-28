% -------------------------------------------------------------------------
%
%    Rasterize iris sectors in pseudo-polar coordinates
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function irim = rasterizeIris(iris_d, iris_scl, im, cxy, min_r, max_r, iris_z, phi_minor, camera_rot)

    [h,w] = size(im);

    n = length(iris_d);
    k = max_r-min_r+1;

    irim = zeros(n,k,'uint8');

    for i=1:n
        dxy0 = iris_scl(i) .* [cos(iris_d(i)) sin(iris_d(i))];

        % correcting for human iris being 11.5x10.5 mm
        % The cornea is transparent and it's size is about 11.5mm h x 10.5mm v (= 0.913 ratio)
        % D. Atchison. 2023. Optics of the Human Eye: 2nd Edition, page 3
        dxys = rotcam( rotcam(dxy0, camera_rot) .* [1 0.95], -camera_rot);

        x = cxy(1) + (min_r:max_r) .* dxys(1);
        y = cxy(2) + (min_r:max_r) .* dxys(2);

        % camera projection correction
        % phi_minor always points away from camera (towards eyeball center)
        % range of iris_d: [-pi/6 .. pi-pi/6]
        % range of phi_minor: [-pi/2 .. 3*pi/2]
        % phi_minor origin is at 12 o'clock; iris_d origin is at 3 o'clock
        % +pi/2 to match origins; +2*pi to prevent negative values in mod
        ang_diff = mod(iris_d(i)-phi_minor+pi/2+2*pi, 2*pi);
        neg = ang_diff < pi/2 | ang_diff > 2*pi-pi/2;
        % negative height - this point on iris is closer to the camera than the imaging plane
        height = (min_r:max_r) * sqrt(1-iris_scl(i).^2) .* (1-neg*2);
        % no /2 here in h,w since im is 2x subsampled
        dx = height .* (x-1-w) ./ (height+iris_z);  % -1 = Matlab 1-based
        dy = height .* (y-1-h) ./ (height+iris_z);
        x = x - dx;
        y = y - dy;

        % divide coordinates by 2 since im is 2x subsampled
        ind = sub2ind(size(im),clip(round((y-1)/2)+1,1,h), clip(round((x-1)/2)+1,1,w));  % -1 +1 = Matlab 1-based
        irim(i,:) = im(ind);
    end

end
