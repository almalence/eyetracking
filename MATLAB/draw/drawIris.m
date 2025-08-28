% -------------------------------------------------------------------------
%
%    Plot colored iris ellipse sectors
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function drawIris(ir_r, iris_d, cxy, el_ratio, phi_minor, iris_z, img_ctr, used, verbose, cneg, cpos)

    if verbose==0, return; end

    iris_d = iris_d(1:4:end);
    iris_scl = ellipseScales(iris_d, phi_minor, el_ratio);

    if ir_r>0
        x = cxy(1)+ir_r*iris_scl.*cos(iris_d);
        y = cxy(2)+ir_r*iris_scl.*sin(iris_d);

        % camera projection correction
        ang_diff = mod(iris_d-phi_minor+pi/2+2*pi, 2*pi);
        neg = ang_diff < pi/2 | ang_diff > 2*pi-pi/2;
        height = ir_r * sqrt(1-iris_scl.^2) .* (1-neg*2);
        dx = height .* (x-img_ctr(1)) ./ (height+iris_z);
        dy = height .* (y-img_ctr(2)) ./ (height+iris_z);
        x = x - dx;
        y = y - dy;

        figure(1);
        subplot(1,2,2);
        hold on;
        style = ['.', '-', '.'];
        plot (x(neg), y(neg), [style(used) cneg]);
        plot (x(~neg), y(~neg), [style(used) cpos]);
        hold off;
    end

end
