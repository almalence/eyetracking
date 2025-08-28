% -------------------------------------------------------------------------
%
%    Estimate eyeball center location in 3d
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function eyeball3d = EyeballCenter3d(iris_z, old_z, iris2d, eyeball2d, iris2ctr_pxl, pxl2mm)

    % Note: here we are working in parallel-projection mode (camera perspective projection is unrelated)

    % squared distance between iris and eyeball center projection on an image plane
    iris2eyeball2d_sq = sum((iris2d-eyeball2d).^2);

    % eyeball z distance
    if iris2eyeball2d_sq > iris2ctr_pxl^2
        % eyeball center too far may happen due to position smoothing in Estimate2dEyeCenter
        eyeball3d = [pxl2mm*eyeball2d old_z];
    else
        eyeball3d = pxl2mm * [eyeball2d iris_z+sqrt(iris2ctr_pxl^2-iris2eyeball2d_sq)];
    end

end
