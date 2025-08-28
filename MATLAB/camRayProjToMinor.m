% -------------------------------------------------------------------------
%
%    Compute angle between ellipse minor axis and projected camera ray
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function ang = camRayProjToMinor(xyc, xye, iris_z)

    % angle between center of pupil/iris ellipse towards camera center and towards eyeball center
    alpha = atan2(xyc(2), xyc(1)) - atan2(xye(2), xye(1));

    % angle of camera ray hitting imaging plane at (x,y)
    % cra = atan2(norm(xyc), iris_z);

    % projected onto plane coinciding with minor ellipse axis and perpendicular to the imaging plane
    % ang = atan(tan(cra) * cos(alpha));
    ang = atan(norm(xyc) / iris_z * cos(alpha));
end
