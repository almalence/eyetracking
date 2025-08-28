% -------------------------------------------------------------------------
%
%    Image brightening towards corners
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

function imbr = cornerBrighten(im, valid)

    [h,w] = size(im);
    [y,x] = ndgrid((0:int32(h)-1)-valid(2), (0:int32(w)-1)-valid(1));
    r2 = valid(3)*valid(3)+valid(4)*valid(4);
    % matching to C code
    brc = idivide(valid(6)*65536, r2);    % brightening coeff
    imbr = uint8( int32(im) + idivide(idivide(int32(im).*(x.*x+y.*y), 65536) * brc + 128, 256) );

end
