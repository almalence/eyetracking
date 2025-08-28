% -------------------------------------------------------------------------
%
%    Speckle removal filter
% 
%    SPDX-FileCopyrightText: Copyright (c) 2024-2025, Almalence Inc.
%    SPDX-License-Identifier: GPL-3.0-only
%    All rights reserved.
%
%    Author: Dmitry Shmunk
%
% -------------------------------------------------------------------------

% will not keep thin features, but will keep strong step-like changes in brightness
function imf = removeSpeckles2x(im)

    [h,w] = size(im);

    % d = [4 4];
    imfv = zeros(h,w, 'uint16');
    im = uint16(im);
    for y=1:2:h
        % up = int16(mean( im(max(1,y-d(2)):y-1,:), 1 ));
        % dn = int16(mean( im(y+1:min(h,y+d(2)),:), 1 ));
        % matching to optimized C version
        up = ( (im(max(1,y-4),:) + im(max(1,y-3),:))/2 + (im(max(1,y-2),:) + im(max(1,y-1),:))/2 ) / 2;
        dn = ( (im(min(h,y+4),:) + im(min(h,y+3),:))/2 + (im(min(h,y+2),:) + im(min(h,y+1),:))/2 ) / 2;

        imfv(y,:) = clip(im(y,:), min(up,dn), max(up,dn));
    end

    imf = zeros(h/2,w/2, 'uint8');
    for x=1:2:w
        % lt = int16(mean( imfv(1:2:end,max(1,x-d(1):x-1)), 2 ));
        % rt = int16(mean( imfv(1:2:end,x+1:min(w,x+d(1))), 2 ));
        lt = ( (imfv(1:2:end,max(1,x-4)) + imfv(1:2:end,max(1,x-3)))/2 + (imfv(1:2:end,max(1,x-2)) + imfv(1:2:end,max(1,x-1)))/2 ) / 2;
        rt = ( (imfv(1:2:end,min(w,x+4)) + imfv(1:2:end,min(w,x+3)))/2 + (imfv(1:2:end,min(w,x+2)) + imfv(1:2:end,min(w,x+1)))/2 ) / 2;

        imf(:,(x-1)/2+1) = clip(imfv(1:2:end,x), min(lt,rt), max(lt,rt));
    end
end
