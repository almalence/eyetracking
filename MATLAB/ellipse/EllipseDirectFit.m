%
% License notice for this function:
%
% Copyright (c) 2009, Nikolai Chernov
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Alabama at Birmingham nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function A = EllipseDirectFit(XY)
%
%  Direct ellipse fit, proposed in article
%    A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%     "Direct Least Squares Fitting of Ellipses"
%     IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%  Our code is based on a numerically stable version
%  of this fit published by R. Halir and J. Flusser
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%  This is a fast non-iterative ellipse fit.
%
%  It returns ellipses only, even if points are
%  better approximated by a hyperbola.
%  It is somewhat biased toward smaller ellipses.
%
    A = [0 0 0 0 0 0]';
    
    centroid = mean(XY);   % the centroid of the data set
    
    D1 = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
          (XY(:,2)-centroid(2)).^2];
    D2 = [XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
    S1 = D1'*D1;
    S2 = D1'*D2;
    S3 = D2'*D2;
    T = -inv(S3)*S2';
    M = S1 + S2*T;
    M = [M(3,:)./2; -M(2,:); M(1,:)./2];
    [evec,eval] = eig(M);
    cond = 4*evec(1,:).*evec(3,:)-evec(2,:).^2;
    poscond = find(cond>0);
    if ~isempty(poscond)
        A1 = evec(:,poscond);
        A = [A1; T*A1];
        A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
        A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
        A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
             A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
        A(4) = A4;  A(5) = A5;  A(6) = A6;
        A = A/norm(A);
    end

end  %  EllipseDirectFit


