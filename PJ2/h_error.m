function [err inliers] = h_error(h, x1, y1, x2, y2, thresh)

% computes error as the sum of square of norm
%
% h represents the transformation matrix [xp yp 1].' = h*[x y 1].'
%
% Inputs:
%   x,y -   reference image coordinates
%   xp,yp - transformed image coordinates
%   h -     transformation matrix
%
% assumes h is in the form h = 3x3 transformation matrix

% get transformed points, (xt,yt)
for i = length(x1):-1:1
    p = h*[x1(i) y1(i) 1].';
    xt(i) = p(1)./p(3);
    yt(i) = p(2)./p(3);
end

% get difference between transformed points and detected points
x_diff = xt - x2;
y_diff = yt - y2;

% get norm of difference
err = sqrt(x_diff.^2 + y_diff.^2);

% threshold for inliers/outliers
inliers = err < median(err);

% get total error
err = sum(err.^2);



