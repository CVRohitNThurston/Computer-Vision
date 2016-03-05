function err = h_error(h, x, y, xp, yp)

% computes error as the sum of square of norm
%
% h represents the transformation matrix [xp yp 1].' = h*[x y 1].'
%
% Inputs:
%   x,y -   reference image coordinates
%   xp,yp - transformed image coordinates
%   h -     transformation matrix
%
% assumes h is in the form [h11 h12 h13 h21 h22 h23 h31 h32].'

% h = 3x3 transformation matrix
h = reshape([h; 1],3,3);

% get transformed points, (xt,yt)
for i = length(x):-1:1
    p = h*[x(i) y(i) 1].';
    xt(i) = p(1);
    yt(i) = p(2);
end

% get difference between transformed points and detected points
x_diff = xt - xp;
y_diff = yt - yp;

% get square of norm of difference
err = x_diff.^2 + y_diff.^2;

% get total error
err = sum(err.^2);



