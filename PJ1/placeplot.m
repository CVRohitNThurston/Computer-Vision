function [varargout] = placeplot(varargin)
% Sizes and positions a plot.
%   placeplot() centers the plot, 
%   placeplot(position) puts the current plot in one of the following positions:
%   1-upper left, 2-upper right, 3-lower left, 4-lower right, 
%   5- full screen, 6-left side, 7-right side, 8-top, 9-botom
%
%   placeplot(position,fighandle) references the figure passed in by
%   fighandle. 
%   dim = placplot(...) returns the coordinates used for the plot size.
%   [dim mon] = placeplot(...) returns the coordinates and monitor size
%   [dim mon scn] =  placeplot(...) returns the dimensions of the plot,
%   the monitor size and screensize.

% interpret inputs/ set defaults
narginchk(0,2); nargoutchk(0,3); % check correct number of inputs/outputs in call
if (nargin<2); fighandle = gcf; else fighandle = varargin{2}; end;
if (nargin<1); position = 0; else position = varargin{1}; end;

% set figure position/size
set(0,'Units','pixels') % make pixel the units
scn=get(0,'ScreenSize');
mon=get(0,'Monitor');

% p(1) = left; p(2)= bottom; p(3)=width; p(4)= height;
% set width/height to section of screen
switch (position)
    case 1 % upper left
        dim = ([1, mon(1,4)/2-scn(1,2), mon(1,3:4)./2]);
    case 2 % upper right
        dim = ([mon(1,3)/2, mon(1,4)/2-scn(1,2), mon(1,3:4)./2]);
    case 3 % lower left
        dim = ([1, 1-scn(1,2), mon(1,3:4)./2]);
    case 4 % lower right
        dim = ([mon(1,3)/2, 1-scn(1,2), mon(1,3:4)./2]);
    case 5 % full screen
        dim = ([1, 1-scn(1,2), mon(1,3) mon(1,4)]);
    case 6 % left side
        dim = ([1, 1-scn(1,2), mon(1,3)/2, mon(1,4)]);
    case 7 % right side
        dim = ([mon(1,3)/2, 1-scn(1,2), mon(1,3)/2, mon(1,4)]);
    case 8 % top
        dim = ([1, mon(1,4)/2-scn(1,2), mon(1,3), mon(1,4)./2]);
    case 9 % bottom
        dim = ([1, 1-scn(1,2), mon(1,3), mon(1,4)./2]);
    otherwise % simply move to center
        movegui(fighandle,'center'); dim = get(fighandle,'OuterPosition');
end
set(fighandle,'OuterPosition',dim);

varargout{1} = dim;
varargout{2} = mon;
varargout{3} = scn;
end