function h = gcurrmap(m)
%GCURRMAP   dark grey-white-Red colormap 
%   prepare for network_sims (in which graph default background is gray)
%   its designed to not draw edge currents when this are low
%   For example, to reset the colormap of the current figure:
%
%             colormap(gcurrmap)
%
%   See also HSV, PARULA, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

%   A.Diaz 01/11/2018

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

n = round(1/3*m);
n2=m-n;

% dark grey=[0.35 0.35 0.35];
% white=[1 1 1];
% red=[1 0 0];
r = [linspace(0.35,1,n)';ones(1,n2)'];
g = [linspace(0.35,1,n)';linspace(1,0,n2)'];
b = [linspace(0.35,1,n)';linspace(1,0,n2)'];

h = [r g b];