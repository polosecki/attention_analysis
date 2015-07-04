function [im, clim] = nantowhite(cvals, clim, cmap)
% nantowhite - Convert a 2D matrix into an image with NaNs replaced by white
%
% USAGE:
%   [im, clim] = nantowhite(cvals)
%   [im, clim] = nantowhite(cvals, clim)
%   [im, clim] = nantowhite(cvals, clim, cmap)
%
% nantowhite takes in an N x M matrix of values ('cvals') and outputs an
% N x M x 3 uint8 image that can be displayed on screen using 'image'. This
% function replaces 'imagesc' for situations where there are NaN values in
% the input matrix. Instead of displaying 'NaN' values using the first
% value in the colormap, it displays them using white. You can replace what
% would originally be 'imagesc(cvals)' with 'image(nantowhite(cvals))'.
%
% This function was created because alpha layers used to make 'NaN'
% values transparent were ignored when images were saved in a vector format.
% See the demo code below to get 'colorbar' to show the correct limits when
% using this function.
% 
% INPUT:
%   cvals - A 2D matrix to be converted to an image. Any NaN values in the
%           matrix will be converted into white pixels in the output image.
%
%   clim - The minimum and maximum value on the colorscale used for
%          converting values in 'cvals' into colors in the output image.
%          An empty maxtix is ignored. Positive or negative 'inf' values
%          are replaced by the minimum or maximum value in 'cvals'.
%
%   cmap - The color map to use for conversion from 'cvals' into colors.
%
% OUTPUT:
%   im   - An N x M x 3 uint8 image matrix. [N,M] = size(cvals)
%
%   clim - The color limits used when converting values to color pixels.
%
% DEMO:
%    [X,Y] = meshgrid(-2*pi:pi/16:2*pi);
%    cvals = 20*cos(X).*sin(Y);
%    cvals(1:20,1:20) = NaN;
%    cvals(randperm(numel(cvals),400)) = NaN;
%
%    Original Code:
%       figure; ax = axes();
%       imagesc(cvals,'Parent',ax);
%       colorbar('Peer',ax);
%
%   Modified Code:
%       figure; ax = axes();
%       [im, cl] = nantowhite(cvals);
%       image(im,'Parent',ax,'CDataMapping','scaled');
%       set(ax,'CLim',cl);
%       colorbar('Peer',ax);
%
% See also image, imagesc, colormap
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2013, Benjamin Kraus
% $Id: nantowhite.m 4862 2013-05-27 03:42:40Z bkraus $

% Get the default color map.
if(nargin < 3); cmap = get(0,'DefaultFigureColormap');
elseif(ischar(cmap)); cmap = feval(cmap,64);
end
ncolors = size(cmap,1);

% Determine the range of the color values provided.
finites = isfinite(cvals(:));
if(~any(finites))
    cvalmin = 0;
    cvalmax = 0;
else
    cvalmin = min(cvals(finites));
    cvalmax = max(cvals(finites));
end

% Process the CLim input, if one is provided.
% If clim isn't provided, default to [-inf inf], which will will be set
% later to match the range of cavls.
if(nargin < 2 || isempty(clim)); clim = [-inf inf];
    
% If clim is a single value, treat it as a maximum (which cannot be -inf).
elseif(isscalar(clim) && clim~=-inf); clim = [-inf clim];
elseif(isscalar(clim));
    error('nantowhite:InvalidCLim','CLim(2) cannot be -inf.');
    
% Otherwise, make sure it is a two value vector.
elseif(numel(clim)~=2);
    error('nantowhite:InvalidCLim','CLim must be a scalar or 2 element vector.');

% Special case, the clim(1) and clim(2) are equal to each other, and equal
% to all values in cvals.
elseif(clim(1)==clim(2) && clim(1)==cvalmin && cvalmin==cvalmax); clim = [-inf inf];
    
% Otherwise, clim(2) must be greater than clim(1).
elseif(any(isnan(clim)) || clim(2)<=clim(1));
    error('nantowhite:InvalidCLim','CLim values must be increasing and non-NaN.');
end

% If clim(1) is -inf, replace with actual minimum.
% Three possibilities:
%   clim(2) == inf, we want clim(1) = cvalmin
%   clim(2) > cvalmin, we want clim(1) = cvalmin
%   clim(2) <= cvalmin, we want clim(1) < clim(2)
if(isinf(clim(1))); clim(1) = min([clim(2)-eps(clim(2))*64, cvalmin]); end

% If clim(2) is inf, replace with actual maximum.
% Three possibilities:
%   clim(1) == -inf, set above to be cvalmin, we want clim(2) > clim(1)
%   clim(1) < cvalmax, we want clim(2) = cvalmax
%   clim(1) >= cvalmax, we want clim(2) > clim(1)
if(isinf(clim(2))); clim(2) = max([clim(1)+eps(clim(1))*64, cvalmax]); end

% Distribute the bins equally among the CLim.
cbins = linspace(clim(1), clim(2), ncolors+1);

% Make the first and last bin '-inf' and 'inf' to include all values.
cbins(1) = -inf;
cbins(end) = inf;

% Bin the color values into the appropriate color bin.
[~,cvals] = histc(cvals,cbins);

% Move 'inf' values to the last bin.
cvals(cvals==ncolors+1) = ncolors;

% Initialize an entirely white image.
im = ones([size(cvals) 3]);

% Copy the appropriate values into the red, green, and blue layers.
for ii = 1:3;
    imlayer = im(:,:,ii);
    imlayer(cvals~=0) = cmap(cvals(cvals~=0),ii);
    im(:,:,ii) = imlayer;
end

end
