function [pts, err] = dotField2GKA(radii, fieldRadius, showPlot)
% Generate a field of non-overlapping dots of various radii
% Usage:
%     pts = dotField2(radii, fieldRadius)
%     pts = dotField2(radii, fieldRadius, showPlot)
% "radii" is a vector of N dot radii, in minimal discrete units (e.g.,
% pixels). "fieldRadius" is the radius of the whole field (in the same
% units). "showPlot", if true, shows a little plot of what the dots will
% look like; defaults to false.
%
% GKA August 2012

% Field coordinates
[xx, yy] = meshgrid(-fieldRadius:fieldRadius);

% Distance to nearest object
d = fieldRadius - sqrt(xx.^2 + yy.^2);

% Coordinates of each point
x = zeros(numel(radii), 1);
y = zeros(numel(radii), 1);

for i = 1:numel(radii)
    is_valid = d > radii(i);
    xo = xx(is_valid);
    if isempty(xo)
        warning('Ran out of space before placing all dots');
        err = 1;
        break;
    else
        err = 0;
    end
    yo = yy(is_valid);
    ind = randi(numel(xo));
    x(i) = xo(ind);
    y(i) = yo(ind);
    
    new_d = sqrt((x(i)-xx).^2 + (y(i)-yy).^2) - radii(i);
    d = min(new_d, d);
end

pts = [x,y];

if exist('showPlot', 'var') && showPlot
    clf
    for i=1:numel(radii)
        rectangle('Position', [pts(i,:)-radii(i), radii(i)*[2 2]], ...
            'Curvature', [1 1]);
    end
    rectangle('Position', fieldRadius * [-1 -1 2 2], 'Curvature', [1 1]);
    axis equal
end