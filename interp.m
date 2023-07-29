% Read the data from the text file into MATLAB using the importdata function.
data = importdata('data.txt');
x = data(:, 1);
y = data(:, 2);
w = data(:, 3);

% Create a grid of points to define the domain of interpolation. 
% You can use the meshgrid function to generate a grid of x and y values. 
[xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% This will create a grid of 100x100 points within the range of x and y values.

wq = griddata(x, y, w, xq, yq, 'cubic');
% This will interpolate the w values at the points defined by xq and yq using cubic interpolation.

surf(xq, yq, wq);
xlabel('x');
ylabel('y');
zlabel('w');
title('Interpolation of the approximate value by MATLBA', FontSize = 15);