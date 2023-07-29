% Surface 1
data = importdata('data.txt');
x = data(:, 1);
y = data(:, 2);
w = data(:, 3);

[xq1, yq1] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
wq1 = griddata(x, y, w, xq1, yq1, 'cubic');

% Surface 2
x2 = linspace(0, 2, 100);
y2 = linspace(0, 1, 100);
[xq2, yq2] = meshgrid(x2, y2);
zq2 = xq2 .* exp(yq2);

% Plotting
figure;
subplot(1, 2, 1);
surf(xq1, yq1, wq1);
xlabel('x');
ylabel('y');
zlabel('w');
title('Surface 1');

subplot(1, 2, 2);
surf(xq2, yq2, zq2);
xlabel('x');
ylabel('y');
zlabel('z');
title('Surface 2');