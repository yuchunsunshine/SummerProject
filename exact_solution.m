x = linspace(0, 2, 100);
y = linspace(0, 1, 100);
[xq, yq] = meshgrid(x, y);

zq = xq .* exp(yq);

figure;
surf(xq, yq, zq);
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Surface Plot of the exact solution  xe^{y} by MATLAB',FontSize = 15);