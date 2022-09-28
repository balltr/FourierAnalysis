close all
clear
clc


tantheta = [2.5, 5, 7.5, 10];

imag_eig_upwind = [2.3838, 5.4452, 8.1973, 10.9442];
FEM_upwind = [2.3111, 5.5851, 8.1733, 11.0387];

imag_eig_central = [2.7898, 5.4687, 8.2171, 10.9624];
FEM_central = [2.9139, 5.5851, 8.3256, 11.0387];

figure
plot(tantheta, imag_eig_upwind, 'ko-')
hold on
plot(tantheta, FEM_upwind, 'rx-')
hold off
xlabel('$tan(\theta)$','Interpreter','latex')
% ylabel('$Im(\lambda)$','Interpreter','latex')
title('Upwind Stabilization On')
legend('$Im(\lambda)$','FEM Oscillation','location','northwest','Interpreter','latex')
axis([2 10.5 2 12])
grid on

figure
plot(tantheta, imag_eig_central, 'ko-')
hold on
plot(tantheta, FEM_central, 'rx-')
hold off
xlabel('$tan(\theta)$','Interpreter','latex')
title('Upwind Stabilization Off')
legend('$Im(\lambda)$','FEM Oscillation','location','northwest','Interpreter','latex')
axis([2 10.5 2 12])
grid on