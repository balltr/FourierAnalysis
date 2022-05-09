clc
clear
close all

kdx = 0:pi/20:pi;

% 1D heat equation with central diff scheme
h_num = -2*(1-cos(kdx));
h_an = -(kdx).^2;
figure
plot(kdx,h_num,'linewidth',1)
hold on
grid on
plot(kdx,h_an,'linewidth',1)
legend('Numerical','Analytical')
xlabel('k \Delta x','fontsize',14), ylabel('\lambda\Delta x^2 / \nu','fontsize',14)
title('1D Heat Diffusion Equation With Central Difference Scheme','fontsize', 14)
hold off

% 1D wave equation with central diff scheme
w1_num = sin(kdx);
w1_an = kdx;
figure
plot(kdx,w1_num,'linewidth',1)
hold on
grid on
plot(kdx,w1_an,'linewidth',1)
legend('Numerical','Analytical','location','northwest')
xlabel('k \Delta x','fontsize',14), ylabel('Im(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation With Central Difference Scheme','fontsize', 14)
hold off

% 1D wave equation with upwind scheme
w2_num_Im = sin(kdx);
w2_an_Im = kdx;
figure
subplot(1,2,1)
plot(kdx,w2_num_Im,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Im,'linewidth',1)
legend('Numerical','Analytical','location','northwest')
xlabel('k \Delta x','fontsize',14), ylabel('Im(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation With Upwind Scheme','fontsize', 14)
hold off

w2_num_Re = -(1 - cos(kdx));
w2_an_Re = 0*kdx;
subplot(1,2,2)
plot(kdx,w2_num_Re,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Re,'linewidth',1)
legend('Numerical','Analytical','location','northeast')
xlabel('k \Delta x','fontsize',14), ylabel('Re(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation With Upwind Scheme','fontsize', 14)
