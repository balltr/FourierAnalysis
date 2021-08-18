clc
clear
close all

% 1D wave equation with upwind flux
kdx = -pi:pi/60:pi;

M = [2/3, 1/3; 1/3, 2/3];
B = [-.5, -.5; .5, -.5];
A = [0, 1; 0, 0];
a = 1;
dx = 1/12;

for z = 1:length(kdx)
    eigval = eig(M\(B+A*exp(1i*kdx(z))) * (2*a/dx)) * dx/a;
    [tmp,ind] = sort(abs(eigval));
    eigval = eigval(ind);
    h_num(z,1) = eigval(1);
    h_num(z,2) = eigval(2);
end
h_num %= h_num/2


% imaginary plane plot
figure
plot(h_num(:,1),'b+')
hold on
plot(conj(h_num(:,2)),'r+')
hold off


figure
for i = 1:size(h_num,2)
    plot (kdx, imag(h_num(:,1)),'+', kdx, -imag(h_num(:,2)), '+')
end

%%
L1_num_Im = imag(h_num(:,1));
w2_an_Im = kdx;
figure
subplot(1,2,1)
plot(kdx,L1_num_Im,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Im,'linewidth',1)
legend('Numerical','Analytical','location','northwest')
xlabel('k \Delta x','fontsize',14), ylabel('Im(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation DG P=1','fontsize', 14)
hold off

L1_num_Re = real(h_num(:,1));
w2_an_Re = 0*kdx;
subplot(1,2,2)
plot(kdx,L1_num_Re,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Re,'linewidth',1)
legend('Numerical','Analytical','location','northeast')
xlabel('k \Delta x','fontsize',14), ylabel('Re(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation DG P=1','fontsize', 14)

L1_num_Im = imag(h_num(:,2));
w2_an_Im = kdx;
figure
subplot(1,2,1)
plot(kdx,L1_num_Im,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Im,'linewidth',1)
legend('Numerical','Analytical','location','northwest')
xlabel('k \Delta x','fontsize',14), ylabel('Im(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation DG P=1','fontsize', 14)
hold off

L1_num_Re = real(h_num(:,2));
w2_an_Re = 0*kdx;
subplot(1,2,2)
plot(kdx,L1_num_Re,'linewidth',1)
hold on
grid on
plot(kdx,w2_an_Re,'linewidth',1)
legend('Numerical','Analytical','location','northeast')
xlabel('k \Delta x','fontsize',14), ylabel('Re(\lambda \Delta x / a)','fontsize',14)
title('1D Wave Equation DG P=1','fontsize', 14)