
clear
clc


v0 = 3.55;
eps = 0.01;

x0 = 0:0.01:1;
x1 = 0:-160;

v = v0*ones(length(x0),length(x1)) + eps*sin(2*pi*x0).*exp(x1);
figure
plot(x0,v(1,:))



