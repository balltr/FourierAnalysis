clc
clear
% close all


%%%%%%%%%%% 1D Wave Equation %%%%%%%%%%%%%%%
figure

a = 2.0;
dx = 2.0;
adis = 1.0;

for kdx=-pi:pi/100:pi
  ddx = 1i*sin(kdx)/dx;
  d2dx = -(2-2*cos(kdx))/dx^2;
  
  % du/dt = -(a du/dx -adis*|a|/2 d^2u/dx^2)
  rhs = -(a*ddx-adis*dx*abs(a)/2*d2dx);
%   rhs = -a/dx*( 1 - (cos(kdx)+1i*sin(kdx)) );
  
  E = eig(rhs);
  
  plot(real(E),imag(E),'x')
  hold on
end
title('1D Wave Eigenvalues')
grid on
xlabel('Re')
ylabel('Im')
axis equal








%% %%%%%%%%%%%%%%%%% 2D Wave Equation %%%%%%%%%%%%%%%%%%%
clc
clear
% close all

figure

ax = 1.0;
ay = 0.2;

dx = 1.0;
dy = 1.0;

y = 0:dy:30;
Ny = length(y);

adis = 1.0;

for kdx = 0-pi:pi/20:pi
  ddx = 1i*sin(kdx)/dx;
  d2dx = -(2-2*cos(kdx))/dx^2;
  
  % du/dt = -(ax du/dx -adis*|ax|/2 d^2u/dx^2 -(ay du/dy -adis*|ay|/2 d^2u/dy^2 )
  % (u_{j+1} - u_{j-1})/(2*dy) + adis|ay|/2 (u_{j+1}-2u_{j}+u_{j-1})/dy^2
  rhsp1 = -ay/(2*dy) + (adis*dy*abs(ay)/2)/dy^2;
  rhs = -(ax*ddx-adis*dx*abs(ax)/2*d2dx) - 2*(adis*dy*abs(ay)/2)/dy^2;
  rhsm1 = ay/(2*dy) + (adis*dy*abs(ay)/2)/dy^2;
  Mat = diag(rhs*ones(1,Ny)) + diag(rhsp1*ones(1,Ny-1),1) + diag(rhsm1*ones(1,Ny-1),-1);
  
  % Periodic BC's
  Mat(1,end) = rhsm1;
  Mat(end,1) = rhsp1;
  
  [V,E_mat] = eig(Mat);
  E_vec = diag(E_mat);
  
  plot(real(E_vec),imag(E_vec),'x')
  hold on
end
title('2D Wave Eigenvalues')
grid on
xlabel('Re')
ylabel('Im')
axis equal

% sort eigenvalues in descending order
[~, ind] = sort(real(E_vec),'descend');
E_mat = E_mat(ind, ind);
V = V(:, ind);

% for col = Ny:-1:1
% figure
% plot(real(V(1:end,col)))
% hold on
% plot(imag(V(1:end,col)))
% title('eigenvalue = ', num2str(E_mat(col,col)))
% end








%% %%%%%%%%%%%%%% 1D Euler Equations %%%%%%%%%%%%%%%%%%%%
clc
clear
close all

figure

dx = 2.0;
adis = 1.0;

syms w1 w2 w3 gamma
F = [w2;
     w2^2/w1+(gamma-1)*(w3-1/2*w2^2/w1);
     w2/w1*(w3+(gamma-1)*(w3-1/2*w2^2/w1))];
J = jacobian(F,[w1,w2,w3]);
syms rho u c
A = subs(J,[w1,w2,w3],[rho,rho*u,rho*(c^2/(gamma*(gamma-1))+1/2*u^2)]);
A = simplify(A);
[V,E] = eig(A);

gamma = 1.4;
pu = 1.015e5; % Pa
Tu = 300; % K
Rbar = 8314; % j/(kmol K)
MW = 28.9;
R = Rbar/MW;
c = sqrt(gamma*R*Tu);
rho = pu/(R*Tu);
eu = R*Tu/(gamma-1);
u = 0.1*c;

A = double(subs(A));
V = double(subs(V));
E = double(subs(E));

absA = V*abs(E)*inv(V); %#ok<*MINV>

for kdx=-pi:pi/100:pi
    
  ddx = 1i*sin(kdx)/dx;
  d2dx = -(2-2*cos(kdx))/dx^2;

  % du/dt = -(a du/dx -adis*|a|/2 d^2u/dx^2)
  rhs = -(A*ddx-adis*dx*absA/2*d2dx);
  
  E = eig(rhs);
  
  plot(real(E),imag(E),'x')
  hold on
end
title('1D Euler Eigenvalues')
grid on
xlabel('Re')
ylabel('Im')
axis equal








%% %%%%%%%%%%%%%% 2D Euler Equations %%%%%%%%%%%%%%%%%%%%%%
% Conservative variable w = rho, rho u, rho v, rho E
clc
clear
close all

figure

syms w1 w2 w3 w4;
syms gamma;
e = [w2;
     w2^2/w1+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w2*w3/w1;
     w2/w1*(w4+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1))];
f = [w3;
     w2*w3/w1;
     w3^2/w1+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w3/w1*(w4+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1))];
 
numVar = length(e);

Ae = jacobian(e,[w1,w2,w3,w4]);
Af = jacobian(f,[w1,w2,w3,w4]);
syms rho u v gamma c
Ae = subs(Ae,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gamma*(gamma-1)) +1/2*(u^2+v^2))});
Af = subs(Af,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gamma*(gamma-1)) +1/2*(u^2+v^2))});

gamma = 1.4;
tantheta = 1;

% flow variables
p = 1;
RT = 1;
c = sqrt(gamma*RT);
M = 3;
v = -M*c; % velocity normal to shock
u = -v*tantheta;
rho = gamma*p/(c^2);

Ae = double(subs(Ae));
[Ve,Ee] = eig(Ae);
Ve = double(subs(Ve));
Ee = double(subs(Ee));
Af = double(subs(Af));
[Vf,Ef] = eig(Af);
Vf = double(subs(Vf));
Ef = double(subs(Ef));

absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);

dx = 1.0;
dy = 1.0;
y = 0:dy:59;
Ny = length(y);
adis = 1.0;

for kdx = pi/2
% for kdx = -pi:pi/5:pi

  ddx = 1i*sin(kdx)/dx;
  d2dx = -(2-2*cos(kdx))/dx^2;
  
  rhsp1 = -Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
  rhs = -(Ae*ddx-adis*dx*absAe/2*d2dx) - 2*(adis*dy*absAf/2)/dy^2;
  rhsm1 = Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
  
  Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
  
  % Periodic BC's
  Mat(1:numVar,end-numVar+1:end) = rhsm1;
  Mat(end-numVar+1:end,1:numVar) = rhsp1;
  
  E = eig(Mat);
  
%   plot(real(E),imag(E),'x')
  scatter(real(E),imag(E),'k','filled')
  hold on
end
% title('Fully Periodic 2D Euler Eigenvalues')
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
axis equal