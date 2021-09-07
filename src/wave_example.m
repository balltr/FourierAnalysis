clc
clear
close all


%%%%%%%%%%% 1D Wave Equation %%%%%%%%%%%%%%%
figure

a = 1.0;
dx = 1.0;
adis = 1.0;

for k=-pi:pi/100:pi
  ddx = 1i*sin(k*dx)/dx;
  d2dx = -(2-2*cos(k*dx))/dx^2;
  
  % du/dt = -(a du/dx -adis*|a|/2 d^2u/dx^2)
  rhs = -(a*ddx-adis*abs(a)/2*d2dx);
  
  [V,E] = eig(rhs);
  
  plot(real(E),imag(E),'x')
  hold on
end
axis equal





%% %%%%%%%%%%%%%%%%% 2D Wave Equation %%%%%%%%%%%%%%%%%%%
figure

ax = 1.0;
ay = 0.2;

dx = 1.0;
dy = 1.0;

y = 0:dy:30;
Ny = length(y);

adis = 1.0;

for k = -pi:pi/20:pi
  ddx = 1i*sin(k*dx)/dx;
  d2dx = -(2-2*cos(k*dx))/dx^2;
  
  % du/dt = -(ax du/dx -adis*|ax|/2 d^2u/dx^2 -(ay du/dy -adis*|ay|/2 d^2u/dy^2 )
  % (u_{j+1} - u_{j-1})/(2*dy) + adis|ay|/2 (u_{j+1}-2u_{j}+u_{j-1})/dy^2
  rhsp1 = -ay/(2*dy) + (adis*abs(ay)/2)/dy^2;
  rhs = -(ax*ddx-adis*abs(ax)/2*d2dx) - 2*(adis*abs(ay)/2)/dy^2;
  rhsm1 = ay/(2*dy) + (adis*abs(ay)/2)/dy^2;
  Mat = diag(rhs*ones(1,Ny)) + diag(rhsp1*ones(1,Ny-1),1) + diag(rhsm1*ones(1,Ny-1),-1);
  
  % Periodic BC's
  Mat(1,end) = rhsm1;
  Mat(end,1) = rhsp1;
  
  [V,E] = eig(Mat);
  
  plot(real(E),imag(E),'x')
  hold on
end
axis equal





%% %%%%%%%%%%%%%% 1D Euler Equations %%%%%%%%%%%%%%%%%%%%
figure

dx = 1.0;
adis = 1.0;

syms w1 w2 w3 ws gam real
F = [w2;
     w2^2/w1+(gam-1)*(w3-1/2*w2^2/w1);
     w2/w1*(w3+(gam-1)*(w3-1/2*w2^2/w1))];
J = jacobian(F,[w1,w2,w3]);
syms rho u c
A = subs(J,[w1,w2,w3],[rho,rho*u,rho*(c^2/(gam*(gam-1))+1/2*u^2)]);
A = simplify(A);
[V,E] = eig(A);

gam = 1.4;
P0 = 1.015e5; % Pa
T0 = 300; % K
Rbar = 8314; % j/(kmol K)
MW = 28.9;
R = Rbar/MW;
c = sqrt(gam*R*T0);
rho = P0/(R*T0);
e0 = R*T0/(gam-1);
u = 0.1*c;

A = double(subs(A));
V = double(subs(V));
E = double(subs(E));

absA = V*abs(E)*inv(V);


for k=-pi:pi/100:pi
  ddx = 1i*sin(k*dx)/dx;
  d2dx = -(2-2*cos(k*dx))/dx^2;
  
  % du/dt = -(a du/dx -adis*|a|/2 d^2u/dx^2)
  rhs = -(A*ddx-adis*absA/2*d2dx);
  
  [V,E] = eig(rhs);
  
  plot(real(E),imag(E),'x')
  hold on
end
title('1D Euler Eigenvalues')
xlabel('Re')
ylabel('Im')
axis equal






%% %%%%%%%%%%%%%% 2D Euler Equations %%%%%%%%%%%%%%%%%%%%%%
% Conservative variable w = rho, rho u, rho v, rho E
figure

syms w1 w2 w3 w4 positive;
syms gam positive;
e = [w2;
     w2^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w2*w3/w1;
     w2/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
f = [w3;
     w2*w3/w1;
     w3^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w3/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
 
numVar = length(e);

Ae = jacobian(e,[w1,w2,w3,w4]);
Af = jacobian(f,[w1,w2,w3,w4]);
syms rho u v gam c positive
Ae = subs(Ae,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
Af = subs(Af,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});

[Ve,Ee] = eig(Ae);
[Vf,Ef] = eig(Af);

gam = 1.4;
P0 = 1.015e5; % Pa
T0 = 300; % K
Rbar = 8314; % j/(kmol K)
MW = 28.9;
R = Rbar/MW;
c = sqrt(gam*R*T0);
rho = P0/(R*T0);
e0 = R*T0/(gam-1);
u = 4*c;
v = -0.7*c;

Ae = double(subs(Ae));
Ve = double(subs(Ve));
Ee = double(subs(Ee));
Af = double(subs(Af));
Vf = double(subs(Vf));
Ef = double(subs(Ef));

absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);

dx = 1.0;
dy = .5;
y = 0:dy:10;
Ny = length(y);
adis = 1.0;

% for k = pi/30
for k = -pi:pi/30:pi

  ddx = 1i*sin(k*dx)/dx;
  d2dx = -(2-2*cos(k*dx))/dx^2;
  
  rhsp1 = -Af/(2*dy) + (adis*absAf/2)/dy^2;
  rhs = -(Ae*ddx-adis*absAe/2*d2dx) - 2*(adis*absAf/2)/dy^2;
  rhsm1 = Af/(2*dy) + (adis*absAf/2)/dy^2;
  
  Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
  
  % Periodic BC's
  Mat(1:numVar,end-numVar+1:end) = rhsm1;
  Mat(end-numVar+1:end,1:numVar) = rhsp1;
  
  [V,E] = eig(Mat);
  
  plot(real(E),imag(E),'x')
  hold on
end
title('2D Euler Eigenvalues')
xlabel('Re')
ylabel('Im')
axis equal






%% %%%%%%%%%%%%%% 2D Euler Equations With Shock On Top Boundary %%%%%%%%%%%%%%%%%%%%%%
% Conservative variable w = rho, rho u, rho v, rho E
clear
figure

syms w1 w2 w3 w4 positive;
syms gam positive;
e = [w2;
     w2^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w2*w3/w1;
     w2/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
f = [w3;
     w2*w3/w1;
     w3^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w3/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
w = [w1; w2; w3; w4];
 
numVar = length(e);

Ae = jacobian(e,[w1,w2,w3,w4]);
Af = jacobian(f,[w1,w2,w3,w4]);
syms rho1 u v gam c1 positive
Ae = subs(Ae,{w1,w2,w3,w4},{rho1,rho1*u,rho1*v,rho1*(c1^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
Af = subs(Af,{w1,w2,w3,w4},{rho1,rho1*u,rho1*v,rho1*(c1^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
w = subs(w,{w1,w2,w3,w4},{rho1,rho1*u,rho1*v,rho1*(c1^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
eEval = subs(e,{w1,w2,w3,w4},{rho1,rho1*u,rho1*v,rho1*(c1^2/(gam*(gam-1)) +1/2*(u^2+v^2))});

[Ve,Ee] = eig(Ae);
[Vf,Ef] = eig(Af);

gam = 1.4;
MW = 28.9;
Rbar = 8314; % j/(kmol K)
R = Rbar/MW;

% upstream properties
Ms = 2;
P0 = 1.197e3; % Pa
T0 = 250; % K
rho0 = P0/(R*T0);
e0 = R*T0/(gam-1);
c0 = sqrt(gam*R*T0);

% downstream properties
P1 = P0*(2*gam*Ms^2-(gam-1))/(gam+1);
T1 = T0*((2*gam*Ms^2 - (gam-1)) * ((gam-1)*Ms^2+2) / ((gam+1)^2*Ms^2));
rho1 = P1/(R*T1);
e1 = R*T1/(gam-1);
c1 = sqrt(gam*R*T1);

% downstream velocities
u = c1 * 4;
v = c1 * -0.7;

Ae = double(subs(Ae));
Ve = double(subs(Ve));
Ee = double(subs(Ee));
Af = double(subs(Af));
Vf = double(subs(Vf));
Ef = double(subs(Ef));
w = double(subs(w));
eEval = double(subs(eEval));

absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);

dx = 1.0;
dy = 0.5;
y = 0:dy:10;
Ny = length(y);
adis = 1.0;

% 0 for central diff or 1 for upwind
fdSwitch = 1;
% 0 for full eig spectrum, 1 for eig pairs
eigSwitch = 0;
% 0 for shock off, 1 for shock on
shockSwitch = 0;

if eigSwitch == 0
    maxWaveNum = pi;
    minWaveNum = -pi;
    incWaveNum = pi/20;
    waveNumArray = minWaveNum:incWaveNum:maxWaveNum;
elseif eigSwitch == 1
    eigPair = pi/10; % +\-
    waveNumArray = [-eigPair, eigPair];
end

CM = turbo(length(waveNumArray));
colorIdx = 0;

for ii = 1:length(waveNumArray)
    k = waveNumArray(ii);
    colorIdx = colorIdx + 1;
    
    ddx = 1i*sin(k*dx)/dx;
    d2dx = -(2-2*cos(k*dx))/dx^2;
    
    if shockSwitch == 0
        rhsp1 = -Af/(2*dy) + (adis*absAf/2)/dy^2;
        rhs = -(Ae*ddx-adis*absAe/2*d2dx) - 2*(adis*absAf/2)/dy^2;
        rhsm1 = Af/(2*dy) + (adis*absAf/2)/dy^2;
  
        Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
  
        % Periodic BC's
        Mat(1:numVar,end-numVar+1:end) = rhsm1;
        Mat(end-numVar+1:end,1:numVar) = rhsp1;
    elseif shockSwitch == 1

        % interior
        rhsp1 = -Af/(2*dy) + (adis*absAf/2)/dy^2;
        rhs = -(Ae*ddx-adis*absAe/2*d2dx) - 2*(adis*absAf/2)/dy^2;
        rhsm1 = Af/(2*dy) + (adis*absAf/2)/dy^2;
        Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
        Mat(end-numVar+1:end,:) = 0; % zeros bottom entries for shock boundary

        % Moving shock BC on top
        wRHSM1 = 1/dy*(Af+absAf);
        wRHS = -2*(Ae*ddx-adis*absAe*dx/2*d2dx) - 1/dy*(Af+absAf);
        %   rhsp1 = zeros(numVar);

        Mat(end-numVar+1:end,end-numVar+1:end) = wRHS;
        Mat(end-numVar+1:end,end-2*numVar+1:end-numVar) = wRHSM1;
        %   Mat(end-numVar+1:end,numVar+1:2*numVar) = wRHSP1;

        if fdSwitch == 0
            % central difference
            ddx_y = 1i*sin(k*dx)/dx;
%             d2dx = -(2-2*cos(k*dx))/dx^2;
        elseif fdSwitch == 1
            % upwind (backward difference)
            ddx_y = (1i*sin(k*dx) - cos(k*dx) + 1)/dx;
%             d2dx = (1i*(2*sin(k*dx)-sin(2*k*dx)) - 2*cos(k*dx)+cos(2*k*dx) + 1)/dx^2;
        end

        % for y'
        yRHS = -4*eEval/dy*ddx;
        Mat((Ny-1)*numVar+1:Ny*numVar, Ny*numVar+1) = yRHS;
        % this comes from wave equation of y'
        Mat(Ny*numVar+1,Ny*numVar+1) = -u/dx*ddx_y;
    
    end

    [V,E] = eig(Mat);

    plot(real(E),imag(E),'color',CM(colorIdx,:),'marker','x','LineStyle','none')
    hold on
end
if shockSwitch == 0
    title('2D Euler With Shock Eigenvalues')
elseif shockSwitch == 1
    title('2D Euler Eigenvalues')
end
xlabel('Re')
ylabel('Im')
axis equal
grid on
colormap('turbo')
colorbar('Ticks',[])