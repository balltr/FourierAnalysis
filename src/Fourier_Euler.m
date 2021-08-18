close all
clear
clc

addpath(genpath('/Users/Tristan/Documents/Research/functions/'))

%%%%%%%%% FIND EIGENVALUES AND EIGENVECTORS %%%%%%%%%%%%%%%%%

syms w1 w2 w3 gam
F = [w2;
     w2^2/w1+(gam-1)*(w3-1/2*w2^2/w1);
     w2/w1*(w3+(gam-1)*(w3-1/2*w2^2/w1))];
J = jacobian(F,[w1,w2,w3]);

syms rho u c
A = subs(J,[w1,w2,w3],[rho,rho*u,rho*(c^2/(gam*(gam-1))+1/2*u^2)]);
A = simplify(A);
[V1,E1] = eig(A)

% In simplest form
E = [u, 0, 0; 0, u+c,  0; 0, 0, u-c]
V = [1, 1, 1; u, u+c, u-c; u^2/2, (u^2/2+c*u)+c^2/(gam-1), (u^2/2 -c*u)+c^2/(gam-1)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set FluxSwitch to 0 for characteristic upwinding, else for Lax-Friedrich's
FluxSwitch = 0;

gam = 1.4;
P0 = 1.015e5; % Pa
T0 = 300; % K
Rbar = 8314; % j/(kmol K)
MW = 28.9;
R = Rbar/MW;
c0 = sqrt(gam*R*T0);
rho0 = P0/(R*T0);
e0 = R*T0/(gam-1);
nVar = 3;

% Numerical parameters
order = 1;
dx = L/N;
nodes = lglnodes(order);    % determine Gauss-Lobatto nodal points
[intPts,wgts] = GaussQuad(order);   % determine Gauss integration points and weights
phi = evalBasis(nodes,intPts); % evaluate basis functions at nodes
massBlock = getMass(order+1,phi,intPts,wgts)


kdx = -pi:pi/60:pi;


if FluxSwitch == 0
    % Characteristic Upwinding
    
else
    % Lax-Friedrich's
    
end




