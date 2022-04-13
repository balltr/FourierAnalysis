clc
clear
close all

syms gam M cu cd uu ud

f = 1/(gam-1) * sqrt( (2*gam*M^2 - (gam-1)) * ((gam-1)*M^2 + 2) / M^2) + (M^2-1)/M - ...
    (gam+1)/(2*cu) * (uu - ud + (2*cd)/(gam-1));
1/diff(f,M)
dMduu = diff(f,uu) * 1/diff(f,M);

dMdud = diff(f,ud)% * 1/diff(f,M);

dMdcd = diff(f,cd)% * 1/diff(f,M);
pretty(dMdcd)

gam = 1.4;
M = 3;
cu = 1;
cd = .2;
uu = M*cu;

Md = sqrt( ((gam-1)*M^2+2)/(2*gam*M^2 - (gam-1)) );
ud = Md*cd;

double(subs(dMduu))
% double(subs(dMdud))

%%
clc
clear all

syms wp1 wp2 wp3 wp4 dypdxi invdfdM Mu

M = 3;
gamma = 1.4;
c_u = gamma^(1/2);
dx = 1;
v_u = M*c_u;
u_u = 0.1*v_u;
w1 = 3.8571;
w2 = 1.3691;
w3 = 3.5496;
w4 = 27.7097;


q_u = u_u * (1/dx * dypdxi) - v_u;
q_d = w2/w1 * (1/dx * dypdxi) + w3/w1 + (wp3/w1 - w3*wp1/w1^2);

barc_d = sqrt( gamma*(gamma-1)/w1 * (w4 - (w2^2 + w3^2)/(2*w1)) );
a = gamma*(gamma-1)/w1 * ( wp4 - (w2*wp2 + w3*wp3)/w1 + wp1*(w2+w3)/(2*w1^2) ) - wp1/w1*barc_d^2;
c_d = barc_d + a/(2*barc_d);

f = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*c_u) * (q_u - q_d + (2*c_d)/(gamma-1));
diff(f,Mu)
invdfdM = 1/diff(f,Mu)

dq_u = u_u * (1/dx * dypdxi);
dq_d = w2/w1 * (1/dx * dypdxi) + (wp3/w1 - w3*wp1/w1^2);
dc_d = a/(2*barc_d);

dyp1dtau = -1/dx*dypdxi + c_u * invdfdM * ( -(gamma+1)/(2*c_u)*dq_u + (gamma+1)/(2*c_u)*dq_d - (gamma+1)/(c_u*(gamma-1))*dc_d );

% coeffs(dyp1dtau,dypdxi)
% coeffs(f,u_un);

%%
clc
clear
syms w1 w2 w3 w4 wp1 wp2 wp3 wp4 eps dx dypdxi gamma

c_dexact = sqrt( gamma*(gamma-1)/(w1+wp1) * (w4+wp4 - ((w2+wp2)^2 + (w3+wp3)^2)/(2*(w1+wp1))) );

barc_d = sqrt( gamma*(gamma-1)/w1 * (w4 - (w2^2 + w3^2)/(2*w1)) );
a = gamma*(gamma-1)/(2*barc_d) * (  1/w1^2*((w2^2+w3^2)/w1 - w4)*wp1 - w2/w1^2*wp2 - w3/w1^2*wp3 + 1/w1*wp4);
c_d = barc_d + a;

q_dexact = (w2+wp2)/(w1+wp1) * 1/dx * dypdxi - (w3+wp3)/(w1+wp1);
q_d = w2/w1 * 1/dx * dypdxi - w3/w1 - ( wp3/w1 - w3*wp1/w1^2 );

gamma = 1.4;
w1 = 3.857152857142857;
w2 = 1.369149892660197;
w3 = 3.549647869859768;
w4 = 27.709666666666671;
wp1 = 0.1;
wp2 = 0.0;
wp3 = 0.0;
wp4 = 0.0;
dx = 0.2;
dypdxi = 1;

double(subs(c_dexact)) / double(subs(c_d))

double(subs(q_dexact)) / double(subs(q_d))

%%
clc
clear
close all

syms w1 w2 w3 w4 wp1 wp2 wp3 wp4 eps dx dypdxi gamma u_u v_u Mu c_u

barc_d = sqrt( gamma*(gamma-1)/w1 * (w4 - (w2^2 + w3^2)/(2*w1)) );
a = gamma*(gamma-1)/(2*barc_d) * (  1/w1^2*((w2^2+w3^2)/w1 - w4)*wp1 - w2/w1^2*wp2 - w3/w1^2*wp3 + 1/w1*wp4);
c_d = barc_d + eps*a;

q_u = u_u * eps * (1/dx * dypdxi) - v_u;
q_d = w2/w1 * eps*(1/dx * dypdxi) - w3/w1 - eps*(wp3/w1 - w3*wp1/w1^2);

f = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*c_u) * (q_u - q_d + (2*c_d)/(gamma-1));
diff(f,Mu);
invdfdM = 1/diff(f,Mu);

dq_u = u_u * eps * (1/dx * dypdxi);
dq_d = w2/w1 * eps*(1/dx * dypdxi) - eps*(wp3/w1 - w3*wp1/w1^2);
dc_d = eps*a;

eqn57 = invdfdM * (gamma+1)/2 * (dq_u - dq_d + 2/(gamma-1)*dc_d) - eps*(1/dx * dypdxi)*u_u;


gamma = 1.4;
Mu = 3;
u_u = 3.55;
v_u = -2;
c_u = 2;
w1 = 10;
w2 = 15;
w3 = 21;
w4 = 40;
wp1 = 20;
wp2 = 20;
wp3 = 20;
wp4 = 20;
eps = 0.01;
dx = 0.2;
dypdxi = 0.5;

double(subs(eqn57))


%%%%%% coefficients %%%%%%

dypdxicoeff = dypdxi * 1/dx * ( invdfdM * (gamma+1)/2 * (eps*u_u - eps*w2/w1) - eps*u_u );
wp1coeff = wp1 * ( invdfdM * (gamma+1)/2 * ( eps*gamma/(c_d*w1^2)*( (w2^2+w3^2)/w1 - w4 ) - eps*w3/w1^2) );
wp2coeff = wp2 * ( invdfdM * (gamma+1)/2 * -eps*(gamma*w2)/(c_d*w1^2) );
wp3coeff = wp3 * ( invdfdM * (gamma+1)/2 * (eps/w1 - eps*gamma*w3/(c_d*w1^2)) );
wp4coeff = wp4 * ( invdfdM * (gamma+1)/2 * (eps*gamma/(c_d*w1)) );

eqnxx = dypdxicoeff + wp1coeff + wp2coeff + wp3coeff + wp4coeff;
double(subs(eqnxx))

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear
close all

syms M_u

gamma = 1.4;
tantheta = 0.1;

%%%%% upstream %%%%%
p_u = 1;
E_u = 10;
c_u = gamma^(1/2);
rho_u = gamma*p_u/(c_u^2);
v_u = 3.5496;
u_u = v_u*tantheta;

w_u = zeros(4,1);
w_u(1) = rho_u;
w_u(2) = rho_u * u_u;
w_u(3) = rho_u * v_u;
w_u(4) = rho_u * E_u;

%%%%% downstream %%%%%
v_s = 0.0;
M_d = (v_u-v_s)/c_u; % relative to shock
p_d = p_u*(1.0+((2.0*gamma)/(gamma+1.0))*(M_d^2-1.0));
E_d = 10;
rho_d = rho_u*(((gamma+1.0)*M_d^2)/((gamma-1.0)*M_d^2+2.0));
v_d = v_s+(rho_u/rho_d)*(v_u-v_s);
u_d = u_u;
RT_d = p_d/rho_d;
c_d = (gamma*RT_d)^(1/2);

w_d = zeros(4,1);
w_d(1) = rho_d;
w_d(2) = rho_d * u_d;
w_d(3) = rho_d * v_d;
w_d(4) = rho_d * E_d;
%%%%%%%%%%%%%%%%%%%%%%


n = [0,-1];
U_u = [u_u, v_u];
U_d = [u_d, v_d];

% c_d = sqrt( gamma*(gamma-1)/w_d(1) * (w_d(4) - (w_d(2)^2 + w_d(3)^2)/(2*w_d(1))) );


R = 2*c_d/(gamma-1) - dot(U_d,n);

f = 1/(gamma-1) * sqrt( (2*gamma*M_u^2 - (gamma-1)) * ((gamma-1)*M_u^2 + 2) / M_u^2) + (M_u^2-1)/M_u - ...
    (gamma+1)/(2*c_u) * (dot(U_d,n) - dot(U_d,n) + (2*c_d)/(gamma-1));
df = diff(f,M_u);
invdfdM = 1/diff(f,M_u);
tol = 1e-5;
M_u = 10;
nmax = 100;
if double(subs(df)) ~= 0
   n = 1;
   while n < nmax
       n = n+1;
       M_u1 = M_u - double(subs(f))/double(subs(df));
       if abs(M_u1-M_u) < tol
           M_u = M_u1;
           break
       end
       M_u = M_u1;
   end
end
M_u


us = v_u - M_u*c_u





%%%%%% coefficients %%%%%%
wp1 = 0.1;
wp2 = 0.0;
wp3 = 0.0;
wp4 = 0.0;
w_d(1) = w_d(1) + wp1;
w_d(2) = w_d(2) + wp2;
w_d(3) = w_d(3) + wp3;
w_d(4) = w_d(4) + wp4;

% dypdxicoeff = dypdxi * 1/dx * eps*( invdfdM * (gamma+1)/2 * (w2/w1 - u_u) - u_u );
wp1coeff = wp1 * eps*( invdfdM * (gamma+1)/2 * (w_d(3)/w_d(1)^2 - gamma*(w_d(2)+w_d(3))/(2*w_d(1)^3*c_d) - c_d/(w_d(1)*(gamma-1)) ) );
% wp2coeff = wp2 * eps*( invdfdM * (gamma+1)/2 * (gamma*w2)/(c_d*w1^2) );
% wp3coeff = wp3 * eps*( invdfdM * (gamma+1)/2 * (-1/w1 - gamma*w3/(c_d*w1^2)) );
% wp4coeff = wp4 * eps*( invdfdM * (gamma+1)/2 * (gamma/(c_d*w1)) );

double(subs(wp1coeff));

%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clc
clear

syms xxi yN ypN yNhalf w wp eps

eqn1 = xxi * ( yN + eps*ypN - yNhalf ) * ( w + eps*wp );

eqn2 = xxi*yN*w + eps*xxi*ypN*w - xxi*yNhalf*w + eps*xxi*yN*wp - eps*xxi*yNhalf*wp;

xxi = 10;
yN = 5;
ypN = 2;
yNhalf = 3;
w = 5;
wp = 2;
eps = 0.01;

double(subs(eqn1))
double(subs(eqn2))

%%

clc
clear
close all

dx = 1.0;

% upstream
gamma = 1.4;
pu = 1;
RTu = 1;
cu = gamma^(1/2);
vu = 3*cu;
rhou = gamma*pu/(cu^2);
tantheta = 0.1;
uu = vu*tantheta;
Eu = pu/(rhou*(gamma-1))+0.5*uu*vu;

wu = [rhou
      rhou*uu
	  rhou*vu
	  rhou*Eu];

% downstream
pd = 10.333333333333334;
rhod = 3.857142857142857;
vd = 0.920279077371051;
ud = uu;
RTd = 2.679012345679012;
Ed = 7.183987654320989;

wd = [rhod
      rhod*ud
	  rhod*vd
	  rhod*Ed];

% perturbations
wdprime = [0.0001
           0.0000
           0.0000
           0.0000];
dypdxi = 0.000;
wd = wd+wdprime;

syms Mu
f = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*cu) * (qu - qd + (2*cd)/(gamma-1));
df = diff(f,Mu);
dfdMu_inv = 1/df;


dqu = uu/dx*dypdxi;
dqd = barwd(2)/barwd(1) * 1/dx*dypdxi + barwd(3)*wdprime(1)/barwd(1)^2 - wdprime(3)/barwd(1);
barcd = sqrt( gamma*(gamma-1)/barwd(1) * (barwd(4) - ((barwd(2))^2 + (barwd(3))^2)/(2*(barwd(1)))) );
w1coeff = wdprime(1) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1)^2 * ( (barwd(2)^2+barwd(3)^2)/barwd(1) - barwd(4) ) );
w2coeff = wdprime(2) * gamma*(gamma-1)/(2*barcd) * ( -barwd(2)/barwd(1)^2 );
w3coeff = wdprime(3) * gamma*(gamma-1)/(2*barcd) * ( -barwd(3)/barwd(1)^2 );
w4coeff = wdprime(4) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1) );
dcd = w1coeff + w2coeff + w3coeff + w4coeff;

eqn57 = (dqu - dqd + 2/(gamma-1)*dcd) * (gamma+1)/2 * dfdMuinv - 1/dx*dypdxi*uu
