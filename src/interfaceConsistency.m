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
Eu = pu/(rhou*(gamma-1))+0.5*uu*vu

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
wdprime = [0.000001
           0.0000
           0.000000
           0.0000];
dypdxi = 0.000;
barwd = wd;
wd = wd+wdprime;

% eqn 38
norm = [1/dx*dypdxi,-1];

% eqn 34 and eqn 57
[Mu,dydt1,dydt57,qu] = fun(gamma,cu,dx,norm,wu,wd,dypdxi,wdprime,barwd);
Mu
dydt1
dydt57
dydt1/dydt57

% eqn 36
x_tdotn = qu - cu*Mu;

% eqn 39 dot 38
dydt2 = -x_tdotn

coeffRatio = dydt2 / dydt1





function [Mu,dydt,dydt2,qu] = fun(gamma,cu,dx,norm,wu,wd,dypdxi,wdprime,barwd)
cd = sqrt( gamma*(gamma-1)/wd(1) * (wd(4) - ((wd(2))^2 + (wd(3))^2)/(2*wd(1))) );
barcd = sqrt( gamma*(gamma-1)/barwd(1) * (barwd(4) - ((barwd(2))^2 + (barwd(3))^2)/(2*barwd(1))) );
qu = dot([wu(2)/wu(1),-wu(3)/wu(1)],norm);
qd = dot([wd(2)/wd(1),-wd(3)/wd(1)],norm);

syms Mu
f = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*cu) * (qu - qd + (2*cd)/(gamma-1));
df = diff(f,Mu);
dfdMu_inv = 1/df;
tol = 1e-5;
Mu = 5;
nmax = 50;
if double(subs(df)) ~= 0
   n = 1;
   while n < nmax
       n = n+1;
       M_u1 = Mu - double(subs(f))/double(subs(df));
       if abs(M_u1-Mu) < tol
           Mu = M_u1;
           break
       end
       Mu = M_u1;
       if n == nmax
           fprintf('nmax reached \n')
       end
   end
end
dfdMu_inv = double(subs(dfdMu_inv));

% wcoeff1 = wdprime(1) * dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)^2)*( (wd(2)^2+wd(3)^2)/wd(1) - wd(4)) - wd(3)/wd(1)^2 );
% wcoeff2 = wdprime(2) * dfdMu_inv*(gamma+1)/2 * ( -gamma*wd(2)/(cd*wd(1)^2) );
% wcoeff3 = wdprime(3) * dfdMu_inv*(gamma+1)/2 * ( 1/wd(1) - gamma*wd(3)/(cd*wd(1)^2) );
% wcoeff4 = wdprime(4) * dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)) );
% dydxicoeff = dypdxi * 1/dx * ( dfdMu_inv*(gamma+1)/2*(wu(2)/wu(1) - wd(2)/(wd(1))) - wu(2)/wu(1) );

wcoeff1 = wdprime(1) * dfdMu_inv*(gamma+1)/2 * ( gamma/(barcd*barwd(1)^2)*( (barwd(2)^2+barwd(3)^2)/barwd(1) - barwd(4)) - barwd(3)/barwd(1)^2 );
wcoeff2 = wdprime(2) * dfdMu_inv*(gamma+1)/2 * ( -gamma*barwd(2)/(barcd*barwd(1)^2) );
wcoeff3 = wdprime(3) * dfdMu_inv*(gamma+1)/2 * ( 1/barwd(1) - gamma*barwd(3)/(barcd*barwd(1)^2) );
wcoeff4 = wdprime(4) * dfdMu_inv*(gamma+1)/2 * ( gamma/(barcd*barwd(1)) );
dydxicoeff = dypdxi * 1/dx * ( dfdMu_inv*(gamma+1)/2*(wu(2)/wu(1) - barwd(2)/(barwd(1))) - wu(2)/wu(1) );

dydt = wcoeff1 + wcoeff2 + wcoeff3 + wcoeff4 + dydxicoeff;






% test eqn 57 directly
dqu = wu(2)/wu(1) * 1/dx*dypdxi;
dqd = barwd(2)/barwd(1) * 1/dx*dypdxi + barwd(3)*wdprime(1)/barwd(1)^2 - wdprime(3)/barwd(1);
w1coeffdcd = wdprime(1) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1)^2 * ( (barwd(2)^2+barwd(3)^2)/barwd(1) - barwd(4) ) );
w2coeffdcd = wdprime(2) * gamma*(gamma-1)/(2*barcd) * ( -barwd(2)/barwd(1)^2 );
w3coeffdcd = wdprime(3) * gamma*(gamma-1)/(2*barcd) * ( -barwd(3)/barwd(1)^2 );
w4coeffdcd = wdprime(4) * gamma*(gamma-1)/(2*barcd) * ( 1/barwd(1) );
dcd = w1coeffdcd + w2coeffdcd + w3coeffdcd + w4coeffdcd;

dydt2 = (dqu - dqd + 2/(gamma-1)*dcd) * (gamma+1)/2 * dfdMu_inv - 1/dx*dypdxi*wu(2)/wu(1);
    
end

