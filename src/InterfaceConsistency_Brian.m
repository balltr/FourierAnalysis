% clc
clear
close all

dx = 1;

% upstream
gamma = 1.4;
pu = 1;
RTu = 1;
cu = sqrt(gamma*RTu);
Mu = 3;
vu = -Mu*cu;
rhou = gamma*pu/(cu^2);
tantheta = 0.1;
uu = -vu*tantheta;
rhoeu = pu/(gamma-1)
rhohu = gamma*rhoeu
rhoHu = rhohu +rhou*(uu^+0*vu^2)/2
Hu = rhoHu/rhou
Eu = rhoeu/rhou+(uu^2+vu^2)/2

wu = [rhou
      rhou*uu
	  rhou*vu
	  rhou*Eu];

% Stationary shock relations
Md = sqrt((1+(gamma-1)/2*Mu^2)/(gamma*Mu^2-(gamma-1)/2))
pd = pu*(1+gamma*Mu^2)/(1+gamma*Md^2)
RTd = RTu*(1+(gamma-1)/2*Mu^2)/(1+(gamma-1)/2*Md^2)
rhod = pd/RTd
cd = sqrt(gamma*RTd)
vd = -Md*cd
ud = uu
rhoed = pd/(gamma-1)
% Checking to make sure total enthalpy is conserved
rhohd = gamma*rhoed
rhoHd = rhohd +rhod*(ud^+0*vd^2)/2
Hd = rhoHd/rhod
Ed = rhoed/rhod+(ud^2+vd^2)/2

wd = [rhod
      rhod*ud
	  rhod*vd
	  rhod*Ed];
dydx = 0.0;

% Consistency Test (should be zero)
[dydt] = fun(gamma,wu,wd,dydx)

eps = 0.00001;
for test = 1:4
    wd(test) = wd(test)+eps;
    [dydt] = fun(gamma,wu,wd,dydx);
    coeff(test) = dydt/eps;
    wd(test) = wd(test)-eps;
end
dydx = eps;
[dydt] = fun(gamma,wu,wd,dydx);
coeff(5) = dydt/eps
double(coeff)

syms M
fM = 1/(gamma-1) * sqrt( (2*gamma*M^2 - (gamma-1)) * ((gamma-1)*M^2 + 2) / M^2) + (M^2-1)/M;
df = diff(fM,M);
M = Mu;
dfdMu_inv = 1/double(subs(df));

wcoeff(1) = dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)^2)*( (wd(2)^2+wd(3)^2)/wd(1) - wd(4)) - wd(3)/wd(1)^2 );
wcoeff(2) = dfdMu_inv*(gamma+1)/2 * ( -gamma*wd(2)/(cd*wd(1)^2) );
wcoeff(3) = dfdMu_inv*(gamma+1)/2 * ( 1/wd(1) - gamma*wd(3)/(cd*wd(1)^2) );
wcoeff(4) = dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)) );
wcoeff(5) = dfdMu_inv*(gamma+1)/2*(wu(2)/wu(1) - wd(2)/(wd(1))) - wu(2)/wu(1);

wcoeff





function [dydt] = fun(gamma,wu,wd,dydx)

norm = [dydx,-1]/sqrt(1+dydx^2)
qu = dot([wu(2)/wu(1),wu(3)/wu(1)],norm);
qd = dot([wd(2)/wd(1),wd(3)/wd(1)],norm);
cd = sqrt(gamma*(gamma-1)/wd(1) * (wd(4) - ((wd(2))^2 + (wd(3))^2)/(2*wd(1))));
cu = sqrt(gamma*(gamma-1)/wu(1) * (wu(4) - ((wu(2))^2 + (wu(3))^2)/(2*wu(1))));


syms Mu
f = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*cu) * (qu - qd + (2*cd)/(gamma-1));
df = diff(f,Mu);


% First guess is stationary shock
Mu = abs(wu(3)/(wu(1)*cu));
tol = 1e-5;
nmax = 50;
for n=1:nmax
    fval = subs(f);
    df = subs(df);
    Muold = Mu;
    Mu = Mu -fval/df;
    if (abs(Mu-Muold) < tol)
        break;
    end
end
if n >= nmax-1
   fprintf('nmax reached \n')
end

% eqn 36
x_tdotn = qu - cu*Mu;

% eqn 39 dot 38
dydt = x_tdotn/norm(2);

end


