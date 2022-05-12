function [dxdt] = get_dxdt(wdFull)
global gam wu meshType
Nx = size(wdFull,2);
dxdt = zeros(1,Nx+1);
wd = wdFull(:,1);
norm = 1;
qu = dot(wu(2)/wu(1),norm);
qd = dot(wd(2)/wd(1),norm);
cu = sqrt(gam*(gam-1)/wu(1) * (wu(3) - wu(2)^2/(2*wu(1))) );
cd = sqrt(gam*(gam-1)/wd(1) * (wd(3) - wd(2)^2/(2*wd(1))) );

syms Mu
f_sym = 1/(gam-1) * sqrt( (2*gam*Mu^2 - (gam-1)) * ((gam-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gam+1)/(2*cu) * (qu - qd + (2*cd)/(gam-1));
df_sym = diff(f_sym,Mu);
myf = matlabFunction(f_sym);
mydf = matlabFunction(df_sym);

% First guess is stationary shock
Mu = abs(wu(2)/(wu(1)*cu));
tol = 1e-5;
nmax = 50;
for n=1:nmax
    fval = myf(Mu);
    df = mydf(Mu);
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
dxdt(1) = x_tdotn/norm;

% calculate dxdt for all other faces
if meshType == 0
    % do nothing
elseif meshType == 1
    dxdt(2) = 0.5*dxdt(1);
elseif meshType == 2
    for i = 2:Nx+1
        dxdt(i) = dxdt(i-1)/2;
    end
else
    fprintf('Unrecognized meshType. Value of 0 selected.\n')
end

end