% Solves nonlinear 1D Euler equations with moving mesh using FVM
clc
clear
close all

global gam wu wdinit meshType

% mesh parameters
dx = 0.05;
x = [0, dx/2:dx:1.0+dx/2]; % initial locations of cell faces
Nx = length(x)-1; % total number of cells
dxArray = x(2:Nx+1) - x(1:Nx); % initial cell widths

% set mesh velocity type
    % 0: only move 1st element face at x=0
    % 1: move 1st element face and linearly move 2nd element face
    % 2: move 1st element face and linearly move all other cells and faces
meshType = 0;

% Run time stepping yes or no
    % 1: yes
    % 0: no
RK4switch = 1;

% upstream conditions
Ms = 3;
gam = 1.4;
pu = 1; % non-dimensionalized in P and RT
RTu = 1;
cu = gam^(1/2);
uu = Ms*cu;
rhou = gam*pu/(cu^2);
rhoeu = pu/(gam-1);
Eu = rhoeu/rhou+(uu^2)/2;

wu = [rhou
      rhou*uu
	  rhou*Eu];

% downstream initial conditions
vs = 0.0;
Md = (uu-vs)/cu; % relative to shock
pd = pu*(2*gam*Md^2-(gam-1))/(gam+1);
rhod = rhou*(((gam+1.0)*Md^2)/((gam-1.0)*Md^2+2.0));
ud = vs+(rhou/rhod)*(uu-vs);
RTd = pd/rhod;
cd = (gam*RTd)^(1/2);
rhoed = pd/(gam-1);
Ed = rhoed/rhod+(ud^2)/2;

wdinit = [rhod
          rhod*ud
          rhod*Ed];
wd = repmat(wdinit,1,Nx); % array to store downstream cons. vars



% perturbations
if RK4switch == 1
    
    eps = 0.01;
    for i=1:size(wd,2)
        wd(:,i) = wd(:,i)*(1+eps*sin(i/size(wd,2)*pi));
    end
    
else
    
    eps = 0.001;
    Nvar = 4;
    M2 = zeros(4*Nx+1,4*Nx+1); % sized with 4 variables to compare to Fourier Analysis
    dwdt4var = zeros(4,size(wd,2));
    for i=1:size(wd,2) % loop over elements in reverse (to match Fourier Analysis)
        for var = 1:4 % loop over conserved variables
            if var ~= 2
                % adds perturbation skipping second var
                if var == 1
                    wd(var,i) = wd(var,i) + eps;
                else
                    wd(var-1,i) = wd(var-1,i) + eps;
                end
                [dwdt,dxdt] = get_dwdt_dxdt(wd,x);
                % calculate numerical derivative
                wnum = dwdt/eps;
                if i == 1
                    xnum = dxdt(1)/eps;
                end
                % removes perturbation skipping second var
                if var == 1
                    wd(var,i) = wd(var,i) - eps;
                else
                    wd(var-1,i) = wd(var-1,i) - eps;
                end
                wnum4var = [wnum(1,:);zeros(1,size(wnum,2));wnum(2:end,:)]; % add row of zeros for 4th var
            elseif var == 2
                % skip a column for rho*u
                wnum4var = zeros(4,size(wd,2));
                xnum = 0;
            end
            wnum4var = reshape(wnum4var,[],1); % reshape matrix into column array
            M2(:,Nvar*(i-1)+var) = [wnum4var; xnum];
            
        end
    end
    % perturb x-pos after loop over variables
    x(1) = x(1) + eps;
    [dwdt,dxdt] = get_dwdt_dxdt(wd,x);
    xnum = dxdt(1)/eps;
    x(1) = x(1) - eps;
    M2(end,end) = xnum;
    
    % rearrange matrix to match Fourier Analysis
    Inds = [];
    for n=Nx:-1:1
        Inds = [Inds n*4-3 n*4-2 n*4-1 n*4];
    end
    Inds = [Inds 4*Nx+1];


    M2 = M2(Inds,Inds);
    
end




% Set time parameters
t = 0;
tfinal = 2*x(end)/ud;
CFL = 0.5;
dt = CFL*min(dxArray)/ud; % use smallest dx for time step

% run RK4 loop
if RK4switch == 1
    while (t < tfinal)

        % RK4 time step
        [wd, x] = RK4(wd,x,dt);

        % increment time
        t = t + dt;

        % plot solution
        figure(1)
        tiledlayout(2,2);
        nexttile
        stairs(x,cat(2,wd(1,:),wd(1,Nx)))
        title('rho')
        nexttile
        stairs(x,cat(2,wd(2,:),wd(2,Nx)))
        title('rho*u')
        nexttile
        stairs(x,cat(2,wd(3,:),wd(3,Nx)))
        title('rho*E')

        figure(2)
        hold on
        plot(t,x(1),'ro')
        title('shock position')
        xlabel('t')
        ylabel('x')

        drawnow
    end
end


function [wd, x] = RK4(wd,x,dt)
[K1w,K1x] = get_dwdt_dxdt(wd,x);
[K2w,K2x] = get_dwdt_dxdt(wd + dt*K1w/2,x + dt*K1x/2);
[K3w,K3x] = get_dwdt_dxdt(wd + dt*K2w/2,x + dt*K2x/2);
[K4w,K4x] = get_dwdt_dxdt(wd + dt*K3w,x + dt*K3x);

wd = wd + 1/6*dt*(K1w + 2*K2w + 2*K3w + K4w);
x = x + 1/6*dt*(K1x + 2*K2x + 2*K3x + K4x);
end


function [dwdt,dxdt] = get_dwdt_dxdt(wd,x)
Nx = size(wd,2);
dwdt = zeros(3,Nx);
dxdt = get_dxdt(wd);
F = flux(wd,dxdt);
for i = 1:Nx
    dwdt(:,i) = 1/(x(i+1)-x(i)) * ( -wd(:,i)*(dxdt(i+1)-dxdt(i)) + F(:,i)-F(:,i+1) );
end
end


function [F] = flux2(w)
global gam
F = [w(2);
     w(2)^2/w(1)+(gam-1)*(w(3)-0.5*w(2)^2/w(1));
     w(2)/w(1)*(w(3)+(gam-1)*(w(3)-0.5*w(2)^2/w(1)))];
end


function [Fi] = flux(wd,dxdt)
global gam wu wdinit
wdBC = cat(2, wd, wdinit); % add ghost cell on right side with value of wdinit
Nx = size(wd,2);
Fi = zeros(3,Nx+1);

% %%%%%%%%%% symbolic version is slow %%%%%%%%%%
% syms w1 w2 w3 
% F = [w2;
%      w2^2/w1+(gam-1)*(w3-0.5*w2^2/w1);
%      w2/w1*(w3+(gam-1)*(w3-0.5*w2^2/w1))];
% J = jacobian(F,[w1,w2,w3]);
% [V,E] = eig(J);
% absA = V*abs(E)*inv(V);
% myA = matlabFunction(absA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fu = flux2(wu);

for i = 1:Nx+1
    if i == 1
        Fi(:,i) = Fu - dxdt(i)*wu;
    else
        % symbolic
        % absAn = myA(wdBC(1,i),wdBC(2,i),wdBC(3,i)); 
        
        %%%%%%%%% analytic version is fast %%%%%%%%%
        cd = sqrt(gam*(gam-1)/wdBC(1,i) * (wdBC(3,i) - ((wdBC(2,i))^2)/(2*wdBC(1,i))));
        ud = wdBC(2,i)/wdBC(1,i);
        E = [ud,    0,       0;
             0,     ud+cd,   0;
             0,     0,       ud-cd];
        V(:,1) = [2/ud^2; 2/ud; 1];
        V(:,2) = [(2*(gam - 1))/(gam*ud^2 - 2*cd*ud + 2*cd^2 - ud^2 + 2*cd*gam*ud);
                  (2*(cd + ud)*(gam - 1))/(gam*ud^2 - 2*cd*ud + 2*cd^2 - ud^2 + 2*cd*gam*ud); 1];
        V(:,3) = [(2*(gam - 1))/(2*cd*ud + gam*ud^2 + 2*cd^2 - ud^2 - 2*cd*gam*ud);
                  -(2*(cd - ud)*(gam - 1))/(2*cd*ud + gam*ud^2 + 2*cd^2 - ud^2 - 2*cd*gam*ud); 1];
        absAn = V*abs(E)*inv(V);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Fi(:,i) = 0.5*( flux2(wdBC(:,i))+flux2(wdBC(:,i-1)) ) ...
                    - 0.5*absAn*( wdBC(:,i)-wdBC(:,i-1) ) - dxdt(i)*0.5*( wdBC(:,i)+wdBC(:,i-1) );
    end
    end
end


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

