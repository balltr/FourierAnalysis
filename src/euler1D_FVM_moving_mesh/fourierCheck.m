% Create a matrix using numerical derivatives to compare with Fourier Analysis code
clc
clearvars

% mesh parameters
dx = 1; % dx for uniform cells
L = 4;
x = [0, dx/2:dx:L+dx/2]; % initial locations of cell faces
Nx = length(x)-1; % total number of cells

% upstream conditions
Mu = 3;
gam = 1.4;
pu = 1; % non-dimensionalized in P and RT
RTu = 1;
cu = gam^(1/2);
uu = Mu*cu;
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


M2 = M2(Inds,Inds)

