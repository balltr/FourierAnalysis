% Solves nonlinear 1D Euler equations on a moving mesh using FVM and RK4
clc
clearvars

global gam wu wdinit meshType

% MESH PARAMETERS
dx = 1; % dx for uniform cells
L = 80;
x = [0, dx/2:dx:L+dx/2]; % initial locations of cell faces
Nx = length(x)-1; % total number of cells
dxArray = x(2:Nx+1) - x(1:Nx); % initial cell widths

% SET MESH VELOCITY TYPE
    % 0: only move 1st element face at x=0
    % 1: move 1st element face and linearly move 2nd element face
    % 2: move 1st element face and linearly move all other cells and faces
meshType = 0;

% EXPORT VARIABLE TYPE
    % 0: conservative
    % 1: primitive
varType = 1;

% UPSTREAM CONDITIONS
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

% DOWNSTREAM INITIAL CONDITIONS
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

% SET INITIAL CONDITION
[x,wd] = initialCond(x,wd);

% TIME PARAMETERS
t = 0;

% CFL = 0.5;
% dt = CFL*min(dxArray)/ud; % use smallest dx for time step
% tfinal = 2*x(end)/ud;

dt = 0.1;
ntstep = 200;

% MAKE RESULTS DIRECTORY
if isfolder('../euler1D_FVM_moving_mesh/Results')
    delete('../euler1D_FVM_moving_mesh/Results/*')
else
    mkdir '../euler1D_FVM_moving_mesh' Results
end

% EXPORT BASE FLOW STATE AND INITIAL CONDITION
baseflowExport(wdinit,varType)
dataExport(0,x,wd,varType);

% RUN RK4 LOOP
itr = 0;
% while (t < tfinal)
for i = 1:ntstep
    itr = itr + 1;

    % RK4 time step
    [wd, x] = RK4(wd,x,dt);

    % increment time
    t(itr+1) = t(itr) + dt;

    % export files
    dataExport(itr,x,wd,varType);
end
tExport(t)

