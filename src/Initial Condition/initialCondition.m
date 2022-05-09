% Creates initial condition for Oblique_Shock_Periodic test case
% Adds perturbation from eigenvector to base flow state

close all
clearvars
clc

% change these to match desired input files
V_filename = 'V_pi_10_upwind.txt';
regular_b0_filename = 'pi_10_b0.grd';
regular_b1_filename = 'pi_10_b1.grd';

% IC upstream
gamma = 1.4;
pu = 1; % non-dimensionalized in P and RT
RTu = 1;
Mu = 3.0;
cu = sqrt(gamma*RTu);
vu = -Mu*cu; % velocity normal to shock
tantheta = 0.1;
uu = -vu*tantheta;
rhou = gamma*pu/(cu^2);
rhoeu = pu/(gamma-1);
Eu = rhoeu/rhou+(uu^2+vu^2)/2;

% IC downstream
vs = 0.0;
Md = (vu-vs)/cu; % relative to shock
pd = pu*(2*gamma*Md^2-(gamma-1))/(gamma+1);
rhod = rhou*(((gamma+1.0)*Md^2)/((gamma-1.0)*Md^2+2.0));
vd = vs+(rhou/rhod)*(vu-vs);
ud = uu;
RTd = pd/rhod;
cd = (gamma*RTd)^(1/2);
rhoed = pd/(gamma-1);
Ed = rhoed/rhod+(ud^2+vd^2)/2;

% read eigenvectors from file
% this V has conservative vars rho, rhou, rhov, rhoE
V = readmatrix(V_filename);

% read uniform grid files
% MAKE SURE rstrt1_b*.grd IS CORRECT FOR DESIRED GRID SIZE
copyfile(regular_b0_filename,'rstrt1_b0.txt');
copyfile(regular_b1_filename,'rstrt1_b1.txt');
rstrt0_b0 = readcell('rstrt1_b0.txt');
rstrt0_b1 = readcell('rstrt1_b1.txt');
rstrt0_b0( cellfun( @(c) isa(c,'missing'), rstrt0_b0 ) ) = {[]};
rstrt0_b1( cellfun( @(c) isa(c,'missing'), rstrt0_b1 ) ) = {[]};
delete rstrt1_b0.txt
delete rstrt1_b1.txt

% get x and y for blocks 0 and 1
[x0,y0] = getGrid(rstrt0_b0);
[x1,y1] = getGrid(rstrt0_b1);

% perturb mesh at interface
[y0new,rstrt0_b0] = interfaceFun(rstrt0_b0);
[y1new,rstrt0_b1] = interfaceFun(rstrt0_b1);

% perturb flow variables
rstrt0_d0_b0 = initCond(0,rhou,rhou*uu,rhou*vu,rhou*Eu,V,x0,y0);
rstrt0_d0_b1 = initCond(1,rhod,rhod*ud,rhod*vd,rhod*Ed,V,x1,y1);

% convert conservative variables to primitive variables
rstrt0_d0_b0 = varConvert(rstrt0_d0_b0,gamma);
rstrt0_d0_b1 = varConvert(rstrt0_d0_b1,gamma);

% create grid files
writecell(rstrt0_b0,'rstrt0_b0.txt','Delimiter','tab')
writecell(rstrt0_b1,'rstrt0_b1.txt','Delimiter','tab')
movefile('rstrt0_b0.txt','rstrt0_b0.grd');
movefile('rstrt0_b1.txt','rstrt0_b1.grd');

% create flow variable files
dExport(0,rstrt0_b0,'rstrt0_d0_b0.txt',rstrt0_d0_b0)
dExport(1,rstrt0_b1,'rstrt0_d0_b1.txt',rstrt0_d0_b1)
copyfile('rstrt0_d0_b0.txt','rstrt0_d1_b0.txt')
copyfile('rstrt0_d0_b1.txt','rstrt0_d1_b1.txt')

% create x and y files
vExport('rstrt0_v1_b0.txt',x0,y0new)
vExport('rstrt0_v1_b1.txt',x1,y1new)





% plot a variable at target x and y to see shape
xtarget = 10;
ytarget = 0;
ix = find(abs(x1-xtarget)<1e-5);
iy = find(abs(y1-ytarget)<1e-5);

figure
tileplot = tiledlayout(2,2);
title(tileplot,'initial condition at x = ',xtarget)
for var = 1:4
    nexttile
    plot(y1(ix),rstrt0_d0_b1(ix,var),'ko')
    title('var ',var)
    xlabel('y')
end

figure
tileplot = tiledlayout(2,2);
title(tileplot,'initial condition at y = ',ytarget)
for var = 1:4
    nexttile
    plot(x1(iy),rstrt0_d0_b1(iy,var),'ko')
    title('var ',var)
    xlabel('x')
end






function [x,y] = getGrid(table)
npnt = cell2mat(table(1,2));
x = zeros(npnt,1);
y = zeros(npnt,1);
for i = 1:npnt
    x(i) = cell2mat(table(i+1,2));
    y(i) = cell2mat(table(i+1,3));
end
end


function [ynew,table] = interfaceFun(table)
ynew = zeros(cell2mat(table(1,2)),1);
% get dx and L from file
dx = cell2mat(table(6,2));
L = cell2mat(table(3,2));
% interface purturbation amplitude
A = dx/10; 
for i = 1:cell2mat(table(1,2))
    ynew(i) = cell2mat(table(i+1,3));
    if ynew(i) == 0
        % uniform perturbation
%         ynew(i) = A;
        
        % sin perturbation
        ynew(i) = A * sin(2*pi*(cell2mat(table(i+1,2)))/L);
    end
    table(i+1,3) = {ynew(i)};
end
end


function [IC] = initCond(block,rho,rhou,rhov,rhoE,V,x,y)
% get dx and L from file
dx = x(5);
L = x(2);
IC = zeros(length(x),4);
IC(:,1) = rho;
IC(:,2) = rhou;
IC(:,3) = rhov;
IC(:,4) = rhoE;

Vcol = 13;
% perturb block 1
if block == 1
    figure
    tileplot = tiledlayout(2,2);
    title(tileplot,'conservative variable eigenvalues')
    for var = 1:4
        realVtemp = real(V(var:4:end-1,Vcol)); % stores eigenvector for current variable
        nexttile
        plot(realVtemp)
        title('var ',var)
        xlabel('y')
        realVtemp = flip(realVtemp); % positions shock at index 1
        for i = 1:length(x)
            vidx = round(abs(y(i))/dx + 1); % get index for eigenvector
            IC(i,var) = IC(i,var) + 2*realVtemp(vidx)*cos(2*pi*x(i)/L);
        end
    end
end
end


function [data] = varConvert(data,gamma)
for i = 1:size(data,1)
    rho = data(i,1);
    rhou = data(i,2);
    rhov = data(i,3);
    rhoE = data(i,4);
    p = rho*(gamma-1)*(rhoE/rho - ( (rhou/rho)^2+(rhov/rho)^2  )/2);
    u = rhou/rho;
    v = rhov/rho;
    RT = p/rho;
    data(i,1) = p;
    data(i,2) = u;
    data(i,3) = v;
    data(i,4) = RT;
end
end

function [] = dExport(block,grdFile,filename,data)
npnt = cell2mat(grdFile(1,2));
nseg = cell2mat(grdFile(1,4));
ntri = cell2mat(grdFile(1,6));

fid = fopen(filename,'w');
% create header
fprintf(fid,'p0 = 1\n');
fprintf(fid,'npnt = %d, nseg = %d, ntri = %d\n',npnt,nseg,ntri);
fprintf(fid,'END OF HEADER\n');
% add data
for i = 1:size(data,1)
    fprintf(fid,'%e %e %e %e\n',data(i,1),data(i,2),data(i,3),data(i,4));
end
% create footer
if block == 0
    fprintf(fid,'b0_s1 shock\n');
    fprintf(fid,'p0: 1\n');
    fprintf(fid,'b0_s2 plain\n');
    fprintf(fid,'b0_s3 inflow\n');
    fprintf(fid,'b0_s2 plain\n');
    fprintf(fid,'b0_v1 hp_deformable_free_pnt\n');
    fprintf(fid,'b0_v1 hp_deformable_free_pnt\n');
elseif block == 1
    fprintf(fid,'b1_s5 characteristic\n');
    fprintf(fid,'b1_s4 plain\n');
    fprintf(fid,'b1_s1 shock\n');
    fprintf(fid,'p0: 1\n');
    fprintf(fid,'b1_s4 plain\n');
    fprintf(fid,'b1_v1 hp_deformable_free_pnt\n');
    fprintf(fid,'b1_v1 hp_deformable_free_pnt\n');
end
fclose(fid);

end

function [] = vExport(filename,x,y)
npnt = length(x);
fid = fopen(filename,'w');
fprintf(fid,'%i:\n',npnt);
for i = 1:npnt
    fprintf(fid,'%i:%e %e\n',i-1,x(i),y(i));
end
fclose(fid);
end



