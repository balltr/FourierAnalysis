%% %%%%%%%%%%%%%% 2D Euler Equations With Shock On Top Boundary %%%%%%%%%%%%%%%%%%%%%%
% Conservative variable w = rho, rho*u, rho*v, rho*E
clc
clearvars
% close all
% figure
% hold on

syms w1 w2 w3 w4 gamma Mu rhod ud vd cd rhou uu vu cu 
e = [w2;
     w2^2/w1+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w2*w3/w1;
     w2/w1*(w4+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1))];
f = [w3;
     w2*w3/w1;
     w3^2/w1+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1);
     w3/w1*(w4+(gamma-1)*(w4-1/2*(w2^2+w3^2)/w1))];
w = [w1; w2; w3; w4];
 
numVar = length(e);

Ae = jacobian(e,[w1,w2,w3,w4]);
Af = jacobian(f,[w1,w2,w3,w4]);

Ae = subs(Ae,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
Af = subs(Af,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wd = subs(w,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wu = subs(w,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
eEval = subs(e,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
eEvalUpstream = subs(e,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
fEval = subs(f,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
fEvalUpstream = subs(f,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});

fM = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*cu) * (vu - vd + (2*cd)/(gamma-1));
dfdMu_inv = 1/diff(fM,Mu);

%%%%%%%%%%%%%%% Use analytic solution to avoid eigenvectors dividing by zero %%%%%%%%%%%%%%%
Ee = [ud,  0,       0,       0;
       0, ud,       0,       0;
       0,  0, cd + ud,       0;
       0,  0,       0, ud - cd];
Ve(:,1) = [2*vd;2*ud*vd;vd^2 - ud^2;0];
Ve(:,2) = [2;2*ud;0;ud^2 - vd^2];
Ve(:,3) = [2*gamma - 2;2*(cd + ud)*(gamma - 1);2*vd*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 + 2*cd*ud) + 2*cd^2];
Ve(:,4) = [2*gamma - 2;-2*(cd - ud)*(gamma - 1);2*vd*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 - 2*cd*ud) + 2*cd^2];

Ef = [vd,  0,       0,       0;
       0, vd,       0,       0;
       0,  0, cd + vd,       0;
       0,  0,       0, vd - cd];
Vf(:,1) = [2*ud;ud^2 - vd^2;2*ud*vd;0];
Vf(:,2) = [2;0;2*vd;vd^2 - ud^2];
Vf(:,3) = [2*gamma - 2;2*ud*(gamma - 1);2*(cd + vd)*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 + 2*cd*vd) + 2*cd^2];
Vf(:,4) = [2*gamma - 2;2*ud*(gamma - 1);-2*(cd - vd)*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 - 2*cd*vd) + 2*cd^2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% IC downstream
vs = 0.0;
Md = (vu-vs)/cu; % relative to shock
pd = pu*(2*gamma*Md^2-(gamma-1))/(gamma+1);
rhod = rhou*(((gamma+1.0)*Md^2)/((gamma-1.0)*Md^2+2.0));
vd = vs+(rhou/rhod)*(vu-vs);
ud = uu;
RTd = pd/rhod;
cd = (gamma*RTd)^(1/2);

Ae = double(subs(Ae));
Ve = double(subs(Ve));
Ee = double(subs(Ee));
Af = double(subs(Af));
Vf = double(subs(Vf));
Ef = double(subs(Ef));
wd = double(subs(wd));
wu = double(subs(wu));
eEval = double(subs(eEval));
eEvalUpstream = double(subs(eEvalUpstream));
fEval = double(subs(fEval));
fEvalUpstream = double(subs(fEvalUpstream));
dfdMu_inv = double(subs(dfdMu_inv));

adis = 1.0;
%%%%%%%%%%%%%%%%%% nondimensionalize by dx %%%%%%%%%%%%%%%%%%
dx = 1.0;
dy = 1.0;
y = 0:dy:80;
Ny = length(y);
absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 for central diff or 1 for upwind
fdSwitch = 1;
% 0 for full eig spectrum, 1 for eig pairs, 2 for single eig
eigSwitch = 2;
% 0 for shock off, 1 for shock on
shockSwitch = 1;
% 1 to show eigenvector plots
vecSwitch = 1;

if eigSwitch == 0
    kdxArray = -pi : pi/10 : pi;

%     % integer divisions of pi if nondimensionalized by kdx
%     denom = -10;
%     kdxArray = zeros(2*abs(denom),1);
%     for i = 1:length(waveNumArray)
%         kdxArray(i) = pi/denom;
%         denom = denom+1;
%         if denom == 0
%             denom = denom+1;
%         end
%     end
elseif eigSwitch == 1
    eigPair = pi/4; % +\-
%     eigPair2 = pi/2;
%     kdxArray = [-eigPair, eigPair, -eigPair2, eigPair2];
    kdxArray = [-eigPair, eigPair];
elseif eigSwitch == 2
    kdxArray = pi/10;
end

CM = turbo(length(kdxArray));
colorIdx = 0;

for ii = 1:length(kdxArray)
    kdx = kdxArray(ii);
%     %%%%%%%%%%%%%%%%%%% nondimensionalize by kdx %%%%%%%%%%%%%%%%%%
%     dx = abs(kdx/(2*pi));
%     dy = dx;
%     Ny = round(4/dy+1,0);
%     absAe = dx*Ve*abs(Ee)*inv(Ve);
%     absAf = dy*Vf*abs(Ef)*inv(Vf);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    colorIdx = colorIdx + 1;
    
    ddx = 1i*sin(kdx)/dx;
    d2dx = -(2-2*cos(kdx))/dx^2;
    
    if shockSwitch == 0
        rhsp1 = -Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
        rhs = -(Ae*ddx-adis*dx*absAe/2*d2dx) - 2*(adis*dy*absAf/2)/dy^2;
        rhsm1 = Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
  
        Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
  
        % Periodic BC's
        Mat(1:numVar,end-numVar+1:end) = rhsm1;
        Mat(end-numVar+1:end,1:numVar) = rhsp1;
        
        Mat2 = eye(Ny*numVar);
        
    elseif shockSwitch == 1
        
        %%%%%%%%%%% Matrix for d(w,y)/dt %%%%%%%%%%%%%%%%%%%
        Ndof = Ny*numVar+1;
        Mat2 = eye(Ndof);
        Mat2(numVar*(Ny-1)+1:numVar*Ny,end) = 2/dy*(wd-wu);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mat = zeros(Ndof);
        
        % interior
        rhsp1 = -Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
        rhs = -(Ae*ddx-adis*dx*absAe/2*d2dx) - 2*(adis*dy*absAf/2)/dy^2;
        rhsm1 = Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
        Mat_temp = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
        Mat(1:numVar*(Ny-1),1:numVar*Ny) = Mat_temp(1:numVar*(Ny-1),1:numVar*Ny);

        % Moving shock BC on top
        wRHSM1 = 1/dy*(Af+absAf);
        wRHS = -(Ae*ddx-adis*dx*absAe/2*d2dx) + 1/dy*(Af-absAf);
        %   wRHSP1 = zeros(numVar);

        Mat((Ny-1)*numVar+1:Ny*numVar,(Ny-1)*numVar+1:Ny*numVar) = wRHS;
        Mat((Ny-1)*numVar+1:Ny*numVar,(Ny-2)*numVar+1:(Ny-1)*numVar) = wRHSM1;

        if fdSwitch == 0
            % central difference
            ddx_y = ddx;
        elseif fdSwitch == 1
            % upwind (backward difference) based off of tan flow direction
            if uu >= 0
                ddx_y = ddx - dx/2*d2dx;
            else
                ddx_y = ddx + dx/2*d2dx;
            end
          
        end

        % for y'
        yRHS = 2/dy*(eEvalUpstream-eEval)*ddx;
        Mat((Ny-1)*numVar+1:Ny*numVar, Ny*numVar+1) = yRHS;
        
        % equation of interface motion
        wpshock = zeros(5,1);
        % w1'
        wpshock(1) = dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)^2)*( (wd(2)^2+wd(3)^2)/wd(1) - wd(4)) - wd(3)/wd(1)^2 );
        % w2'
        wpshock(2) = dfdMu_inv*(gamma+1)/2 * ( -gamma*wd(2)/(cd*wd(1)^2) );
        % w3'
        wpshock(3) = dfdMu_inv*(gamma+1)/2 * ( 1/wd(1) - gamma*wd(3)/(cd*wd(1)^2) );
        % w4'
        wpshock(4) = dfdMu_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)) );
        % dy'dxi
        wpshock(5) = ddx_y * ( dfdMu_inv*(gamma+1)/2*(uu - wd(2)/wd(1)) - uu );
        Mat(Ny*numVar+1,Ny*numVar-3:Ny*numVar+1) = wpshock;
    end
    
    [V,E_mat,W] = eig(Mat2\Mat);
    E_vec = diag(E_mat);

%     plot(real(E_vec),imag(E_vec),'color',CM(colorIdx,:),'marker','x','LineStyle','none')
%     hold on
end
% if shockSwitch == 1
%     if fdSwitch == 1
%         title('2D Euler Eigenvalues With Shock, Upwind')
%     elseif fdSwitch == 0
%         title('2D Euler Eigenvalues With Shock, Central Difference')
%     end
% elseif shockSwitch == 0
%     title('2D Euler Eigenvalues, No Shock')
% end
% xlabel('Re'), ylabel('Im')
% axis equal
% grid on
% colormap('turbo')
% colorbar('Ticks',[])

% normalize eigenvectors
if shockSwitch == 1
    Vnorm = zeros(1,size(V,2));
    for i = 1:size(Vnorm,2)
        Vnorm(i) = abs(V(end,i))/norm(V(:,i),2);
    end
    % sort eigenvectors in descending order
    [~, ind] = sort((Vnorm),'descend');
%     [~, ind] = sort(real(E_vec),'descend');
    Vnorm = Vnorm(ind);
    V = V(:, ind);
    E_mat = E_mat(ind, ind);
    E_vec = E_vec(ind);
    
    figure
    plot(Vnorm,'x-')
    title('Shock Position Normalized Eigenvectors')
    xlabel('DOF')
end
% plot eigenvalues colored by magnitude of corresponding normalized eigenvector
figure
scatplot = scatter(real(E_vec),imag(E_vec),25,Vnorm,'filled');
scatplot.MarkerEdgeColor = 'k';
axis equal
grid on

% plot eigenvectors
if vecSwitch == 1
    for col = 1:2
        figure
        tiledlayout(2,2)

        nexttile
        plot(real(V(1:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(1:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho')
        xlabel('y')

        nexttile
        plot(real(V(2:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(2:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho u')
        xlabel('y')

        nexttile
        plot(real(V(3:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(3:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho v')
        xlabel('y')

        nexttile
        plot(real(V(4:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(4:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho E')
        xlabel('y')
        
        if shockSwitch == 1
            txt=['eigenvalue = ',num2str(E_mat(col,col))];
        else
            txt=['eigenvalue = ',num2str(E_mat(col,col))];
        end
        sgtitle(txt)
    end
end

% % weighted average of 25 most dominant eigenvectors
% Nvec = Ndof
% VnormSum = sum(Vnorm(1:Nvec));
% Vweight = Vnorm(1:Nvec)/VnormSum;
% weightedOmega = 0;
% for i = 1:Nvec
%     weightedOmega = weightedOmega + abs(imag(E_vec(i)))*Vweight(i);
% end
% weightedOmega





% 
% % wp0 = [0,0,0,0];
% % wp0 = repmat(wp0,length(y),1);
% % wp0(end+1) = 0;
% IC = zeros(size(V,1),1);
% IC(end,1) = 1.0;
% 
% C = V\IC;
% 
% % Pick a time
% t = 0;
% 
% Elam = diag(exp(E_vec*t));
% w = V*Elam*C;
% 
% figure
% tiledlayout(2,2)
% 
% nexttile
% plot(real(w(1:4:end-shockSwitch)),'LineWidth',1.5)
% hold on
% plot(imag(w(1:4:end-shockSwitch)),'LineWidth',1.5)
% hold off
% title('\rho')
% xlabel('y')
% 
% nexttile
% plot(real(w(2:4:end-shockSwitch)),'LineWidth',1.5)
% hold on
% plot(imag(w(2:4:end-shockSwitch)),'LineWidth',1.5)
% hold off
% title('\rho u')
% xlabel('y')
% 
% nexttile
% plot(real(w(3:4:end-shockSwitch)),'LineWidth',1.5)
% hold on
% plot(imag(w(3:4:end-shockSwitch)),'LineWidth',1.5)
% hold off
% title('\rho v')
% xlabel('y')
% 
% nexttile
% plot(real(w(4:4:end-shockSwitch)),'LineWidth',1.5)
% hold on
% plot(imag(w(4:4:end-shockSwitch)),'LineWidth',1.5)
% hold off
% title('\rho E')
% xlabel('y')


