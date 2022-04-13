%% %%%%%%%%%%%%%% 2D Euler Equations With Shock On Top Boundary %%%%%%%%%%%%%%%%%%%%%%
% Conservative variable w = rho, rho*u, rho*v, rho*E
clc
% close all
clear


hold on

syms w1 w2 w3 w4 gamma
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
syms gamma rhod ud vd cd rhou uu vu cu 
Ae = subs(Ae,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
Af = subs(Af,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wd = subs(w,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wu = subs(w,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
eEval = subs(e,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
eEvalUpstream = subs(e,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
fEval = subs(f,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(cd^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
fEvalUpstream = subs(f,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});

[Ve,Ee] = eig(Ae);
[Vf,Ef] = eig(Af);

syms gamma Ms cu cd vu vd
fff = 1/(gamma-1) * sqrt( (2*gamma*Ms^2 - (gamma-1)) * ((gamma-1)*Ms^2 + 2) / Ms^2) + (Ms^2-1)/Ms - ...
    (gamma+1)/(2*cu) * (vu - vd + (2*cd)/(gamma-1));
dfdMs_inv = 1/diff(fff,Ms);

% IC upstream
Ms = 3;
gamma = 1.4;
pu = 1; % non-dimensionalized in P and RT
RTu = 1;
cu = gamma^(1/2);
vu = -Ms*cu; % velocity normal to shock
rhou = gamma*pu/(cu^2);
tantheta = 0.1;
uu = -vu*tantheta;

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
dfdMs_inv = double(subs(dfdMs_inv));

adis = 1.0;
%%%%%%%%%%%%%%%%%% nondimensionalize by dx %%%%%%%%%%%%%%%%%%
dx = 1.0;
dy = 1.0;
y = 0:dy:30;
Ny = length(y);
absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 for central diff or 1 for upwind
fdSwitch = 1;
% 0 for full eig spectrum, 1 for eig pairs, 2 for single eig
eigSwitch = 2;
% 0 for shock off, 1 for shock on
shockSwitch = 0;
% 1 to show eigenvector plots
vecSwitch = 0;

if eigSwitch == 0
    maxWaveNum = pi;
    minWaveNum = -pi;
    incWaveNum = pi/30;
    waveNumArray = minWaveNum:incWaveNum:maxWaveNum;

%     % integer divisions of pi if nondimensionalized by kdx
%     denom = -10;
%     waveNumArray = zeros(2*abs(denom),1);
%     for i = 1:length(waveNumArray)
%         waveNumArray(i) = pi/denom;
%         denom = denom+1;
%         if denom == 0
%             denom = denom+1;
%         end
%     end
elseif eigSwitch == 1
    eigPair = pi/4; % +\-
%     eigPair2 = pi/2;
%     waveNumArray = [-eigPair, eigPair, -eigPair2, eigPair2];
    waveNumArray = [-eigPair, eigPair];
elseif eigSwitch == 2
    waveNumArray = pi/100;
%     waveNumArray = pi/8;
end

CM = turbo(length(waveNumArray));
colorIdx = 0;

for ii = 1:length(waveNumArray)
    kdx = waveNumArray(ii);
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
        Mat2(end-4:end-1,end) = 2/dy*(wd-wu);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % interior
        rhsp1 = -Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
        rhs = -(Ae*ddx-adis*dx*absAe/2*d2dx) - 2*(adis*dy*absAf/2)/dy^2;
        rhsm1 = Af/(2*dy) + (adis*dy*absAf/2)/dy^2;
        Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
        Mat(end-numVar+1:end,:) = 0; % zeros bottom entries for shock boundary elements

        % Moving shock BC on top
        wRHSM1 = 1/dy*(Af+absAf);
        wRHS = -(Ae*ddx-adis*absAe*dx/2*d2dx) - 1/dy*(Af+absAf);
        %   wRHSP1 = zeros(numVar);

        Mat(end-numVar+1:end,end-numVar+1:end) = wRHS;
        Mat(end-numVar+1:end,end-2*numVar+1:end-numVar) = wRHSM1;
        %   Mat(end-numVar+1:end,numVar+1:2*numVar) = wRHSP1;

        if fdSwitch == 0
            % central difference
            ddx_y = ddx;
        elseif fdSwitch == 1
            % upwind (backward difference) based off of tan flow direction
            if uu >= 0
                ddx_y = ddx - dx/2 * d2dx;
            else
                ddx_y = ddx + dx/2 * d2dx;
            end
          
        end

        % for y'
        yRHS = 2/dy*(eEvalUpstream-eEval)*ddx;
        Mat((Ny-1)*numVar+1:Ny*numVar, Ny*numVar+1) = yRHS;
        
        % equation of interface motion
        wpshock = zeros(5,1);
        % w1'
        wpshock(1) = dfdMs_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)^2)*( (wd(2)^2+wd(3)^2)/wd(1) - wd(4)) - wd(3)/wd(1)^2 );
        % w2'
        wpshock(2) = dfdMs_inv*(gamma+1)/2 * ( -gamma*wd(2)/(cd*wd(1)^2) );
        % w3'
        wpshock(3) = dfdMs_inv*(gamma+1)/2 * ( 1/wd(1) - gamma*wd(3)/(cd*wd(1)^2) );
        % w4'
        wpshock(4) = dfdMs_inv*(gamma+1)/2 * ( gamma/(cd*wd(1)) );
        % dy'dxi
        wpshock(5) = ddx_y * 1/dx * ( dfdMs_inv*(gamma+1)/2*(uu - wd(2)/(wd(1))) - uu );
        Mat(Ny*numVar+1,Ny*numVar-3:Ny*numVar+1) = wpshock;
    end
    
    [V,E_mat] = eig(Mat,Mat2);
    E_vec = diag(E_mat);

    plot(real(E_vec),imag(E_vec),'color',CM(colorIdx,:),'marker','x','LineStyle','none')
    hold on
end
if shockSwitch == 1
    if fdSwitch == 1
        title('2D Euler Eigenvalues With Shock, Upwind')
    elseif fdSwitch == 0
        title('2D Euler Eigenvalues With Shock, Central Difference')
    end
elseif shockSwitch == 0
    title('2D Euler Eigenvalues, No Shock')
end
xlabel('Re'), ylabel('Im')
axis equal
grid on
colormap('turbo')
colorbar('Ticks',[])

% % sort eigenvalues in descending order
% [~, ind] = sort(real(E_vec),'descend');
% E_mat = E_mat(ind, ind);
% V = V(:, ind);
% 
% if eigSwitch == 2
%     shockEig = E_vec(1)
% end

% plot eigenvectors
if vecSwitch == 1
    for col = length(E_mat):-1:length(E_mat)-10
        figure
        tiledlayout(2,2)

        nexttile
        plot(real(V(1:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(1:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho')

        nexttile
        plot(real(V(2:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(2:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho u')

        nexttile
        plot(real(V(3:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(3:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho v')

        nexttile
        plot(real(V(4:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold on
        plot(imag(V(4:4:end-shockSwitch,col)),'LineWidth',1.5)
        hold off
        title('\rho E')
        
        if shockSwitch == 1
            txt=['eigenvalue = ',num2str(E_mat(col,col)), '     ',num2str(V(end,col))];
        else
            txt=['eigenvalue = ',num2str(E_mat(col,col))];
        end
        sgtitle(txt)
    end
end


% if shockSwitch == 1
%     figure
%     plot(real(V(end,:)))
%     hold on
%     plot(imag(V(end,:)))
%     hold off
%     legend('Re(V)','Im(V)')
% end