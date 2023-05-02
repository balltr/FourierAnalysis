%% Vary tantheta
clc
clearvars
% close all

upwind = 1;
scale_eig_switch = 1;

tanthetaArray = 0.1:0.1:2;

MuArray = [2,3,4];

shockEigArray = zeros(length(MuArray),length(tanthetaArray),2);

% expectedOmega = zeros(length(tanthetaArray),1);
for i = 1:length(MuArray)
    for j = 1:length(tanthetaArray)

        tantheta = tanthetaArray(j);
        Mu = MuArray(i);
        kdxDenom = 10; % only use positive numbers
        kdx = pi/kdxDenom;

        [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind);
        if scale_eig_switch == 1
            % nondim eig value is \lambda / (kdx*sqrt(uu^2+(Mu*cu)^2))
%             E_vec = E_vec/(kdx*sqrt(uu^2+1.4*Mu^2));
            E_vec = E_vec/(kdx*(uu+sqrt(1.4)*Mu));
        end
        expectedOmega = -(2*pi)/( 2*kdxDenom / uu);

        if Vnorm(2)/Vnorm(1) > 0.15 && Vnorm(3)/Vnorm(2) < 0.8
            % two dominant shock modes
            if imag(E_vec(1)) > imag(E_vec(2))
                E1 = E_vec(1);
                E2 = E_vec(2);
            else
                E1 = E_vec(2);
                E2 = E_vec(1);
            end

            % plot 1: upper shock eig
            shockEigArray(i,j,1) = E1;
            
            % plot 2: lower shock eig
            shockEigArray(i,j,2) = E2;
        else
            % one dominant shock mode
            EsortedDecay = sort(real(E_vec),'descend');
            
            if imag(E_vec(1)) > expectedOmega % above expected omega
                if real(E_vec(1)) == real(EsortedDecay(1)) % if mode is slowest decaying
                    % export to plot 2
                    shockEigArray(i,j,1) = NaN + 1i*NaN;
                    shockEigArray(i,j,2) = E_vec(1);
                else % if mode is not slowest decaying
                    
                    % export to plot 1
                    shockEigArray(i,j,1) = E_vec(1);
                    shockEigArray(i,j,2) = NaN + 1i*NaN;
                    
                end
            else % below expected omega
                % export to plot 2
                shockEigArray(i,j,1) = NaN + 1i*NaN;
                shockEigArray(i,j,2) = E_vec(1);
            end
        end
        
    %     expectedOmega(i) = -2*pi/(20/uu);
    end
end

for n = 1:length(MuArray)
    legendInfo{n} = ['Mu = ', num2str(MuArray(n))];
end

% plot imaginary component for E1
figure
for n = 1:length(MuArray)
    plot(tanthetaArray,imag(shockEigArray(n,:,1)))
    hold on
end
if upwind == 1
    title(['Im(upper eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding on'])
else
    title(['Im(upper eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding off'])
end
xlabel('tan(\theta)')
ylabel('Im(eig)')
legend(legendInfo)
grid on
hold off

% plot real component for E1
figure
for n = 1:length(MuArray)
    plot(tanthetaArray,real(shockEigArray(n,:,1)))
    hold on
end
if upwind == 1
    title(['Re(upper eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding on'])
else
    title(['Re(upper eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding off'])
end
xlabel('tan(\theta)')
ylabel('Re(eig)')
legend(legendInfo)
grid on
hold off


% plot imaginary component for E2
figure
for n = 1:length(MuArray)
    plot(tanthetaArray,imag(shockEigArray(n,:,2)))
    hold on
end
if upwind == 1
    title(['Im(lower eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding on'])
else
    title(['Im(lower eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding off'])
end
xlabel('tan(\theta)')
ylabel('Im(eig)')
legend(legendInfo)
grid on
hold off

% plot real component for E2
figure
for n = 1:length(MuArray)
    plot(tanthetaArray,real(shockEigArray(n,:,2)))
    hold on
end
if upwind == 1
    title(['Re(lower eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding on'])
else
    title(['Re(lower eigenvalue) vs tan(\theta),    k\Deltax = \pi/',num2str(kdxDenom),',    upwinding off'])
end
xlabel('tan(\theta)')
ylabel('Re(eig)')
legend(legendInfo)
grid on
hold off


% figure
% plot(tanthetaArray,expectedOmega,'r--')


% plot(tanthetaArray,imag(shockEigArray(:,1)),'k-')
% hold on
% plot(tanthetaArray,imag(shockEigArray(:,2)),'r-')
% plot(tanthetaArray,imag(shockEigArray(:,3)),'b-')

% % plot([10,7.5,5,2.5],[-11.1701,-8.3776,-5.7652,-2.6808],'bd') % upwind
% plot([10,7.5,5,2.5],[-11.1701,-8.2487,-5.2565,-2.3111],'bd') % central
% hold off
% legend('Expected','Fourier','2D FEM')


% figure
% plot(tanthetaArray,real(shockEigArray(:,1)),'k-')
% hold on
% plot(tanthetaArray,real(shockEigArray(:,2)),'r-')
% plot(tanthetaArray,real(shockEigArray(:,3)),'b-')
% title('Re(eig) vs tan(\theta),    k\Deltax = \pi/10,    upwinding off')
% xlabel('tan(\theta)')
% ylabel('Re(eig)')
% grid on
% hold off
% legend(legendInfo)




%% Mu scaling
clc
clearvars
% close all

upwind = 0;

MuArray = 1:1:100;
tantheta = 2.0;
kdx = pi/2;
domainDepth = pi/abs(kdx)*8;

shockEigArray = zeros(length(MuArray),1);
% expectedOmega = zeros(length(MuArray),1);
for i = 1:length(MuArray)
    [E_vec,V,uu] = Fourier2DEulerFunction(tantheta,MuArray(i),kdx,upwind,domainDepth);

    shockEigArray(i) = E_vec(1);
%     expectedOmega(i) = -2*pi/(20/uu);
end


figure;
% plot(MuArray,expectedOmega,'r--')
% hold on
plot(MuArray,abs(imag(shockEigArray)),'ko')
title(['k\Deltax = \pi/10,    tan(\theta) = ',num2str(tantheta)])
xlabel('M')
ylabel('Im(\lambda)')
grid on


%% tantheta scaling
clc
clearvars
% close all

tanThetaArray = 0:0.1:10;
Mu = 3;
kdx = pi/10;
domainDepth = pi/abs(kdx)*8;

shockEigArray = zeros(length(tanThetaArray),1);
for i = 1:length(tanThetaArray)
    [E_vec,V,uu] = Fourier2DEulerFunction(tanThetaArray(i),Mu,kdx,0,domainDepth);
    shockEigArray(i) = E_vec(1);
end

figure;
plot(tanThetaArray,(imag(shockEigArray)),'ko')
title(['k\Deltax = \pi/10,    M = ',num2str(Mu)])
xlabel('tan(\theta)')
ylabel('Im(\lambda)')
grid on

figure;
plot(tanThetaArray,(real(shockEigArray)),'ko')
title(['k\Deltax = \pi/10,    M = ',num2str(Mu)])
xlabel('tan(\theta)')
ylabel('Re(\lambda)')
grid on


%% Vary kdx
clc
clearvars
close all

upwind = 1;

if upwind == 1
    
    % integer divisions of pi
    denom = 100;
    kdxArray = zeros(denom,1);
    for i = 1:length(kdxArray)
        kdxArray(i) = pi/denom;
        denom = denom-1;
    end
    
%     kdxArray = linspace(pi/20,pi,50);
    
else
    kdxArray = [];
end
% kdxArray = pi/10;

shockEigArray = zeros(length(kdxArray),1);
expectedOmega = zeros(length(kdxArray),1);
for i = 1:length(kdxArray)

    tantheta = 1.5;
    Mu = 1.8;
    kdx = kdxArray(i);

    [E_vec,V,uu] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind);

    shockEigArray(i) = E_vec(1);
    expectedOmega(i) = -2*pi/((1/kdx*2*pi)/uu);
end


figure
semilogx(kdxArray,expectedOmega,'r--')
hold on
semilogx(kdxArray,imag(shockEigArray),'k-')
title('Im(eig) vs k\Deltax,    tantheta = 10,    upwinding on')
xlabel('k\Deltax')
ylabel('Im(eig)')
grid on
hold off
legend('Expected','Fourier')

figure
semilogx(kdxArray,real(shockEigArray),'k-')
title('Re(eig) vs k\Deltax,    tantheta = 10,    upwinding on')
xlabel('k\Deltax')
ylabel('Re(eig)')
grid on



%% tantheta and M Contour Plot (scaling comparison)
clearvars
% close all

tanThetaArray = 0:0.1:10;
MuArray = 1:0.1:10;
kdx = pi/10;
domainDepth = pi/abs(kdx)*8;

shockEigArray = zeros(length(tanThetaArray),length(MuArray));
for i = 1:length(tanThetaArray)
    for j = 1:length(MuArray)
        [E_vec,V,uu] = Fourier2DEulerFunction(tanThetaArray(i),MuArray(j),kdx,0,domainDepth);
        shockEigArray(i,j) = E_vec(1);
    end
end

[X,Y] = meshgrid(tanThetaArray,MuArray);

%%

figure;
contourf(X',Y',imag(shockEigArray))
title('Im(\lambda),    Not Scaled,    k\Deltax = \pi/10')
xlabel('tan(\theta)')
ylabel('M')
colorbar


figure;
contourf(X',Y',real(shockEigArray))
title('Re(\lambda),    Not Scaled,    k\Deltax = \pi/10')
xlabel('tan(\theta)')
ylabel('M')
colorbar


% scaling
uMatrix = zeros(length(tanThetaArray),length(MuArray));
shockEigArrayScaled = zeros(length(tanThetaArray),length(MuArray));
cu = sqrt(1.4);

for i = 1:length(tanThetaArray)
    for j = 1:length(MuArray)
        uMatrix(i,j) = -vMatrix(i,j)*tanThetaArray(i);
        shockEigArrayScaled(i,j) = shockEigArray(i,j)/(kdx*(uMatrix(i,j)+MuArray(j)*cu));
    end
end

figure;
contourf(X',Y',imag(shockEigArrayScaled))
title('Im(\lambda),    Scaled,    k\Deltax = \pi/10')
xlabel('tan(\theta)')
ylabel('M')
colorbar

figure;
contourf(X',Y',real(shockEigArrayScaled))
title('Re(\lambda),    Scaled,    k\Deltax = \pi/10')
xlabel('tan(\theta)')
ylabel('M')
colorbar













%% Compare upwinding - vary tantheta

clc
clearvars
% close all

% set these parameters
scale_eig_switch = 1; % scale eigenvalues
threshold = 0.4; % threshold to determine dominant shock motion modes
tanthetaArray = 0.0:0.2:10;
% tanthetaArray = [0:0.05:1, 1.1:0.2:10];
Mu = 3;
kdxDenom = 1;
kdx = pi/kdxDenom;
% domainDepth = pi/abs(kdx)*8;
domainDepth = pi/abs(kdx)*100;

upwindArray = [0,1];
cu = sqrt(1.4);
shockEigArray = (NaN + 1i*NaN)*ones(2,length(tanthetaArray),4);

for i = 1:2
    for j = 1:length(tanthetaArray)

        tantheta = tanthetaArray(j);
        upwind = upwindArray(i);
        
        [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind,domainDepth);
        if scale_eig_switch == 1
            E_vec = E_vec/(kdx*(uu+Mu*cu));
        end
        
        if Mu==5 && tantheta==0.2
            % omit this strange case with bad condition number
        else
            if Vnorm(2) < threshold*Vnorm(1)
                % 1 dominant shock mode
                shockEigArray(i,j,1) = E_vec(1);
            elseif Vnorm(3) < threshold*Vnorm(1)
                % 2 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
            elseif Vnorm(4) < threshold*Vnorm(1)
                % 3 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
            elseif Vnorm(5) < threshold*Vnorm(1)
                % 4 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
                shockEigArray(i,j,4) = E_vec(4);
            else
                % no dominant shock modes
            end
        end
    end
end

markerType = {'ro','k*'};

% plot real components
fReal = figure;
figure(fReal)
for n = 1:2
    plot(tanthetaArray,real(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(tanthetaArray,real(shockEigArray(n,:,2)),markerType{n})
    plot(tanthetaArray,real(shockEigArray(n,:,3)),markerType{n})
    plot(tanthetaArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
if kdxDenom ~= 1
    title(['k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
else
    title(['k\Deltax = \pi,    M = ',num2str(Mu)])
end
xlabel('tan($\theta$)','Interpreter','Latex')
ylabel('Re($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on
ylim([-inf 0])

% plot imaginary components
fImag = figure;
figure(fImag)
for n = 1:2
    plot(tanthetaArray,imag(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(tanthetaArray,imag(shockEigArray(n,:,2)),markerType{n})
    plot(tanthetaArray,imag(shockEigArray(n,:,3)),markerType{n})
    plot(tanthetaArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
if kdxDenom ~= 1
    title(['k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
else
    title(['k\Deltax = \pi,    M = ',num2str(Mu)])
end
xlabel('tan($\theta$)','Interpreter','Latex')
ylabel('Im($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on




%% Compare upwinding - vary M

clc
clearvars
% close all

% set these parameters
scale_eig_switch = 1; % scale eigenvalues
threshold = 0.1; % threshold to determine dominant shock motion modes
tantheta = 0.1;
MuArray = 1.2:0.2:20;
% MuArray = 1.5;
kdxDenom = 10;
kdx = pi/kdxDenom;
domainDepth = pi/abs(kdx)*8;
% domainDepth = pi/abs(kdx)*100;

upwindArray = [0,1];
cu = sqrt(1.4);
shockEigArray = (NaN + 1i*NaN)*ones(2,length(MuArray),4);

for i = 1:2
    for j = 1:length(MuArray)

        Mu = MuArray(j);
        upwind = upwindArray(i);
        
        [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind,domainDepth);
        if scale_eig_switch == 1
            E_vec = E_vec/(kdx*(uu+Mu*cu));
        end
        
        if Mu==5 && tantheta==0.2
            % omit this strange case with bad condition number
        else
            if Vnorm(2) < threshold*Vnorm(1)
                % 1 dominant shock mode
                shockEigArray(i,j,1) = E_vec(1);
            elseif Vnorm(3) < threshold*Vnorm(1)
                % 2 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
            elseif Vnorm(4) < threshold*Vnorm(1)
                % 3 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
            elseif Vnorm(5) < threshold*Vnorm(1)
                % 4 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
                shockEigArray(i,j,4) = E_vec(4);
            else
                % no dominant shock modes
            end
        end
    end
end

markerType = {'ro','k*'};

% plot real components
fReal = figure;
figure(fReal)
for n = 1:2
    plot(MuArray,real(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(MuArray,real(shockEigArray(n,:,2)),markerType{n})
    plot(MuArray,real(shockEigArray(n,:,3)),markerType{n})
    plot(MuArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
title(['k\Deltax = \pi/',num2str(kdxDenom),',    tan(\theta) = ',num2str(tantheta)])
xlabel('M','Interpreter','Latex')
ylabel('Re($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on
ylim([-inf 0])

% plot imaginary components
fImag = figure;
figure(fImag)
for n = 1:2
    plot(MuArray,imag(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(MuArray,imag(shockEigArray(n,:,2)),markerType{n})
    plot(MuArray,imag(shockEigArray(n,:,3)),markerType{n})
    plot(MuArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
title(['k\Deltax = \pi/',num2str(kdxDenom),',    tan(\theta) = ',num2str(tantheta)])
xlabel('M','Interpreter','Latex')
ylabel('Im($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on






%% Compare upwinding - vary kdx

clc
clearvars
% close all

% set these parameters
scale_eig_switch = 1; % scale eigenvalues
threshold = 0.1; % threshold to determine dominant shock motion modes
tantheta = 0.0;
Mu = 3;

% % integer divisions of pi
% denom = 20;
% kdxArray = zeros(1,denom);
% for i = 1:denom
%     kdxArray(i) = pi/denom;
%     denom = denom-1;
% end

kdxArray = pi/20 : pi/40 : pi;


domainDepth = 160;
% domainDepth = pi/abs(kdx)*8; 

upwindArray = [0,1];
cu = sqrt(1.4);
shockEigArray = (NaN + 1i*NaN)*ones(2,length(kdxArray),4);

for i = 1:2
    for j = 1:length(kdxArray)

        kdx = kdxArray(j);
        upwind = upwindArray(i);
        
        [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind,domainDepth);
        if scale_eig_switch == 1
            E_vec = E_vec/(kdx*(uu+Mu*cu));
        end
        
        if Mu==5 && tantheta==0.2
            % omit this strange case with bad condition number
        else
            if Vnorm(2) < threshold*Vnorm(1)
                % 1 dominant shock mode
                shockEigArray(i,j,1) = E_vec(1);
            elseif Vnorm(3) < threshold*Vnorm(1)
                % 2 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
            elseif Vnorm(4) < threshold*Vnorm(1)
                % 3 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
            elseif Vnorm(5) < threshold*Vnorm(1)
                % 4 dominant shock modes
                shockEigArray(i,j,1) = E_vec(1);
                shockEigArray(i,j,2) = E_vec(2);
                shockEigArray(i,j,3) = E_vec(3);
                shockEigArray(i,j,4) = E_vec(4);
            else
                % no dominant shock modes
            end
        end
    end
end

markerType = {'ro','k*'};

% plot real components
fReal = figure;
figure(fReal)
for n = 1:2
    plot(kdxArray,real(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(kdxArray,real(shockEigArray(n,:,2)),markerType{n})
    plot(kdxArray,real(shockEigArray(n,:,3)),markerType{n})
    plot(kdxArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
title(['M = ',num2str(Mu),',    tan(\theta) = ',num2str(tantheta)])
xlabel('$k \Delta x$','Interpreter','Latex')
ylabel('Re($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on
ylim([-inf 0])

% plot imaginary components
fImag = figure;
figure(fImag)
for n = 1:2
    plot(kdxArray,imag(shockEigArray(n,:,1)),markerType{n})
    hold on
    plot(kdxArray,imag(shockEigArray(n,:,2)),markerType{n})
    plot(kdxArray,imag(shockEigArray(n,:,3)),markerType{n})
    plot(kdxArray,real(shockEigArray(n,:,4)),markerType{n})
end
hold off
title(['M = ',num2str(Mu),',    tan(\theta) = ',num2str(tantheta)])
xlabel('$k \Delta x$','Interpreter','Latex')
ylabel('Im($\lambda$)','Interpreter','Latex')
legend('upwind off','','','','upwind on','','','','location','best')
grid on






%% upwinding on or off - vary tantheta

clc
clearvars
% close all

% set these parameters
scale_eig_switch = 1; % nondimensionalize eigenvalues by u+Mc
threshold = 0.1; % threshold to determine dominant shock motion modes
tanthetaArray = 1.0:0.5:10;
% tanthetaArray = [0:0.05:1, 1.1:0.2:10];
Mu = 3;
kdxDenom = 10;
kdx = pi/kdxDenom;
domainDepth = pi/abs(kdx)*8;
% domainDepth = pi/abs(kdx)*100;

upwind = 0;
cu = sqrt(1.4);
shockEigArray = (NaN + 1i*NaN)*ones(length(tanthetaArray),4);

for j = 1:length(tanthetaArray)

    tantheta = tanthetaArray(j);

    [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,Mu,kdx,upwind,domainDepth);
    if scale_eig_switch == 1
        E_vec = E_vec/(kdx*(uu+Mu*cu));
    end

    if Mu==5 && tantheta==0.2
        % omit this strange case with bad condition number
    else
        if Vnorm(2) < threshold*Vnorm(1)
            % 1 dominant shock mode
            shockEigArray(j,1) = E_vec(1);
        elseif Vnorm(3) < threshold*Vnorm(1)
            % 2 dominant shock modes
            shockEigArray(j,1) = E_vec(1);
            shockEigArray(j,2) = E_vec(2);
        elseif Vnorm(4) < threshold*Vnorm(1)
            % 3 dominant shock modes
            shockEigArray(j,1) = E_vec(1);
            shockEigArray(j,2) = E_vec(2);
            shockEigArray(j,3) = E_vec(3);
        elseif Vnorm(5) < threshold*Vnorm(1)
            % 4 dominant shock modes
            shockEigArray(j,1) = E_vec(1);
            shockEigArray(j,2) = E_vec(2);
            shockEigArray(j,3) = E_vec(3);
            shockEigArray(j,4) = E_vec(4);
        else
            % no dominant shock modes
        end
    end
end

% plot real components
fReal = figure;
figure(fReal)
plot(tanthetaArray,real(shockEigArray(:,1)),'ko')
hold on
plot(tanthetaArray,real(shockEigArray(:,2)),'ko')
plot(tanthetaArray,real(shockEigArray(:,3)),'ko')
plot(tanthetaArray,real(shockEigArray(:,4)),'ko')
hold off
if kdxDenom ~= 1
    if scale_eig_switch == 1
        title(['scaled,    k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
    else
        title(['not scaled,    k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
    end
else
    if scale_eig_switch == 1
        title(['scaled,    k\Deltax = \pi,    M = ',num2str(Mu)])
    else
        title(['not scaled,    k\Deltax = \pi,    M = ',num2str(Mu)])
    end
end
xlabel('tan($\theta$)','Interpreter','Latex')
ylabel('Re($\lambda$)','Interpreter','Latex')
grid on
ylim([-inf 0])

% plot imaginary components
fImag = figure;
figure(fImag)
plot(tanthetaArray,imag(shockEigArray(:,1)),'ko')
hold on
plot(tanthetaArray,imag(shockEigArray(:,2)),'ko')
plot(tanthetaArray,imag(shockEigArray(:,3)),'ko')
plot(tanthetaArray,imag(shockEigArray(:,4)),'ko')
hold off
hold off
if kdxDenom ~= 1
    if scale_eig_switch == 1
        title(['scaled,    k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
    else
        title(['not scaled,    k\Deltax = \pi/',num2str(kdxDenom),',    M = ',num2str(Mu)])
    end
else
    if scale_eig_switch == 1
        title(['scaled,    k\Deltax = \pi,    M = ',num2str(Mu)])
    else
        title(['not scaled,    k\Deltax = \pi,    M = ',num2str(Mu)])
    end
end
xlabel('tan($\theta$)','Interpreter','Latex')
ylabel('Im($\lambda$)','Interpreter','Latex')
grid on

