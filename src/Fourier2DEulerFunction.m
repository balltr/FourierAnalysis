function [E_vec,V,uu,Vnorm] = Fourier2DEulerFunction(tantheta,MuINPUT,kdx,upwindSwitch,domainDepth)
% nondimensionalized by dx

% inputs:
    % 1) tan of incoming flow angle
    % 2) upstream Mach number
    % 3) nondimensional wave number
    % 4) upwinding of shock motion equation on or off (1 or 0)
    % 5) depth of domain. Use whole numbers

% returns:
    % 1) eigenvalues sorted by magnitude of shock motion contribution
    % 2) eigenvectors normalized and sorted by magnitude of shock motion contribution
    % 3) tangential velocity
    % 4) magnitude of shock motion contribution

gamma = 1.4;
syms w1 w2 w3 w4 Mu rhod ud vd c_d rhou uu vu cu 
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

Ae = subs(Ae,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(c_d^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
Af = subs(Af,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(c_d^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wd = subs(w,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(c_d^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
wu = subs(w,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
eEval = subs(e,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(c_d^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
eEvalUpstream = subs(e,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});
    % fEval is not used, but here for completeness
% fEval = subs(f,{w1,w2,w3,w4},{rhod,rhod*ud,rhod*vd,rhod*(c_d^2/(gamma*(gamma-1)) +1/2*(ud^2+vd^2))});
% fEvalUpstream = subs(f,{w1,w2,w3,w4},{rhou,rhou*uu,rhou*vu,rhou*(cu^2/(gamma*(gamma-1)) +1/2*(uu^2+vu^2))});

fM = 1/(gamma-1) * sqrt( (2*gamma*Mu^2 - (gamma-1)) * ((gamma-1)*Mu^2 + 2) / Mu^2) + (Mu^2-1)/Mu - ...
    (gamma+1)/(2*cu) * (vu - vd + (2*c_d)/(gamma-1));
dfdMu_inv = 1/diff(fM,Mu);

%%%%%%%%%%%%%%% Use analytic solution to avoid division by zero in eigenvectors %%%%%%%%%%%%%%%
Ee = [ud,  0,       0,       0;
       0, ud,       0,       0;
       0,  0, c_d + ud,       0;
       0,  0,       0, ud - c_d];
Ve(:,1) = [2*vd;2*ud*vd;vd^2 - ud^2;0];
Ve(:,2) = [2;2*ud;0;ud^2 - vd^2];
Ve(:,3) = [2*gamma - 2;2*(c_d + ud)*(gamma - 1);2*vd*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 + 2*c_d*ud) + 2*c_d^2];
Ve(:,4) = [2*gamma - 2;-2*(c_d - ud)*(gamma - 1);2*vd*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 - 2*c_d*ud) + 2*c_d^2];

Ef = [vd,  0,       0,       0;
       0, vd,       0,       0;
       0,  0, c_d + vd,       0;
       0,  0,       0, vd - c_d];
Vf(:,1) = [2*ud;ud^2 - vd^2;2*ud*vd;0];
Vf(:,2) = [2;0;2*vd;vd^2 - ud^2];
Vf(:,3) = [2*gamma - 2;2*ud*(gamma - 1);2*(c_d + vd)*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 + 2*c_d*vd) + 2*c_d^2];
Vf(:,4) = [2*gamma - 2;2*ud*(gamma - 1);-2*(c_d - vd)*(gamma - 1);(gamma - 1)*(vd^2 + ud^2 - 2*c_d*vd) + 2*c_d^2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% IC upstream
pu = 1; % non-dimensionalized in P and RT
RTu = 1;
cu = sqrt(gamma*RTu);
Mu = MuINPUT;
vu = -Mu*cu; % velocity normal to shock
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
c_d = (gamma*RTd)^(1/2);


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
% fEval = double(subs(fEval));
% fEvalUpstream = double(subs(fEvalUpstream));
dfdMu_inv = double(subs(dfdMu_inv));

adis = 1.0;
%%%%%%%%%%%%%%%%%% nondimensionalize by dx %%%%%%%%%%%%%%%%%%
dx = 1.0;
dy = 1.0;
y = 0:dy:domainDepth;
Ny = length(y);
absAe = Ve*abs(Ee)*inv(Ve);
absAf = Vf*abs(Ef)*inv(Vf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ddx = 1i*sin(kdx)/dx;
d2dx = -(2-2*cos(kdx))/dx^2;


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
%   wRHSP1 = zeros(numVar); % not necessary but here for completeness

Mat((Ny-1)*numVar+1:Ny*numVar,(Ny-1)*numVar+1:Ny*numVar) = wRHS;
Mat((Ny-1)*numVar+1:Ny*numVar,(Ny-2)*numVar+1:(Ny-1)*numVar) = wRHSM1;

if upwindSwitch == 0
    % central difference
    ddx_y = ddx;
elseif upwindSwitch == 1
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
wpshock(1) = dfdMu_inv*(gamma+1)/2 * ( gamma/(c_d*wd(1)^2)*( (wd(2)^2+wd(3)^2)/wd(1) - wd(4)) - wd(3)/wd(1)^2 );
% w2'
wpshock(2) = dfdMu_inv*(gamma+1)/2 * ( -gamma*wd(2)/(c_d*wd(1)^2) );
% w3'
wpshock(3) = dfdMu_inv*(gamma+1)/2 * ( 1/wd(1) - gamma*wd(3)/(c_d*wd(1)^2) );
% w4'
wpshock(4) = dfdMu_inv*(gamma+1)/2 * ( gamma/(c_d*wd(1)) );
% dy'dxi
wpshock(5) = ddx_y * ( dfdMu_inv*(gamma+1)/2*(uu - wd(2)/wd(1)) - uu );
Mat(Ny*numVar+1,Ny*numVar-3:Ny*numVar+1) = wpshock;


[V,E_mat] = eig(Mat2\Mat);
E_vec = diag(E_mat);




% normalize eigenvectors
Vnorm = zeros(1,size(V,2));
for i = 1:size(Vnorm,2)
    Vnorm(i) = abs(V(end,i))/norm(V(:,i),2);
end

% sort eigenvectors in descending order
[Vnorm, ind] = sort((Vnorm),'descend');


V = V(:, ind);
% E_mat = E_mat(ind, ind);
E_vec = E_vec(ind);




end

