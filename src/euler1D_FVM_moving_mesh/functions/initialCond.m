function [x,wd] = initialCond(x,wd)

% perturb interface and flow using eigenvectors from Fourier2DEuler
V = readmatrix('V_0_upwind_tantheta_0');
Vcol = 64;
epsFactor = 10;
realV1 = flip(real(V(1:4:end-1,Vcol)));
realV2 = flip(-real(V(3:4:end-1,Vcol)));
realV3 = flip(real(V(4:4:end-1,Vcol)));

x(1) = x(1)+epsFactor*real(V(end,Vcol));
for i = 1:size(wd,2)
    wd(1,i) = wd(1,i) + epsFactor*2*realV1(i);
    wd(2,i) = wd(2,i) + epsFactor*2*realV2(i);
    wd(3,i) = wd(3,i) + epsFactor*2*realV3(i);
end

end