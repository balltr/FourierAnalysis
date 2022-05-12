

V = readmatrix('V_0_upwind_tantheta_0');
Vcol = 2;
realV1 = real(V(1:4:end-1,Vcol));
realV2 = real(V(3:4:end-1,Vcol));
realV3 = real(V(4:4:end-1,Vcol));

figure
plot(realV1)

figure
plot(realV2)

figure
plot(realV3)