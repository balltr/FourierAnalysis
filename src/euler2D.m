% 2D Steady Euler equations
% Conservative variable w = rho, rho u, rho v, rho E
syms w1 w2 w3 w4 positive;
syms gam positive;
e = [w2;w2^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);w2*w3/w1;w2/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
f = [w3;w2*w3/w1;w3^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);w3/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];

Ae = jacobian(e,[w1,w2,w3,w4]);
Af = jacobian(f,[w1,w2,w3,w4]);

syms rho u v gam c positive
Ae = subs(Ae,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
Af = subs(Af,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});

[Ve,Ee] = eig(Ae);
Ve = simplify(Ve)
Ee = simplify(Ee)
Ve_inv = simplify(Ve)

%absAe = simplify(Ve*abs(Ee)*inv(Ve))

return;

C = simplify(Ae\Af);
eig(C);

return


[Vf,Ef] = eig(Af);
absAf = simplify(Vf*abs(Ef)*inv(Vf));

% SEMI-DIAGONAL FORM

cross = inv(Ve)*Af*Ve;
cross = simplify(cross)

syms E u0 v0 rho0 c0

Ve = subs(Ve,{rho,u,v,c},{rho0,u0,v0,c0})
charvar = simplify(inv(Ve)*[rho;rho*u;rho*v;rho*(c^2/(gam*(gam-1))+1/2*(u^2+v^2))]);

pretty(charvar)



%charvarhorizontal = simplify(subs(charvar,v0,0));
% CHARVAR uncoupled = rho/gam*(gam-c^2/c0^2 -gam*(gam-1)*((u0-u)^2+(v0-v)^2)/(2*c0^2))

%pretty(charvarhorizontal)

