function [data] = varConvert(data)
global gam
for i = 1:size(data,2)
    rho = data(1,i);
    rhou = data(2,i);
    rhoE = data(3,i);
    p = rho*(gam-1)*(rhoE/rho - ((rhou/rho)^2)/2);
    u = rhou/rho;
    RT = p/rho;
    data(1,i) = p;
    data(2,i) = u;
    data(3,i) = RT;
end
end