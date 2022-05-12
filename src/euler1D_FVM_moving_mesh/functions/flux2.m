function [F] = flux2(w)
global gam
F = [w(2);
     w(2)^2/w(1)+(gam-1)*(w(3)-0.5*w(2)^2/w(1));
     w(2)/w(1)*(w(3)+(gam-1)*(w(3)-0.5*w(2)^2/w(1)))];
end