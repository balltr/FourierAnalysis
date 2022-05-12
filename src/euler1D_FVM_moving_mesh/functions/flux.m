function [Fi] = flux(wd,dxdt)
global gam wu wdinit
wdBC = cat(2, wd, wdinit); % add ghost cell on right side with value of wdinit
Nx = size(wd,2);
Fi = zeros(3,Nx+1);

% %%%%%%%%%% symbolic version is slow %%%%%%%%%%
% syms w1 w2 w3 
% F = [w2;
%      w2^2/w1+(gam-1)*(w3-0.5*w2^2/w1);
%      w2/w1*(w3+(gam-1)*(w3-0.5*w2^2/w1))];
% J = jacobian(F,[w1,w2,w3]);
% [V,E] = eig(J);
% absA = V*abs(E)*inv(V);
% myA = matlabFunction(absA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fu = flux2(wu);

for i = 1:Nx+1
    if i == 1
        Fi(:,i) = Fu - dxdt(i)*wu;
    else
        % symbolic
        % absAn = myA(wdBC(1,i),wdBC(2,i),wdBC(3,i)); 
        
        %%%%%%%%% analytic version is fast %%%%%%%%%
        cd = sqrt(gam*(gam-1)/wdBC(1,i) * (wdBC(3,i) - ((wdBC(2,i))^2)/(2*wdBC(1,i))));
        ud = wdBC(2,i)/wdBC(1,i);
        E = [ud,    0,       0;
             0,     ud+cd,   0;
             0,     0,       ud-cd];
        V(:,1) = [2/ud^2; 2/ud; 1];
        V(:,2) = [(2*(gam - 1))/(gam*ud^2 - 2*cd*ud + 2*cd^2 - ud^2 + 2*cd*gam*ud);
                  (2*(cd + ud)*(gam - 1))/(gam*ud^2 - 2*cd*ud + 2*cd^2 - ud^2 + 2*cd*gam*ud); 1];
        V(:,3) = [(2*(gam - 1))/(2*cd*ud + gam*ud^2 + 2*cd^2 - ud^2 - 2*cd*gam*ud);
                  -(2*(cd - ud)*(gam - 1))/(2*cd*ud + gam*ud^2 + 2*cd^2 - ud^2 - 2*cd*gam*ud); 1];
        absAn = V*abs(E)*inv(V);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Fi(:,i) = 0.5*( flux2(wdBC(:,i))+flux2(wdBC(:,i-1)) ) ...
                    - 0.5*absAn*( wdBC(:,i)-wdBC(:,i-1) ) - dxdt(i)*0.5*( wdBC(:,i)+wdBC(:,i-1) );
    end
    end
end