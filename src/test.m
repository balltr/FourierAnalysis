clc


a = [1-20i, -1+3i]

sort(real(a))

%%
numVar = 4;

rhsm1 = ones(numVar);
rhs = 2*ones(numVar);
rhsp1 = 3*ones(numVar);
Ny = 5;
numVar = length(rhs);

% tmp = repmat({rhs},Ny,1);
% M = blkdiag(tmp{:});

Mat = full(blktridiag(rhs,rhsm1,rhsp1,Ny));
Mat(end-numVar+1:end,:) = 0;

% top BC
s_rhsm1 = 7*ones(numVar);
s_rhs = 8*ones(numVar);
s_rhsp1 = 9*ones(numVar);
Mat(end-numVar+1:end,end-numVar+1:end) = s_rhs
Mat(end-numVar+1:end,end-2*numVar+1:end-numVar) = s_rhsm1
% Mat(end-numVar+1:end,numVar+1:2*numVar) = s_rhsp1;

% Mat(length(Mat)+1,length(Mat)+1) = 0

ytest = 5*ones(numVar,1);
Mat((Ny-1)*numVar+1:Ny*numVar, Ny*numVar+1) = ytest;
Mat(Ny*numVar+1,Ny*numVar+1) = 6;

% Periodic BC
% Mat(1:numVar,end-numVar+1:end) = s_rhsm1;
% Mat(end-numVar+1:end,1:numVar) = rhsp1



