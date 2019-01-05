addpath(genpath('./'))
load('wholewn.mat');
nG = size(yy, 2);
sigmaIJN = zeros(10000,nG*(nG-1)/2);
sigmaIJT = zeros(10000,nG*(nG-1)/2);
psiIJN = zeros(10000,nG*(nG-1)/2);
psiIJT = zeros(10000,nG*(nG-1)/2);
sigmaII = zeros(10000,nG);
lambdam = zeros(10000,2);
gammaSqrm = zeros(10000,1);
Number = 200;
accept_n = 0;
for i = 101:Number
    flnum = i*100;
    load(['reswn' num2str(flnum) '.mat']);
    for j = 1:100
        index = j + (i-101)*100;
        sigmaIJN(index,:) = res{1,j}(1,:);
        sigmaIJT(index,:) = res{1,j}(2,:);
        sigmaII(index,:) = res{2,j};
        psiIJN(index,:) = res{4,j}(1,:);
        psiIJT(index,:) = res{4,j}(2,:);
        lambdam(index,:) = res{5,j};
        gammaSqrm(index,:) = res{6,j};
    end
  
end

res{3,100}/20000
fname = 'reswn.mat';
save(fname, 'sigmaIJN', 'sigmaIJT', 'psiIJN', 'psiIJT', 'sigmaII', 'lambdam', 'gammaSqrm');

