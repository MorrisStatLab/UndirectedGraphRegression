%%%sigma_value can be chosen for any value that gives acceptance rate between 0.2 and 0.3
sigma_value = 0.075;
%%%import data and calculated Mm
addpath(genpath('./')) %%% add the current path or your running path
load('Profile.mat')
%%matlab version of normal gamma graphical regression
N = size(yy, 1);
nG = size(yy, 2);
nX = size(xx, 2);
accept_iter = zeros(nX,1);
idsigmaIJ = zeros(nG, nG);
h = 0;
for i = 1:(nG-1)
    for j = (i+1):nG
        h = h + 1;
        idsigmaIJ(i,j) = h;
        idsigmaIJ(j,i) = h;
    end
end
%%an initialization
sigmaIJ = zeros(nX, nG*(nG-1)/2);
sigmaII = ones(nG,1);
psiIJ = zeros(nX, nG*(nG-1)/2);
psiIJ(:,:) = 0.25;
lambdam = ones(nX,1);
gammaSqrm = 0.5;
sigma_lambda = zeros(nX,1);
sigma_lambda(:,:) = sigma_value;
niteration = 20000;
res = cell(6,niteration/100);
%%finish the initialization
for t = 1:niteration
    t
    disp('update sigmaIJ')
    for i = 1:(nG-1)
        for j = (i+1):nG
        h = idsigmaIJ(i,j);
        psiIJm = diag(1./psiIJ(:,h));
        S1 = diag((yy(:,j).^2)./sigmaII(i,:)+(yy(:,i).^2)./sigmaII(j,:));
        idi = idsigmaIJ(i,[1:i-1, i+1:j-1, j+1:nG]);
        idj = idsigmaIJ(j,[1:i-1, i+1:j-1, j+1:nG]);
        S2 = 2*yy(:,i).*yy(:,j);
        tmp = diag(xx*sigmaIJ(:,idi)*yy(:,[1:i-1, i+1:j-1, j+1:nG]).');
        S2 = S2 + yy(:,j).*tmp./sigmaII(i,:);
        tmp = diag(xx*sigmaIJ(:,idj)*yy(:,[1:i-1, i+1:j-1, j+1:nG]).');
        S2 = S2 + yy(:,i).*tmp./sigmaII(j);
        %%use a cholesky to calculate inverse for Normal_var
        AA = psiIJm + xx.'*S1*xx;
        RR = chol(AA);
        RRI = RR\speye(size(RR));
        Normal_var = RRI*RRI.';
        Normal_mu = -Normal_var*xx.'*S2;
        %%sample
        sigmaIJ(:,h) = mvnrnd(Normal_mu,Normal_var,1);
        end
    end
    disp('update sigmaII')
    gig1_lambda = N/2 + 1;
    for i = 1:nG
        gig1_psi = sum(yy(:,i).^2);
        idg = idsigmaIJ(i, [1:i-1, i+1:nG]);
        gig1_chi = sum(diag(xx*sigmaIJ(:,idg)*yy(:,[1:i-1,i+1:nG]).').^2);
        sigmaII(i,:) = randraw('gig', [gig1_lambda, gig1_chi, gig1_psi]);
    end
    %%this is a precision check for gig
    lambdam_l = find(lambdam <= 0.5); lln = size(lambdam_l);lln = lln(1);
    if lln > 0
       for l = 1:lln
           ll = lambdam_l(l,:);
           temp = sigmaIJ(ll,:);
           temp(abs(temp) <= sqrt(1e-10)) = sign(temp(abs(temp) <= sqrt(1e-10)))*sqrt(1e-10);
           sigmaIJ(ll,:) = temp;
       end
    end
    
    disp('update psiIJ')
    for m = 1:nX
        gig2_lambda = lambdam(m,:) - 0.5;
        gig2_psi = 1/gammaSqrm;
        for i = 1:(nG-1)
            for j = (i+1):nG
                h = idsigmaIJ(i,j);
                gig2_chi = sigmaIJ(m,h)^2;
                tmpp =  randraw('gig', [gig2_lambda, gig2_chi, gig2_psi]);
                if tmpp < sqrt(1e-10);
                   tmpp = sqrt(1e-10)*2;
                end
                psiIJ(m,h) = tmpp;
            end
        end
    end
    
    %%metropolis step
    disp('update lambda')
    p = nG*(nG-1)/2;
    z_lambda = normrnd(0,1,nX,1);
    lambda_pv = exp(sigma_lambda.*z_lambda).*lambdam;
    accept = zeros(nX,1);
    for l = 1:nX
        logR_lambda = log(exppdf(lambda_pv(l), 1)) - log(exppdf(lambdam(l), 1));
        logR_lambda = logR_lambda + p*(log(gamma(lambdam(l,:))) - log(gamma(lambda_pv(l,:))));
        logR_lambda = logR_lambda + (lambda_pv(l,:) - lambdam(l,:))*(-p*log(2*gammaSqrm) + sum(log(psiIJ(l,:))));
        logR_lambda = logR_lambda + log(lambda_pv(l,:)) - log(lambdam(l,:));
        R_lambda = exp(logR_lambda);
        u = rand;
        if u < R_lambda 
           lambdam(l,:) = lambda_pv(l,:);
           accept(l,:) = 1;
        else
            accept(l,:) = 0;
        end
    end
        accept_iter = accept_iter + accept;
       
        disp('update gammaSqr')
        e_star_gamma = 2 + p*sum(lambdam);
        f_star_gamma = 1/(sum(Mm)/(2*sum(lambdam)) + sum(sum(psiIJ))/2);    
        gammaSqrm = 1/gamrnd(e_star_gamma, f_star_gamma, 1, 1);
    
    st = mod(t-1,100)+1; 
    res{1,st} = sigmaIJ;
    res{2,st} = sigmaII;
    res{3,st} = accept_iter;
    res{4,st} = psiIJ;
    res{5,st} = lambdam;
    res{6,st} = gammaSqrm;
    accept_iter/t
    if mod(t,100) == 0
       fname = strcat('reswn',num2str(t),'.mat');
       save(fname, 'res');
    end
    
end







