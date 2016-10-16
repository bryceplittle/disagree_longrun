function [PP] = solve_condvcov(A,B,C1,C2,D,n_shk,n_Y,kbar,n_hist,jlag)

% pre-allocate matrices
AA  = zeros([kbar+1+n_shk,kbar+1+n_shk,n_hist]);
BB  = zeros([kbar+1+n_shk,n_shk+1,n_hist]);
CC1 = zeros([n_Y,kbar+1+n_shk,n_hist]);
CC2 = zeros([n_Y,kbar+1+n_shk,n_hist]);
DD  = zeros([n_Y,n_shk+1,n_hist]);
PP  = zeros([kbar+1+n_shk,kbar+1+n_shk,n_hist]);
pp  = zeros([kbar+1+n_shk,kbar+1+n_shk,n_hist]);
KK  = zeros([kbar+1+n_shk,n_Y,n_hist]);
LL  = zeros([n_Y,n_Y,n_hist]);

for j=1:n_hist
    AA(1:kbar+1,1:kbar+1,j)     = A(:,:,j);
    BB(1:kbar+1,1:n_shk,j)      = B(:,:,j);
    BB(kbar+2:end,1:n_shk,j)    = eye(n_shk,n_shk);
    CC1(:,1:kbar+1,j)           = C1(:,:,j);
    CC1(:,kbar+2:end,j)         = D(:,1:end-1,j);
    CC2(:,1:kbar+1,j)           = C2(:,:,j);
    DD(:,end,j)                 = D(:,end,j);
    pp(:,:,j)                   = BB(:,:,j)*BB(:,:,j)';
    PP(:,:,j)                   = AA(:,:,j)*pp(:,:,j)*AA(:,:,j)'+BB(:,:,j)*BB(:,:,j)';
end

PP_old  = PP;
count   = 1;
dif     = inf;

% kalman filter
while dif>1e-9 && count<10000
    for j=1:n_hist
        PP(:,:,j)    = AA(:,:,j)*pp(:,:,jlag(j)+1)*AA(:,:,j)'+BB(:,:,j)*BB(:,:,j)';
        LL(:,:,j)    = (CC1(:,:,j)*AA(:,:,j)+CC2(:,:,j))*pp(:,:,jlag(j)+1)*(CC1(:,:,j)*AA(:,:,j)+CC2(:,:,j))'+...
            (CC1(:,:,j)*BB(:,:,j)+DD(:,:,j))*(CC1(:,:,j)*BB(:,:,j)+DD(:,:,j))';
        KK(:,:,j)    = (AA(:,:,j)*pp(:,:,jlag(j)+1)*(CC1(:,:,j)*AA(:,:,j)+CC2(:,:,j))'+...
            BB(:,:,j)*BB(:,:,j)'*CC1(:,:,j)'+BB(:,:,j)*DD(:,:,j)')/LL(:,:,j);
        pp(:,:,j)    = PP(:,:,j)-KK(:,:,j)*LL(:,:,j)*KK(:,:,j)';
    end
    dif = max(max(max(abs(PP-PP_old))));
    PP_old = PP;
    count=count+1;
end

% clearvars -except PP;

end

