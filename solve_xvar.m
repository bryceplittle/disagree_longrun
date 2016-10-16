function [Pi] = solve_xvar(A,C1,C2,D,K,n_hist,jlag)

Pi      = zeros(size(A));
Pi_old  = zeros(size(A));
dif     = inf;

while dif>1e-10
    for jj=1:n_hist
        Pi(:,:,jj) = (A(:,:,jj)-K(:,:,jj)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj)))*Pi(:,:,jlag(jj)+1)*(A(:,:,jj)-...
		K(:,:,jj)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj)))'+K(:,:,jj)*D(:,end,jj)*D(:,end,jj)'*K(:,:,jj)';
    end
    dif     = max(max(max(abs(Pi-Pi_old))));
    Pi_old  = Pi;
end

end

