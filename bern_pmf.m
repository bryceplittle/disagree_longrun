function pmf = bern_pmf(T,p,k)
# credit: Kris Nimark

pmf = (p^k)*(1-p)^(T-k);

end


