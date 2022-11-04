function [init_FRF_ac,init_WRF]=generatemorefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf)
m = Nt/Ntrf;
n = Nr/Nrrf;
%%

init_FRF_fc = exp(1i*unifrnd(0,2*pi,Nt,Ntrf));
for ind = 1:5
    library =[16 8 4 2 1];
    di = library(ind);
Lftmp = ones(Nt/di,Ntrf/di);
Lf = kron(Lftmp,eye(di));
init_FRF_ac(:,:,ind) = init_FRF_fc .* Lf;
end


init_WRF= zeros(Nr,Nrrf,K);
for k= 1:K
    init_WRF(:,:,k) = exp(1i*unifrnd(0,2*pi,Nr,Nrrf));%init_fully-connected
end
end
