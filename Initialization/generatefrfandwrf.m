function [init_FRF_pc, init_FRF_fc,init_WRF_fc,init_WRF_pc,init_FRF_ac,init_WRF_ac]=generatefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf)
m = Nt/Ntrf;
n = Nr/Nrrf;
%%
di = 2; %division index
Lftmp = ones(Nt/di,Ntrf/di);
Lf = kron(Lftmp,eye(di));
%%
m = Nt/Ntrf;
n = Nr/Nrrf;
init_WRF_fc = zeros(Nr,Nrrf,K);
init_WRF_pc = zeros(Nr,Nrrf,K);
init_WRF_ac = zeros(Nr,Nrrf,K);
init_FRF_pc = zeros(Nt,Ntrf);

init_FRF_fc = exp(1i*unifrnd(0,2*pi,Nt,Ntrf));%init_fully-connected
for i = 1 : Ntrf %initialize F_RF
    init_FRF_pc((i-1) * m + 1:i*m,i) = exp( 1i*unifrnd(0,2*pi,m,1)); %init_partially-connected
end
init_FRF_ac = init_FRF_fc .* Lf;

for k= 1:K
    init_WRF_fc(:,:,k) = exp(1i*unifrnd(0,2*pi,Nr,Nrrf));%init_fully-connected
    init_WRF_pc(:,:,k) = zeros(Nr,Nrrf);
    for i = 1 : Nrrf %initialize W_RF
        init_WRF_pc((i-1) * n + 1:i*n,i,k) = exp( 1i*unifrnd(0,2*pi,n,1)); %init_partially-connected
    end
    init_WRF_ac(:,:,k) = init_WRF_fc(:,:,k);% .* Lw;
end

end
