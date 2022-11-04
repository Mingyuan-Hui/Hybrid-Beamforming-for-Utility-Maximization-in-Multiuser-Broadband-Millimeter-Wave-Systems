%% Main_SRM_comparison
% Main program for comparing the performance of the proposed HBF-SRM algorithm with that of the FDBF one and some exist algorithms. (Fig.4)
% Author: Hui Mingyuan
% Date: 20221103
clear
close all
disp(datestr(now))
tic

%% Initialization parameters
Nt = 128;
Ntrf = 16;
Nr = 8;
Nrrf = 2;
Ns = 2;
K = 8;
Nk = 64;
Iter = 10;
numMC = 8;
if_fc =1 ;

%% pathloss & power of noise
PON=-156+10*log10(100*1e6/Nk); %power of background noise  dbm
pd=gen_lossmultiple(K);
lossmultiple=pd;
%% POT and channel
POT=0:5:20; % PowerOfTransmit dbm

h=zeros(Nr,Nt,K,Nk,numMC);

for ii = 1:K
    for kk = 1:numMC
        [h(:,:,ii,:,kk),~,~] = OMPHWB(Nt,Nr,Nk);
    end
    h(:,:,ii,:,:)=sqrt(lossmultiple(ii)/mean(pd)).*h(:,:,ii,:,:);
end

H = zeros(K*Nr,Nt,Nk,numMC);
for mc = 1:numMC
    for k = 1:K
        H(Nr*(k-1)+1:Nr*k,:,:,mc)=h(:,:,k,:,mc);
    end
end


%% 
POTnum=length(POT);
Vn =db2pow(PON)/mean(pd); % ¦Ò^2
parfor n = 1 : numMC
    [init_FRF_pc, init_FRF_fc,init_WRF_fc,init_WRF_pc,init_FRF_ac,init_WRF_ac]=generatefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf);
    for POT_idx = 1:POTnum
        Pt = db2pow(POT(POT_idx));
        fprintf('In simulation %d, current POT is %d dBm.\n ', n,POT(POT_idx) )
        %%
        [~,~,B,A,W] = FDBF_SRM(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,Iter,Pt);
        %%
        [HBF_typicalsum(n,POT_idx),~] = HBF_SRM(H(:,:,:,n),h(:,:,:,:,n),init_FRF_fc,if_fc,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter,init_WRF_fc,Nrrf,Pt,B,A,W);
        %%
        [FDBF_typicalsum(n,POT_idx),~,Fopt,Wopt,~] = FDBF_SRM(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,2*Iter,Pt); % "*2" is for fairness
        %% 
        PEAltminBF(n,POT_idx)=cal_PEaltmin(h(:,:,:,:,n),init_FRF_fc,Fopt ,Ntrf,Nk,K,Ns,init_WRF_fc,Wopt,Nrrf,Nr,Nt,Pt,Vn);
        %%
        CTDBF(n,POT_idx)=cal_CTDv2(h(:,:,:,:,n),Ns, Vn,Pt,Nk,K)  ;
        %%
        DZTBF(n,POT_idx)=cal_DZTv2(h(:,:,:,:,n),Ns, Vn,Pt,Nk,K)  ;
    end
end
disp(datestr(now))
toc

HBF_typicalsum_mean = sum(HBF_typicalsum)/numMC;
FDBF_typicalsum_mean = sum(FDBF_typicalsum)/numMC;
PEAltminBF_mean = sum(PEAltminBF)/numMC;
CTDBF_mean = sum(CTDBF)/numMC;
DZTBF_mean = sum(DZTBF)/numMC;

%%
figure
plot(POT, FDBF_typicalsum_mean, 'LineWidth', 2);hold on
plot(POT, HBF_typicalsum_mean, 'LineWidth', 2);hold on
plot(POT, CTDBF_mean, 'LineWidth', 2);hold on
plot(POT, DZTBF_mean, 'LineWidth', 2);hold on
plot(POT, PEAltminBF_mean, 'LineWidth', 2);hold on
grid on
legend('FDBF-WMMSE','HBF-SRM','HBF-CTD','HBF-TU','HBF-PE-Altmin')
xlabel('$P_{\mathrm{t}}$(dBm)','Interpreter','latex')
ylabel('Spectral efficiency (bit/s/Hz)','Interpreter','latex')


