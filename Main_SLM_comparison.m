%% Main_SLM_comparison
% Main program for comparing the performance of the proposed HBF-SLM algorithm with that of the FDBF one.
% Author: Hui Mingyuan
% Date: 20221103
clear
close all
disp(datestr(now))
tic

%% Initialization parameters
Nt = 32; 
Ntrf = 16;
Nr = 8; 
Nrrf = 2; 
Ns = 2; 
K = 8; 
Nk = 64; 
Iter = 10;
numMC = 8; 
if_fc=1;
%% pathloss & power of noise
PON=-156+10*log10(100*1e6/Nk); %power of background noise  dbm 
pd=gen_lossmultiple(K);
lossmultiple=pd;
%% POT and channel
POT=0:2:20; % PowerOfTransmit (dBm)

h=zeros(Nr,Nt,K,Nk,numMC);

for ii = 1:K
    for kk = 1:numMC
        [h(:,:,ii,:,kk),~,~] = OMPHWB(Nt,Nr,Nk);
    end
    h(:,:,ii,:,:)=sqrt(lossmultiple(ii)).*h(:,:,ii,:,:); 
end

H = zeros(K*Nr,Nt,Nk,numMC);
for mc = 1:numMC
    for k = 1:K
        H(Nr*(k-1)+1:Nr*k,:,:,mc)=h(:,:,k,:,mc);
    end
end

%% 
TotalRateofFDBF_sumlog = zeros(numMC,length(POT));
UserRateofFDBF_sumlog = zeros(numMC,length(POT),K);
TotalRateofHBF_sumlog = zeros(numMC,length(POT));
UserRateofHBF_sumlog = zeros(numMC,length(POT),K);

POTnum=length(POT);
Vn =db2pow(PON); % ¦Ò^2

        
parfor n = 1 : numMC
        [init_FRF_pc, init_FRF_fc,init_WRF_fc,init_WRF_pc,init_FRF_ac,init_WRF_ac]=generatefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf);
    for POT_idx = 1:POTnum
        Pt = db2pow(POT(POT_idx));
                fprintf('In simulation %d, current POT is %d dBm.\n ', n,POT(POT_idx) )   
        %% MO utility hybrid receiver log
        %VFD initialization
        [~,~,B,A,W] = FDBF_sumlog(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,Iter,Pt); %for the initialization of HBF-SUM
        %WMMSE based HBF-SUM
        [TotalRateofHBF_sumlog(n,POT_idx),UserRateofHBF_sumlog(n,POT_idx,:)] = ...
            HBF_sumlog(H(:,:,:,n),h(:,:,:,:,n),init_FRF_fc,if_fc,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter,init_WRF_fc,Nrrf,Pt,B,A,W); 
        %WMMSE based FDBF-SUM 
        [TotalRateofFDBF_sumlog(n,POT_idx),UserRateofFDBF_sumlog(n,POT_idx,:),B,A,W] = FDBF_sumlog(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,2*Iter,Pt);% "*2" is for fairness
    end
end
disp(datestr(now))
toc
%% 
UserRateofFDBF_sumlog_mean = squeeze(sum(UserRateofFDBF_sumlog)/numMC);
UserRateofHBF_sumlog_mean = squeeze(sum(UserRateofHBF_sumlog)/numMC);
%% Utility Comparison log
UtilityofHBF_sumlog_mean=log(prod(UserRateofHBF_sumlog_mean'));
UtilityofFDBF_sumlog=log(prod(UserRateofFDBF_sumlog_mean'));
% Fig. 3 in the revised manuscript.
figure 
plot(POT, UtilityofFDBF_sumlog);hold on
plot(POT, UtilityofHBF_sumlog_mean);
legend('FDBF-SLM','HBF-SLM')
