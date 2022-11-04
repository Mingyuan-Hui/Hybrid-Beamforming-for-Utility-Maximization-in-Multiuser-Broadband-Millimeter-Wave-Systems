%% Main_convergence
% Convergence of the proposed HBF-SUM algorithm (taking HBF-SLM as example)
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
Iter = 15;
numMC = 8;
if_fc=1;

%% pathloss & power of noise
PON=-156+10*log10(100*1e6/Nk); %power of background noise  dbm
pd=gen_lossmultiple(K);
lossmultiple=pd;
%% POT and channel
POT = 12; % PowerOfTransmit dbm
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

for FDhelpflag =0:1

    %%
    RateofHBF_sumlog = zeros(numMC,K,Iter);

    POTnum=length(POT);
    Vn =db2pow(PON); % ¦Ò^2
    parfor n = 1 : numMC
        [init_FRF_pc, init_FRF_fc,init_WRF_fc,init_WRF_pc,init_FRF_ac,init_WRF_ac]=generatefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf);
        [~,~,B,A,W] = FDBF_sumlog(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,Iter,db2pow(POT));
        Pt = db2pow(POT);
        %%
        %If VFD is used
        if FDhelpflag ==1
            [~,RateofHBF_sumlog(n,:,:)] = ...
                HBF_sumlog_for_converge_test(H(:,:,:,n),h(:,:,:,:,n),init_FRF_fc,1,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter,init_WRF_fc,Nrrf,Pt,B,A,W); %fully-connected
            %If VFD is not used
        else
            [~,RateofHBF_sumlog(n,:,:),~,~,~,~]=HBF_sumlog_withoutVFD_for_converge_test(H(:,:,:,n),h(:,:,:,:,n),init_FRF_fc,1,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter,init_WRF_fc,Nrrf,Pt);

        end
        disp(n)
    end
    disp(datestr(now))
    toc

    RateofHBF_sumlog_mean = squeeze(sum(RateofHBF_sumlog)/numMC);

    UtilityofHBF_logsum=sum(log(RateofHBF_sumlog_mean));
    datasave{FDhelpflag+1}=UtilityofHBF_logsum;

end

% fig. 2(a)
figure
plot(1:Iter,datasave{2}); hold on
plot(1:Iter,datasave{1}); hold on
legend('HBF-SLM (VFD-ini)','HBF-SLM (Random-ini)')
