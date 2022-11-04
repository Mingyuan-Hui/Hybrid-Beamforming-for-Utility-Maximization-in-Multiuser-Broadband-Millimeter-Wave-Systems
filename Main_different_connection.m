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
%%

%% pathloss & power of noise
PON=-156+10*log10(100*1e6/Nk); %power of background noise  dbm
pd=gen_lossmultiple(K);
% load('userposition_and_lossmultiple1.mat')
lossmultiple=pd;
%% channel
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
Vn =db2pow(PON); % ¦Ò^2

parfor n = 1 : numMC
    [init_FRF_ac,init_WRF]=generatemorefrfandwrf( Nt, Nr,  Ntrf, K, Nrrf);
    pot=[4 12 20];
    for potidx=1:3
        Pt = db2pow(pot(potidx));
        [~,~,B,A,W] = FDBF_SRM(H(:,:,:,n),h(:,:,:,:,n),Nr,Nt, Ns, Nk, K, Vn,Iter,Pt);
        for psidx = 1:5
            fprintf('In simulation %d, current point is %d, Pot = %d dBm. \n ', n,psidx,pot(potidx) )
            [RateofHBF_typicalsum(n,psidx,potidx),~,~,~,~] = ...
                HBF_SRM_for_connectionmethod_test(H(:,:,:,n),h(:,:,:,:,n),init_FRF_ac(:,:,psidx),psidx,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter,init_WRF,Nrrf,Pt,B,A,W);
        end
    end
end
disp(datestr(now))
toc
RateofHBF_typicalsum_mean = squeeze(sum(RateofHBF_typicalsum)*K/numMC);

% %% sum rate figure
figure
semilogx([1/16 1/8 1/4 1/2 1], RateofHBF_typicalsum_mean(:,1),'color',[12/255,127/255,63/255]);hold on
semilogx([1/16 1/8 1/4 1/2 1], RateofHBF_typicalsum_mean(:,2),'color',[185/255,81/255,159/255]);hold on
semilogx([1/16 1/8 1/4 1/2 1], RateofHBF_typicalsum_mean(:,3),'color',[222/255,127/255,39/255]);hold on
legend('4dBm','12dBm','20dBm')
