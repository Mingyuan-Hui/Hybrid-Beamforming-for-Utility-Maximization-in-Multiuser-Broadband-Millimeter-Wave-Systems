function [PEAltminBF]=cal_PEaltmin(h,init_FRF_fc,Fopt ,Ntrf,Nk,K,Ns,init_WRF_fc,Wopt,Nrrf,Nr,Nt,Pt,Vn)
WRFk = zeros(Nr,Nrrf,K);
Fkl = zeros(Nt,Ns,Nk);
F = zeros(Nt,K*Ns,Nk);
phi = zeros(Nrrf,Nrrf,K,Nk);
[ FRFM, FBBM ] = PE_AltMin(init_FRF_fc, Fopt ,Ntrf);
for i = 1:Nk
    FBBM(:,:,i) = sqrt(Pt) * FBBM(:,:,i) / norm(FRFM * FBBM(:,:,i),'fro');
end

for k =1:K
    init_WRF_fc_combined((k-1)*Nr+1:k*Nr,(k-1)*Nrrf+1:k*Nrrf)=init_WRF_fc(:,:,k);
end
[WRF,~]=PE_AltMin(init_WRF_fc_combined,Wopt,K*Nrrf);
for k=1:K
    WRFk(:,:,k)=WRF((k-1)*Nr+1:k*Nr,(k-1)*Nrrf+1:k*Nrrf);
end

for i = 1:Nk
    for k = 1:K
        Fkl(:,:,i)=FRFM*FBBM(:,(k-1)*Ns+1:k*Ns,i);
        F(:,:,i)=FRFM*FBBM(:,:,i);
        phi(:,:,k,i) = Vn *WRFk(:,:,k)' * WRFk(:,:,k) + WRFk(:,:,k)'*h(:,:,k,i) *(F(:,:,i)*F(:,:,i)'-Fkl(:,:,i)*Fkl(:,:,i)')*h(:,:,k,i)'*WRFk(:,:,k);
        rate(i,k)=log2(det( eye(Ns)+phi(:,:,k,i)^(-1)*WRFk(:,:,k)'*h(:,:,k,i)* Fkl(:,:,i)*Fkl(:,:,i)'*h(:,:,k,i)'*WRFk(:,:,k)));
    end
end
PEAltminBF=(sum(squeeze(sum(rate))))/Nk;