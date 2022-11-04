function smy=cal_DZT(h,Ns, sigma2,Pt_true,Nk,K)  
        H=permute(h,[1 2 4 3]);
        Pt = K *Ns;
        rho=(Pt_true/Pt);

%         rho=1;
[ F_my, W_hy] = DZhang2019_TensorUnfolding( H, Ns, sigma2, Pt_true, rho);

%    smy       = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_hy,rho,'total',Pt),W_hy,Ns,rho,sigma2);
% smy       = MUspectralEfficiency(H,F_my,W_my,Ns,rho,sigma2);
% 
%                     HH=permute(h,[1 2 4 3]);
%                     [F_hybrid,W_hy,~,FRF,FD]=myMUOFDMBeamformingRCD(HH,Ns, Vn,Pt, 1); 
%                     F=F_hybrid;



%                     F_my=OFDMWaterFilling(H,F_my,W_hy,rho,'total',Pt);
                    for i = 1:Nk
                        for k = 1:K
                            F(:,(k-1)*Ns+1:k*Ns,i)=F_my(:,:,i,k);

                        end
                        F(:,:,i)= sqrt(Pt_true)*F(:,:,i)/norm(F(:,:,i),'fro');

                        for k = 1:K
                            Fkl(:,:,i) = F(:,(k-1)*Ns+1:k*Ns,i);
                            phi(:,:,k,i) = sigma2 *W_hy(:,:,i,k)' * W_hy(:,:,i,k) + W_hy(:,:,i,k)'*h(:,:,k,i) *(F(:,:,i)*F(:,:,i)'-Fkl(:,:,i)*Fkl(:,:,i)')*h(:,:,k,i)'*W_hy(:,:,i,k);
                            rateresp_allsc(i,k)=log2(det( eye(Ns)+phi(:,:,k,i)^(-1)*W_hy(:,:,i,k)'*h(:,:,k,i)* Fkl(:,:,i)*Fkl(:,:,i)'*h(:,:,k,i)'*W_hy(:,:,i,k)));
                        end
                    end
rateresp = sum(rateresp_allsc,1)/Nk;
smy=sum(rateresp);
end