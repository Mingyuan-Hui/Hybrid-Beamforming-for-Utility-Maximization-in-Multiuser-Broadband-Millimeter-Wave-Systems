function [Rate,rateresp,B,A,W]=FDBF_sumlog(H,h,Nr, Nt, Ns, Nk, K, Vn,Iter,Pt)

KNs=K*Ns;
B=zeros(Nt,KNs,Nk);
for jj=1:Nk
    Hv=[];
    for ii=1:K
        [~,~,V]=svd(h(:,:,ii,jj));
        Hv=[Hv; V(:,1:Ns)'];
    end
    Btmp=Hv'*inv(Hv*Hv');
    B(:,:,jj)=sqrt(Pt)*Btmp/norm(Btmp,'fro');

end


A = zeros(K*Nr,K*Ns,Nk);
W = zeros(K*Ns,K*Ns,Nk);

for n = 1:Iter  
    for jj=1:Nk  
        for uu1 = 1:K  
            Qk=Vn*eye(Nr); 
            Bkn=B(:,(uu1-1)*Ns+1:uu1*Ns,jj);  
            Heff=h(:,:,uu1,jj)*Bkn; 
            for uu2 = setdiff(1:K,uu1)  
                Qk = Qk + h(:,:,uu1,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj)'*h(:,:,uu1,jj)'; 
            end 
            A((uu1-1)*Nr+1:uu1*Nr,(uu1-1)*Ns+1:uu1*Ns,jj) = (Qk+Heff*Heff')^(-1)*Heff;

                                Einv= (eye(Ns)+Heff'*Qk^(-1)*Heff);
            E(:,:,uu1,jj)=(eye(Ns)+Heff'*Qk^(-1)*Heff)^(-1);
            
        end

    end

    for uu1 =1:K
        for jj=1:Nk
            Ek((jj-1)*Ns+1:jj*Ns,(jj-1)*Ns+1:jj*Ns,uu1)=E(:,:,uu1,jj);
            Ekinv((jj-1)*Ns+1:jj*Ns,(jj-1)*Ns+1:jj*Ns,uu1)=E(:,:,uu1,jj)^(-1);
        end
     Tk(:,:,uu1)=(-1/(log(det(Ek(:,:,uu1))))*(Ek(:,:,uu1)^(-1))); 
    end
    
    for jj=1:Nk
        for uu1=1:K
            W((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns,jj)=Tk((jj-1)*Ns+1:jj*Ns,(jj-1)*Ns+1:jj*Ns,uu1);
        end
    end



    
    for jj = 1:Nk
        for uu1 = 1:K 
        v =  Vn * trace(W(:,:,jj) * A(:,:,jj)' * A(:,:,jj))/Pt;
        H1 = H(:,:,jj)'*A(:,:,jj);
        B(:,:,jj) = ( H1 * W(:,:,jj) * H1' +  v*eye(Nt))^(-1)*H1*W(:,:,jj);
        B(:,:,jj) =sqrt(Pt)/norm(B(:,:,jj),'fro')* B(:,:,jj); 
        end
    end
end

    for jj = 1:Nk
        for uu1 = 1:K
                        Qk=Vn*eye(Nr); 
            Bkn=B(:,(uu1-1)*Ns+1:uu1*Ns,jj);  
            Heff=h(:,:,uu1,jj)*Bkn; 
            for uu2 = setdiff(1:K,uu1)  
                Qk = Qk + h(:,:,uu1,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj)'*h(:,:,uu1,jj)'; 
            end 
            E(:,:,uu1,jj)=(eye(Ns)+Heff'*Qk^(-1)*Heff)^(-1);
            rate(jj,uu1)=log2(det(   E(:,:,uu1,jj)^(-1)  ));
        end
    end
rateresp=squeeze(real(sum(rate,1))/Nk);
Rate = sum(rateresp)/K;

