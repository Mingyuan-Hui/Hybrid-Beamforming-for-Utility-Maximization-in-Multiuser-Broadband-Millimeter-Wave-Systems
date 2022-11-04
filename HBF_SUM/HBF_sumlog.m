function [Rate,rateresp] = HBF_sumlog(H,h,F_RF,if_fc,Nk, Nt, Nr, Vn, Ns, Ntrf, K, Iter, W_RF,Nrrf,Pt,B,A,W)

v = zeros(1,Nk);
xi = zeros(1,Nk);
H2 =zeros(Nt,K*Ns,Nk);
Fu = zeros(Ntrf,K*Ns,Nk);
% T = zeros(K*Ns,K*Ns,Nk);

Ekl= zeros(Ns,Ns,K,Nk);
Ek=zeros(Nk*Ns,Nk*Ns,K);
Tk=zeros(Nk*Ns,Nk*Ns,K);
%% initialize FD W_D gamma
FD = zeros(Ntrf,K*Ns,Nk);
WD = zeros(K*Nrrf,K*Ns,Nk); 
R = zeros(Nr,Nr,K,Nk);
G = zeros(Nr,Ns,K,Nk);
S = zeros(Nr,Nr,K,Nk);
%% gamma ini
T = W;
%% Fd ini
for i = 1:Nk
    Weff = A(:,:,i);
    v(i) =  Vn * trace(T(:,:,i) * Weff' * Weff) / Pt;  
    H2(:,:,i) = H(:,:,i)'*Weff;  
end
%% Frf ini
F_RF = ciw_frf_mo_algorithm(F_RF,v ,T, H2,if_fc,Nt,Ntrf,Nk);

for i = 1:Nk
    Fu(:,:,i) = ((F_RF')* H2(:,:,i) * T(:,:,i) * (H2(:,:,i)')* F_RF +  v(i) * (F_RF)'*F_RF)^(-1)*F_RF'*H2(:,:,i)*T(:,:,i);
    xi(i) = sqrt(Pt) * (norm(F_RF*Fu(:,:,i),'fro')^2)^(-0.5);
    FD(:,:,i) = Fu(:,:,i)* xi(i);
    F(:,:,i)=B(:,:,i);
    for m1 = 1:K


        R(:,:,m1,i) = Vn * xi(i)^(-2)*eye(Nr)+  xi(i)^(-2)*h(:,:,m1,i)*F(:,:,i)*F(:,:,i)'*h(:,:,m1,i)';
        G(:,:,m1,i) = xi(i)^(-1) * h(:,:,m1,i) * F(:,(m1-1)*Ns+1:m1*Ns,i);
        S(:,:,m1,i) = G(:,:,m1,i)*T((m1-1)*Ns+1:m1*Ns,(m1-1)*Ns+1:m1*Ns,i)*G(:,:,m1,i)';

    end
end
%% Wrf ini
for m3 = 1:K

    W_RF(:,:,m3) = ciw_wrf_mo_algorithm(R(:,:,m3,:),S(:,:,m3,:),W_RF(:,:,m3),T((m3-1)*Ns+1:m3*Ns,(m3-1)*Ns+1:m3*Ns,:),if_fc,Nr,Nrrf,Nk);
end
WRF = W_RF(:,:,1);
for k =2 :K
    WRF=blkdiag(WRF, W_RF(:,:,k));
end

%% Wd ini
for i = 1: Nk

    for m1 = 1:K
        WD((m1-1)*Nrrf+1:m1*Nrrf,(m1-1)*Ns+1:m1*Ns,i) = ...
            xi(i) * ( Vn * W_RF(:,:,m1)'*W_RF(:,:,m1)+ W_RF(:,:,m1)'*h(:,:,m1,i)*F(:,:,i)*F(:,:,i)'*h(:,:,m1,i)'*W_RF(:,:,m1))^(-1)...
            *W_RF(:,:,m1)'*h(:,:,m1,i) * F(:,(m1-1)*Ns+1:m1*Ns,i);

    end

end

%% begin iteration
for idx = 1:Iter
    for i = 1:Nk
        Weff =WRF * WD(:,:,i);
        v(i) =  Vn * trace(T(:,:,i) * Weff' * Weff) / Pt; 
        H2(:,:,i) = H(:,:,i)'*Weff;  
    end
    %% Frf
    F_RF = ciw_frf_mo_algorithm(F_RF,v ,T, H2,if_fc,Nt,Ntrf,Nk);
    %% Fd
    for i = 1:Nk
        Fu(:,:,i) = ((F_RF')* H2(:,:,i) * T(:,:,i) * (H2(:,:,i)')* F_RF +  v(i) * (F_RF)'*F_RF)^(-1)*F_RF'*H2(:,:,i)*T(:,:,i);
        xi(i) = sqrt(Pt) * (norm(F_RF*Fu(:,:,i),'fro')^2)^(-0.5);
        FD(:,:,i) = Fu(:,:,i)* xi(i);
        F(:,:,i)=F_RF*FD(:,:,i);
        for m1 = 1:K


            R(:,:,m1,i) = Vn * xi(i)^(-2)*eye(Nr)+  xi(i)^(-2)*h(:,:,m1,i)*F(:,:,i)*F(:,:,i)'*h(:,:,m1,i)';
            G(:,:,m1,i) = xi(i)^(-1) * h(:,:,m1,i) * F(:,(m1-1)*Ns+1:m1*Ns,i);
            S(:,:,m1,i) = G(:,:,m1,i)*T((m1-1)*Ns+1:m1*Ns,(m1-1)*Ns+1:m1*Ns,i)*G(:,:,m1,i)';

        end
    end
    %% Wrf
    for m3 = 1:K

        W_RF(:,:,m3) = ciw_wrf_mo_algorithm(R(:,:,m3,:),S(:,:,m3,:),W_RF(:,:,m3),T((m3-1)*Ns+1:m3*Ns,(m3-1)*Ns+1:m3*Ns,:),if_fc,Nr,Nrrf,Nk);
    end
    WRF = W_RF(:,:,1);
    for k =2 :K
        WRF=blkdiag(WRF, W_RF(:,:,k));
    end

    %% Wd
    for i = 1: Nk

        for m1 = 1:K
            WD((m1-1)*Nrrf+1:m1*Nrrf,(m1-1)*Ns+1:m1*Ns,i) = ...
                xi(i) * ( Vn * W_RF(:,:,m1)'*W_RF(:,:,m1)+ W_RF(:,:,m1)'*h(:,:,m1,i)*F(:,:,i)*F(:,:,i)'*h(:,:,m1,i)'*W_RF(:,:,m1))^(-1)...
                *W_RF(:,:,m1)'*h(:,:,m1,i) * F(:,(m1-1)*Ns+1:m1*Ns,i);

        end

    end
    %% Ekl
    for i =1 :Nk
        for m1 = 1:K
            TT = Vn * xi(i)^(-2)*W_RF(:,:,m1)'*W_RF(:,:,m1) +  xi(i)^(-2)*W_RF(:,:,m1)'* h(:,:,m1,i)*F(:,:,i)*F(:,:,i)'*h(:,:,m1,i)' * W_RF(:,:,m1) - xi(i)^(-2)*W_RF(:,:,m1)'* h(:,:,m1,i)*F(:,(m1-1)*Ns+1:m1*Ns,i)*F(:,(m1-1)*Ns+1:m1*Ns,i)'*h(:,:,m1,i)' * W_RF(:,:,m1);
            Ekl(:,:,m1,i) = (eye(Ns)+ G(:,:,m1,i)'* W_RF(:,:,m1) * TT^(-1) * W_RF(:,:,m1)' * G(:,:,m1,i) )^(-1);%hef'*Qk^(-1)*hef)^(-1);
        end
    end

    %% Gamma k
    for m2 =1:K
        for ii=1:Nk
            Ek((ii-1)*Ns+1:ii*Ns,(ii-1)*Ns+1:ii*Ns,m2)=Ekl(:,:,m2,ii);
        end
        Tk(:,:,m2)=(-1/(log(det(Ek(:,:,m2))))*(Ek(:,:,m2)^(-1)));
    end

    for jj=1:Nk
        for m3=1:K
            T((m3-1)*Ns+1:m3*Ns,(m3-1)*Ns+1:m3*Ns,jj)=Tk((jj-1)*Ns+1:jj*Ns,(jj-1)*Ns+1:jj*Ns,m3);
        end
    end

    for ii = 1:Nk
        for uu1 = 1:K
            rate(ii,uu1,idx)=log2(det(   Ekl(:,:,uu1,ii)^(-1)  ));
        end
    end



end



rate=squeeze(rate(:,:,Iter));
rateresp=squeeze(real(sum(rate,1))/Nk);
Rate = sum(rateresp)/K;


