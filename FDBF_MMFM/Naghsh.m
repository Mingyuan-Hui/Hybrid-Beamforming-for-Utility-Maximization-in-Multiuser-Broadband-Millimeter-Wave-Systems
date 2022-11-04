function [Rout,Routresp] = Naghsh(h,Nr,Nt, Ns, Nk, K, Vn,Iter,Pt)
%% 
% simulation codes for M. M. Naghsh, M. Masjedi, A. Adibi, and P. Stoica, "Max–min fairness design for MIMO interference channels: A minorization–maximization approach," IEEE Trans. Signal Process., vol. 67, no. 18, pp. 4707–4719, Jul. 2019.
% Author: Hui Mingyuan
% in this paper
% Mi → Nt
% Lj → Nr
% di → Ns
% Vi → F precoder in Nt × Ns
% H in Nr × Nt
% N → K
for f = 1:Nk
H=h(:,:,:,f);
Ubar = [eye(Ns) zeros(Ns,Nr)].';


%% ini V
Vini = rand(Nt*Ns,K).*exp(1j*2*pi*rand(Nt*Ns,K));
Vini2 = sqrt(Pt)*Vini/norm(Vini,'fro');
for i = 1:K
V(:,:,i)= reshape(Vini2(:,i),Nt,Ns);
% V(:,:,i)=B(:,Ns*(i-1)+1:Ns*i,1);
end




flag = 1;
ite = 1;
    for i =1:K
        NoiseEquiini(:,:,i) = Vn*eye(Nr);
        for jj = setdiff(1:K,i)
            NoiseEquiini(:,:,i) = NoiseEquiini(:,:,i) + H(:,:,i)*V(:,:,jj)*V(:,:,jj)'*H(:,:,i)';
        end
    end
for i = 1:K
    Rini(i,ite) =abs( log2(det(eye(Ns)+V(:,:,i)'*H(:,:,i)'*(NoiseEquiini(:,:,i))^(-1)*H(:,:,i)*V(:,:,i))));
end
R(ite)=min(Rini);
%% in loop
while flag
    ite=ite+1;
    for i =1:K
        Bbar22(:,:,i) = Vn*eye(Nr);
        for jj = 1:K
            Bbar22(:,:,i) = Bbar22(:,:,i) + H(:,:,i)*V(:,:,jj)*V(:,:,jj)'*H(:,:,i)';
        end
    end
    for i = 1:K
       Bbar(:,:,i)=[eye(Ns) V(:,:,i)'*H(:,:,i)'; H(:,:,i)*V(:,:,i) Bbar22(:,:,i)];
        Bbarinv = inv(Bbar(:,:,i));
        Fbar(:,:,i) = Bbarinv*Ubar*(Ubar'*Bbarinv*Ubar)^(-1)*Ubar'*Bbarinv;
    end

    for i = 1:K
        Fbar11(:,:,i)=Fbar(1:Ns,1:Ns,i);
        Fbar12(:,:,i)=Fbar(1:Ns,Ns+1:end,i);
        Fbar21(:,:,i)=Fbar(Ns+1:end,1:Ns,i);
        Fbar22(:,:,i)=Fbar(Ns+1:end,Ns+1:end,i);
    end

    for i = 1:K
        Gbar(:,:,i) = kron(eye(Ns),H(:,:,i)'*Fbar22(:,:,i)*H(:,:,i));
        Cbar(i) = - log(det(Ubar'*Bbar(:,:,i)^(-1)*Ubar))-trace(Fbar(:,:,i)*Bbar(:,:,i))+trace(Fbar11(:,:,i))+trace(Fbar22(:,:,i)*Vn*eye(Nr));
        bbar(:,:,i) = vec(H(:,:,i)'*Fbar12(:,:,i)');
    end
    



    cvx_begin quiet
%     cvx_begin 
    variable t 
    variable xbar(Nt*Ns,K) 
    expressions temp1(K,1) 
    for i = 1:K
    for jj = 1:K
        temp1(i) = temp1(i) + real(quad_form( xbar(:,jj), Gbar(:,:,i)' ));   
    end
        
    end

    maximize( t )
    subject to
    for i = 1:K
        real(Cbar(i)) + 2 * real(bbar(:,:,i)'*xbar(:,i)) + temp1(i) <= -t
    end
    norm(xbar,'fro') <= sqrt(Pt)
    cvx_end
    for i =1:K
        V(:,:,i)=reshape(xbar(:,i),Nt,Ns);
    end

    for i =1:K
        NoiseEqui(:,:,i) = Vn*eye(Nr);
        for jj = setdiff(1:K,i)
            NoiseEqui(:,:,i) = NoiseEqui(:,:,i) + H(:,:,i)*V(:,:,jj)*V(:,:,jj)'*H(:,:,i)';
        end
    end
    for i = 1:K
        Rite(i,ite) = abs(log2(det(eye(Ns)+V(:,:,i)'*H(:,:,i)'*(NoiseEqui(:,:,i))^(-1)*H(:,:,i)*V(:,:,i))));
    end
    R(ite) = min(Rite(:,ite));

        fprintf('In subcarrier %d, current ite is %d. The current min rate is %f.\n', f,ite-1,R(ite) )
    if  abs(1-abs(R(ite-1))/abs(R(ite)))<0.001
        flag = 0;
    end
end

Routall(f)=R(ite);
Routrespall(f,:) = Rite(:,ite);
end
Routresp=squeeze(real(sum(Routrespall,1))/Nk);
Rout = sum(Routresp)/K;





