function egrad = frf_egrad(x,V3 ,O, H,if_fc,Ntrf, Nt, Nk,L)

W = reshape(x,Nt,Ntrf);
if if_fc == 1
    WW = (W'*W)^(-1);
    for i = 1:Nk
        Mexp_m2 = (1/V3(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+(O(:,:,i))^(-1))^(-2);
        B = H(:,:,i)'*W*WW;
        C = B';
        M = 1/V3(i)*W*C*Mexp_m2*B;
        N = 1/V3(i)*H(:,:,i)*Mexp_m2*B;
        egrad(:,:,i) = M-N;
    end
else
    WW = (W'*W)^(-1);
    for i = 1:Nk
        Mexp_m2 = (1/V3(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+(O(:,:,i))^(-1))^(-2);
        egrad(:,:,i) = -Ntrf/Nt*1/V3(i)*H(:,:,i)*Mexp_m2*H(:,:,i)'*W;
    end
end
egrad = sum(egrad,3) .* L ;
egrad = egrad(:);