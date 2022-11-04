function egrad = wrf_egrad(R,S,x,if_fc,Nrrf, Nr, Nk, L)

% global  Nrf Nt Nk;
W = reshape(x,Nr,Nrrf);
if if_fc == 1
%     WW = (W'*W)^(-1);
    for i = 1:Nk
%         A = (1/Vn(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+(O(:,:,i))^(-1))^(-2);
%         B = H(:,:,i)'*W*WW;
%         C = B';
%         M = 1/Vn(i)*W*C*A*B;
%         N = 1/Vn(i)*H(:,:,i)*A*B;
%         egrad(:,:,i) = M-N;
%        egrad(:,:,i)=R(:,:,i)*W*W'*S(:,:,i)*W*(W'*R(:,:,i)*W)^(-2)-S(:,:,i)*W*(W'*R(:,:,i)*W)^(-1);
       egrad(:,:,i)=(R(:,:,i)*W*(W'*R(:,:,i)*W)^(-1)*W'-eye(Nr))*S(:,:,i)*W*(W'*R(:,:,i)*W)^(-1);

    end
else
% %     WW = (W'*W)^(-1);
%     for i = 1:Nk
% %         A = (1/Vn(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+(O(:,:,i))^(-1))^(-2);
% %         egrad(:,:,i) = -Nrf/Nt*1/Vn(i)*H(:,:,i)*A*H(:,:,i)'*W;
%     end
end
egrad = sum(egrad,3) .* L ;
egrad = egrad(:);