function [ FRF,FBB ] = PE_AltMin( FRF, Fopt ,NRF)
%% Reference
% Simulation codes for "Alternating minimization algorithms for hybrid precoding in millimeter wave MIMO systems," by Xianghao Yu, Juei-Chin Shen, Jun Zhang, and Khaled B. Letaief, IEEE J. Sel. Topics Signal Process., 2016.
% https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems
[Nt, Ns, K] = size(Fopt);
mynorm = [];
while (isempty(mynorm) || abs( (mynorm(1) - mynorm(2))/mynorm(1)) > 0.001)
    mynorm = [0,0];
    temp = zeros(Nt, NRF);
    for k = 1:K
        [U,S,V] = svd(Fopt(:,:,k)'*FRF);
        FBB(:,:,k) = V(:,[1:Ns])*U';
        mynorm(1) = mynorm(1) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
        temp = temp + Fopt(:,:,k) * FBB(:,:,k)';
    end

    FRF = exp(1i * angle(temp));
    for k = 1:K
        mynorm(2) = mynorm(2) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
    end
end
end