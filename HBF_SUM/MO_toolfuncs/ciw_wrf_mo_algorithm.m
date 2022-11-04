function [W_RF] = ciw_wrf_mo_algorithm(R,S,W_RF,T,if_fc,Nr,Nrrf,Nk)
R=squeeze(R);
S=squeeze(S);
% global Nr Nrf
manifold = complexcirclefactory(Nr*Nrrf);
m = Nr/Nrrf;


switch(if_fc)
    case 1
        L = ones(Nr,Nrrf);
    case 2
%         L = zeros(Nr,Nrrf);
%         for i = 1 : Nrrf %initialize V_RF
%             L((i-1) * m + 1:i*m,i) = ones(m,1);
%         end
        L = ones(Nr,Nrrf);
    case 3
%         L = zeros(Nr,Nrrf);
%         for i = 1 : Nrrf-1 %initialize V_RF
%             L((i-1) * m + 1:i*m,i) = ones(m,1);
%             L((i-1) * m + 1:i*m,i+1) = ones(m,1);
%         end
%         L((Nrrf-1) * m + 1:Nrrf*m,Nrrf) = ones(m,1);
%         L((Nrrf-1) * m + 1:Nrrf*m,1) = ones(m,1);

        L = ones(Nr,Nrrf);
%         Lwtmp = ones(Nr/2,Nrrf/2);
%         L = kron(Lwtmp,[1 0; 0 1]);
    case 4
        L = ones(Nr,Nrrf);
            case 5
        L = ones(Nr,Nrrf);
end


problem.M = manifold;

problem.cost = @(x)wrf_cost(R,S,x,T,Nrrf,Nr, Nk);
problem.egrad = @(x)wrf_egrad(R,S,x,1,Nrrf, Nr, Nk, L);

% if if_fc == 1
%     L = ones(Nr,Nrrf);
% else
%     L = zeros(Nr,Nrrf);
%     for i = 1 : Nrf %initialize V_RF
%         L((i-1) * m + 1:i*m,i) = ones(m,1);
%     end
% end





L = L(:);
[x, cost, iter] = conjugategradient(problem,W_RF(:),L);

W_RF = reshape(x,Nr,Nrrf);