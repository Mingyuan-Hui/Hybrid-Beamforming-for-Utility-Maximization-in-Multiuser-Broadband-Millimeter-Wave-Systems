function [F_RF] = ciw_frf_mo_algorithm(F_RF,v2 ,T, He2,if_fc,Nt ,Ntrf,Nk)

% global Nt Nrf
manifold = complexcirclefactory(Nt*Ntrf);
m = Nt/Ntrf;

switch(if_fc)
    case 1
        L = ones(Nt,Ntrf);
    case 2
        L = zeros(Nt,Ntrf);
        for i = 1 : Ntrf 
            L((i-1) * m + 1:i*m,i) = ones(m,1);
        end
    case 3
%         L = zeros(Nt,Ntrf);
%         for i = 1 : Ntrf-1 
%             L((i-1) * m + 1:i*m,i) = ones(m,1);
%             L((i-1) * m + 1:i*m,i+1) = ones(m,1);
%         end
%         L((Ntrf-1) * m + 1:Ntrf*m,Ntrf) = ones(m,1);
%         L((Ntrf-1) * m + 1:Ntrf*m,1) = ones(m,1);
        di = 2; % division index
        Ltmp = ones(Nt/di,Ntrf/di);
        L = kron(Ltmp,eye(di));

end

problem.M = manifold;

problem.cost = @(x)frf_cost(x,v2 ,T, He2,Ntrf ,Nt, Nk);
problem.egrad = @(x)frf_egrad(x,v2 ,T, He2,1,Ntrf, Nt, Nk, L);




% if if_fc == 1
%     L = ones(Nt,Ntrf);
% elseif if_fc == 2
%     L = zeros(Nt,Ntrf);
%     for i = 1 : Nrf %initialize V_RF
%         L((i-1) * m + 1:i*m,i) = ones(m,1);
%     end
% end

L = L(:);
[x, cost, iter] = conjugategradient(problem,F_RF(:),L);

F_RF = reshape(x,Nt,Ntrf);