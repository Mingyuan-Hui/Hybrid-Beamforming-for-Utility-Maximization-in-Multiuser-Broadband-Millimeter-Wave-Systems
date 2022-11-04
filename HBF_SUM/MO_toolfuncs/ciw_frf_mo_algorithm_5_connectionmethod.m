function [F_RF] = ciw_frf_mo_algorithm_5_connectionmethod(F_RF,v2 ,T, He2,if_fc,Nt ,Ntrf,Nk)

% global Nt Nrf
manifold = complexcirclefactory(Nt*Ntrf);
m = Nt/Ntrf;

switch(if_fc)
    case 1
        di = 16; % division index
        Ltmp = ones(Nt/di,Ntrf/di);
        L = kron(Ltmp,eye(di));
    case 2
        di = 8; % division index
        Ltmp = ones(Nt/di,Ntrf/di);
        L = kron(Ltmp,eye(di));
    case 3
        %         L = zeros(Nt,Ntrf);
        %         for i = 1 : Ntrf
        %             L((i-1) * m + 1:i*m,i) = ones(m,1);
        %         end
        di = 4; % division index
        Ltmp = ones(Nt/di,Ntrf/di);
        L = kron(Ltmp,eye(di));
    case 4
        di = 2; % division index
        Ltmp = ones(Nt/di,Ntrf/di);
        L = kron(Ltmp,eye(di));
    case 5
        L = ones(Nt,Ntrf);

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