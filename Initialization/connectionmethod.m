function [Lf,Lw] = connectionmethod(inputArg1,inputArg2)
%CONNECTIONMETHOD 此处显示有关此函数的摘要
%   此处显示详细说明

%% Frf
m = Nt/Ntrf;
switch(methodindex)
    case 1
        Lf = ones(Nt,Ntrf);
    case 2
        Lf = zeros(Nt,Ntrf);
        for i = 1 : Ntrf 
            Lf((i-1) * m + 1:i*m,i) = ones(m,1);
        end
    case 3
        Lf = zeros(Nt,Ntrf);
        for i = 1 : Ntrf-1 
            Lf((i-1) * m + 1:i*m,i) = ones(m,1);
            Lf((i-1) * m + 1:i*m,i+1) = ones(m,1);
        end
        Lf((Ntrf-1) * m + 1:Ntrf*m,Ntrf) = ones(m,1);
        Lf((Ntrf-1) * m + 1:Ntrf*m,1) = ones(m,1);
end

%% Wrf
b = Nr/Nrrf;
switch(methodindex)
    case 1
        Lw = ones(Nr,Nrrf);
    case 2
        Lw = zeros(Nr,Nrrf);
        for i = 1 : Nrrf %initialize V_RF
            Lw((i-1) * b + 1:i*b,i) = ones(b,1);
        end
    case 3
        Lw = zeros(Nr,Nrrf);
        for i = 1 : Nrrf-1 %initialize V_RF
            Lw((i-1) * b + 1:i*b,i) = ones(b,1);
            Lw((i-1) * b + 1:i*b,i+1) = ones(b,1);
        end
        Lw((Nrrf-1) * b + 1:Nrrf*b,Nrrf) = ones(b,1);
        Lw((Nrrf-1) * b + 1:Nrrf*b,1) = ones(b,1);
end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

