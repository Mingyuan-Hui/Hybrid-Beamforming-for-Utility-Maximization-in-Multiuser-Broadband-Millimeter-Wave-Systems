function [pd]=gen_lossmultiple(K)

%% gen_location

% position ini
Pt=zeros(K,2);
Lroom=0;
Wroom=0;
k1=[1,0];
k2=[0,1];
R=300;
% Pts are the location of the UEs
for k0=1:(K/2)
    r=(0.1+0.9*rand(1))*R;
    theta=rand(1,1)*2*pi;
    px=r*cos(theta);
    py=r*sin(theta);
    pt=[Lroom,Wroom]+px*k1+py*k2;
    Pt(k0,:)=pt;
end

for k0=((K/2)+1):K
    r=(0.1+0.9*rand(1))*R;
    theta=rand(1,1)*2*pi;
    px=r*cos(theta);
    py=r*sin(theta);
    pt=[Lroom,Wroom]+px*k1+py*k2;
    Pt(k0,:)=pt;
end
%% gen_pathloss
AP=[0,0];
Ld=zeros(1,K);
for k0=1:K
    pt=Pt(k0,:)-AP;
    dk=sqrt(sum(abs(pt).^2));
    Ld(k0)=path_NLOS( dk );
%     Ld_test(k0)=path_LOS( dk );
end

%% gen_dB-power

% PON=-174+10*log10(100*1e6);

path_d=db2pow(-Ld);

% path_d=10.^((-noise-Ld)/10); %AP-UE （dBm）
% 
% pd=sqrt(path_d);
pd = path_d;

% pd=repmat(pd.',1,Nt); %4*4








end