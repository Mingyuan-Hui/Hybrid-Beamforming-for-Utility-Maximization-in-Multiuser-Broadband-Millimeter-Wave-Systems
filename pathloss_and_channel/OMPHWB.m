function [H,AT,AR] = OMPHWB(N_t,N_r,N)
% N=64;
N_c=5;
N_ray=10;
E_aoa = 2*pi* rand(N_c,1);                               
sigma_aoa = 10*pi/180;                                    
b = sigma_aoa/sqrt(2);                                    
a = rand(N_c,N_ray)-0.5;                                  
aoa = repmat(E_aoa,1,N_ray)-b*sign(a).*log(1-2*abs(a));   
aoa = sin(aoa);
%-----------AOD
E_aod = 2*pi* rand(N_c, 1);                               
sigma_aod = 10*pi/180;                                    
b = sigma_aod/sqrt(2);                                    
a = rand(N_c,N_ray)-0.5;                                  
aod = repmat(E_aod,1, N_ray)-b*sign(a).*log(1-2*abs(a));   
aod = sin(aod);


signature_t = [0:(N_t-1)]';
signature_t = 1i*pi* signature_t;                           
signature_r = [0:(N_r-1)]';
signature_r = 1i*pi* signature_r;                           

                                             
H_ray = zeros(N_r, N_t, N_c, N_ray);
H_cl = zeros(N_r, N_t, N_c);
H = zeros(N_r, N_t, N);                                 

    for i= 1: N_c
        for m = 1: N_ray
            H_ray(:,:,i,m)=complex(randn(1),randn(1))/sqrt(2)*exp((aoa(i,m)*signature_r))*exp((aod(i,m)*signature_t))'/sqrt(N_t*N_r); 
        end
    end  
        H_cl = sum(H_ray, 4);    
    
for k = 1: N
    for i = 1: N_c
    H_cltmp(:,:,i) = H_cl(:,:,i)*exp(-1i*2*pi*(i-1)*(k-1)/N); 

    end
    H(:,:,k) = sqrt(N_t*N_r/N_c/N_ray)*sum(H_cltmp,3);       
end


    
    aod = aod(:).';
    aoa = aoa(:).';
    A = kron(aod,signature_t);
    AT = 1/sqrt(N_t)*exp(A);
    A = kron(aoa,signature_r);
    AR = 1/sqrt(N_r)*exp(A);