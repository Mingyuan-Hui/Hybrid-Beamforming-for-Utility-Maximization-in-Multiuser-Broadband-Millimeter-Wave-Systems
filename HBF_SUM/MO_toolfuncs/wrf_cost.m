function cost = wrf_cost(R,S,x,T,Nrrf,Nr, Nk)

% global Ntrf Nt Nk ;


x = reshape(x,Nr,Nrrf);
for i = 1:Nk
%     cost(i) = trace((He2(:,:,i)'*x*(x'*x)^(-1)*x'*He2(:,:,i)/v2(i)+(T(:,:,i))^(-1))^(-1));
    cost(i) = trace(T(:,:,i)-(x'*R(:,:,i)*x)^(-1)*(x'*S(:,:,i)*x));
end
cost = sum(cost);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end