function cost = frf_cost(x,v2 ,T, He2,Ntrf ,Nt, Nk)



x = reshape(x,Nt,Ntrf);
for i = 1:Nk
    cost(i) = trace((He2(:,:,i)'*x*(x'*x)^(-1)*x'*He2(:,:,i)/v2(i)+(T(:,:,i))^(-1))^(-1));
end
cost = sum(cost);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end