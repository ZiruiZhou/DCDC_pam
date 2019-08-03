function x = ALM(Q,q,A,a,ac,B,b,bc,rho,eps,x0,itermax)

x = x0;
eta = 0.1;
lambda = zeros(2,1);
alpha = 1;

for k=1:itermax
    
%     fprintf('QPM: Outer iteration = %d\n',k);

    x_new = ALM_newiter(Q,q,A,a,ac,B,b,bc,rho,eta,eps,x,lambda);
    if norm(x_new - x)<0.01
        x = x_new;
        break
    end
    x = x_new;
    [~,f1,g11,g12,f2,g21,g22] = funv_eval(Q,q,A,a,ac,B,b,bc,x_new);
    lambda = max(0, lambda+rho*[f1 - max(g11,g12); f2 - max(g21,g22)]);
    rho = max(2*rho,norm(lambda)^(1+alpha));
    eta = max(0.5*eta,1e-6);
end

end


