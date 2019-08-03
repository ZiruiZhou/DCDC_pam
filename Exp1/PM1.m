function x = PM1(Q,q,A,a,ac,B,b,bc,rho,eps,x0,itermax)

x = x0;
eta = 0.1;

for k=1:itermax
    
%     fprintf('EPM1: Outer iteration = %d\n',k);

    x_new = PM1_newiter(Q,q,A,a,ac,B,b,bc,rho,eta,eps,x);
    if norm(x_new - x)<0.01
        x = x_new;
        break
    end
    x = x_new;
    rho = 5*rho;
    eta = max(0.5*eta,1e-6);
end

end

