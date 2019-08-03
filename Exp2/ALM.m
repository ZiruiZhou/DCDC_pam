function x = ALM(A,b,x,s,c)
%solve min F_rho(x,y)=f(x,y)+rho([g(x,y)]_+^2+[h(x,y)]_+^2);
%delta = 1e-3;

rho = 10;
eta = 1e-2;
lambda = 0;
alpha = 1;

while 1
    x_new = ALM_newiter(A,b,rho,s,c,x,lambda,eta);
    if norm(x_new-x)<1e-2
        break
    end
    x = x_new;
    lambda = max(0, lambda+rho*(norm(x,1)-sum(max(0,x-s))-sum(max(0,-x-s))-c));
    rho = max(2*rho,norm(lambda)^(1+alpha));    
    eta = max(0.1*eta,1e-6);
end