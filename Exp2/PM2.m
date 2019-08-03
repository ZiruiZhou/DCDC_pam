function x = PM2(A,b,x,s,c)

rho = 10;
eta = 1e-2;

while 1
    x_new = PM2_newiter(A,b,rho,s,c,x,eta);
    if norm(x_new-x)<1e-2
        break
    end
    x = x_new;
    rho = 2*rho;
    eta = max(0.1*eta,1e-6);
end