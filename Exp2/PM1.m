function x = PM1(A,b,x,s,c)

rho = 10;
eta = 0.1;

while 1      
    
    x_new = PM1_newiter(A,b,rho,s,c,x,eta);
    if norm(x_new-x)<1e-2
        break
    end
    x = x_new;
    rho = 2*rho;
    eta = 0.5*eta;
end

