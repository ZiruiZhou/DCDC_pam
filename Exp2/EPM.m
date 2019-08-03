function x = EPM(A,b,x,s,c)

rho = 10;

while 1      
    
    x_new = EPM_newiter(A,b,rho,s,c,x);
    if norm(x_new-x)<1e-2
        break
    end
    x = x_new;
    rho = 2*rho;
end