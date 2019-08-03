function x=NPG(L,rho,delta,c,f,g_res,gx,bx)
M=5;
x=bx;
iter=0;
alpha=1;
gamma=1e-4;
F_nm=zeros(M,1);
penal=max(0,norm(x,1)-c-f-gx'*(x-bx));
grad=g_res+L*(x-bx);
hx=min(grad+2*rho*penal*(1-gx),max(x,grad+2*rho*penal*(-1-gx)));
while 1
    res=norm(hx);
    if res<delta%|iter>1000
        %res
        break
    end
    F=g_res'*(x-bx)+L*norm(x-bx)^2/2+rho*penal^2;%F
    i=rem(iter,M)+1;F_nm(i,1)=F;%max(F_nm)
    while alpha>1e-8 %
        x_n=x-alpha*hx;
        penal=max(0,norm(x_n,1)-c-f-gx'*(x_n-bx));
        F_obj=g_res'*(x_n-bx)+L*norm(x_n-bx)^2/2+rho*penal^2;
        if F_obj<max(F_nm)-gamma*alpha*norm(hx)^2
            break
        end
        alpha=0.2*alpha;
    end
    sx=x_n-x;
    x=x_n;
    sy=hx;
    grad=g_res+L*(x_n-bx);
    hx=min(grad+2*rho*penal*(1-gx),max(x,grad+2*rho*penal*(-1-gx)));
    
    bk=(hx-sy)'*sx;
    ak=sx'*sx;
    if bk<=0
        alpha=1e+4;
    else
        alpha=min(1e+4,max(1e-4,ak/bk));
    end
    iter=iter+1;
    if iter>20000
        fprintf('more than 10000 nonmonotone gradient-type steps.\n');
        break
    end
end