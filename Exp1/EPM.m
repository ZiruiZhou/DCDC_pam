function x = EPM(Q,q,A,a,ac,B,b,bc,rho,eps,x0,itermax,estart)

x = x0;

for k=1:itermax
    
%     fprintf('EPM: Outer iteration = %d\n',k);

    x_new = EPM_newiter(Q,q,A,a,ac,B,b,bc,rho,eps,x,estart);
    if norm(x_new - x)<1e-2
        x = x_new;
        break
    end
    x = x_new;
    rho = 5*rho;
    if toc(estart)>7200
        break
    end
end

end

