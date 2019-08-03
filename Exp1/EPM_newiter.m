function x = EPM_newiter(Q,q,A,a,ac,B,b,bc,rho,eps,x0,estart)

x = x0;
n = length(x);
A1 = A(:,:,1); a1 = a(:,1); ac1 = ac(1);
A2 = A(:,:,2); a2 = a(:,2); ac2 = ac(2);
B11 = B(:,:,1,1); b11 = b(:,1,1); bc11 = bc(1,1);
B12 = B(:,:,2,1); b12 = b(:,2,1); bc12 = bc(2,1);
B21 = B(:,:,1,2); b21 = b(:,1,2); bc21 = bc(1,2);
B22 = B(:,:,2,2); b22 = b(:,2,2); bc22 = bc(2,2);


for k=1:10000

%     fprintf('EPM: Inner iteration = %d\n',k);

%     f1 = x'*A1*x + x'*a1 + ac1;
    g11 = x'*B11*x + x'*b11 + bc11;
    g12 = x'*B12*x + x'*b12 + bc12;
%     g1 = max(g11,g12);
%     f2 = x'*A2*x + x'*a2 + ac2;
    g21 = x'*B21*x + x'*b21 + bc21;
    g22 = x'*B22*x + x'*b22 + bc22;
%     g2 = max(g21,g22);
%     penfunv = x'*Q*x + q'*x + rho*max([0,f1-g1,f2-g2]);
    
    if abs(g11-g12)<eps
        j1 = [1,2];
    elseif g11>g12
        j1 = 1;
    elseif g11<g12
        j1 = 2;
    end
    nj1 = length(j1);
    B1 = zeros(n,n,nj1);
    b1 = zeros(n,nj1);
    for i=1:nj1
        B1(:,:,i) = B(:,:,j1(i),1);
        b1(:,i) = b(:,j1(i),1);
    end
    
    if abs(g21-g22)<eps
        j2 = [1,2];
    elseif g21>g22
        j2 = 1;
    elseif g21<g22
        j2 = 2;
    end
    nj2 = length(j2);
    B2 = zeros(n,n,nj2);
    b2 = zeros(n,nj2);
    for i=1:nj2
        B2(:,:,i) = B(:,:,j2(i),2);
        b2(:,i) = b(:,j2(i),2);
    end
                
    funv = Inf;
    for J1=1:nj1
        for J2=1:nj2
            BB1 = B1(:,:,J1);
            bb1 = b1(:,J1);
            BB2 = B2(:,:,J2);
            bb2 = b2(:,J2);
            
            grad1 = 2*BB1*x+bb1;
            grad2 = 2*BB2*x+bb2;

%             cvx_solver Gurobi
            cvx_begin quiet
%             cvx_precision low
            variables z(n) t;
            minimize(z'*Q*z+z'*q + rho*t - rho*grad1'*z - rho*grad2'*z);
            subject to
                z'*B11*z+z'*b11+bc11  +  z'*B21*z+z'*b21+bc21 <= t;
                z'*B12*z+z'*b12+bc12  +  z'*B21*z+z'*b21+bc21 <= t;
                z'*B11*z+z'*b11+bc11  +  z'*B22*z+z'*b22+bc22 <= t;
                z'*B12*z+z'*b12+bc12  +  z'*B22*z+z'*b22+bc22 <= t;
                z'*A1*z+z'*a1+ac1  +  z'*B21*z+z'*b21+bc21 <= t;
                z'*A1*z+z'*a1+ac1  +  z'*B22*z+z'*b22+bc22 <= t;
                z'*A2*z+z'*a2+ac2  +  z'*B11*z+z'*b11+bc11 <= t;
                z'*A2*z+z'*a2+ac2  +  z'*B12*z+z'*b12+bc12 <= t;
            cvx_end
            
            f1 = z'*A1*z + z'*a1 + ac1;
            g11 = z'*B11*z + z'*b11 + bc11;
            g12 = z'*B12*z + z'*b12 + bc12;
            g1 = max(g11,g12);
            f2 = z'*A2*z + z'*a2 + ac2;
            g21 = z'*B21*z + z'*b21 + bc21;
            g22 = z'*B22*z + z'*b22 + bc22;
            g2 = max(g21,g22);
            penfunv = z'*Q*z + q'*z + rho*max([0,f1-g1,f2-g2]);
            if penfunv<funv
                funv = penfunv;
                x_new = z;
            end
        end
    end
    
    if norm(x_new - x)<1e-4
        x = x_new;
        break
    else
        x = x_new;
    end
    
    if toc(estart)>7200
        break
    end
    
end

end

