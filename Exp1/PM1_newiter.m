function x = PM1_newiter(Q,q,A,a,ac,B,b,bc,rho,eta,eps,x0)

x = x0;
n = length(x);
A1 = A(:,:,1); a1 = a(:,1); ac1 = ac(1);
A2 = A(:,:,2); a2 = a(:,2); ac2 = ac(2);
B11 = B(:,:,1,1); b11 = b(:,1,1); bc11 = bc(1,1);
B12 = B(:,:,2,1); b12 = b(:,2,1); bc12 = bc(2,1);
B21 = B(:,:,1,2); b21 = b(:,1,2); bc21 = bc(1,2);
B22 = B(:,:,2,2); b22 = b(:,2,2); bc22 = bc(2,2);


for k=1:10000
    
%     fprintf('EPM1: Inner iteration = %d\n',k);

%     f1 = x'*A1*x + x'*a1 + ac1;
    g11 = x'*B11*x + x'*b11 + bc11;
    g12 = x'*B12*x + x'*b12 + bc12;
    g1 = max(g11,g12);
%     f2 = x'*A2*x + x'*a2 + ac2;
    g21 = x'*B21*x + x'*b21 + bc21;
    g22 = x'*B22*x + x'*b22 + bc22;
    g2 = max(g21,g22);
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
            
            cvx_begin quiet
            variables z(n) t1 t2;
            minimize(z'*Q*z+z'*q + rho*t1 + rho*t2)
            subject to
                t1 >= 0;
                z'*A1*z+z'*a1+ac1 - g1 - grad1'*(z-x) <= t1;
                t2 >= 0;
                z'*A2*z+z'*a2+ac2 - g2 - grad2'*(z-x) <= t2;
            cvx_end
            if z'*Q*z+z'*q + rho*t1 + rho*t2<funv
                funv = z'*Q*z+z'*q + rho*t1 + rho*t2;
                x_new = z;
            end
            
        end
    end
    
    if norm(x_new - x)<eta
        x = x_new;
        break
    else
        x = x_new;
    end
    
end

end

