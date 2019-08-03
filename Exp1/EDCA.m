function x = EDCA(Q,q,A,a,ac,B,b,bc,eps,x0,itermax)

x = x0;
n = length(x);
A1 = A(:,:,1); a1 = a(:,1); ac1 = ac(1);
A2 = A(:,:,2); a2 = a(:,2); ac2 = ac(2);
B11 = B(:,:,1,1); b11 = b(:,1,1); bc11 = bc(1,1);
B12 = B(:,:,2,1); b12 = b(:,2,1); bc12 = bc(2,1);
B21 = B(:,:,1,2); b21 = b(:,1,2); bc21 = bc(1,2);
B22 = B(:,:,2,2); b22 = b(:,2,2); bc22 = bc(2,2);


[~,feas_x] = funv_feas(Q,q,A,a,ac,B,b,bc,x);
while feas_x>0
    x = randn(n,1);
    [~,feas_x] = funv_feas(Q,q,A,a,ac,B,b,bc,x);
end
if feas_x>0
    error('The input of EDCA must be a feasible point.')
end

for k=1:itermax
    
%     fprintf('EDCA: Iteration = %d\n', k);
    
    
    g11 = x'*B11*x + x'*b11 + bc11;
    g12 = x'*B12*x + x'*b12 + bc12;
    g1 = max(g11,g12);
    g21 = x'*B21*x + x'*b21 + bc21;
    g22 = x'*B22*x + x'*b22 + bc22;
    g2 = max(g21,g22);
    
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
            
%             cvx_begin quiet
%             variable z(n);
%             minimize(z'*Q*z+z'*q)
%             subject to
%                 z'*A1*z+z'*a1+ac1 - g1 - grad1'*(z-x) <= 0;
%                 z'*A2*z+z'*a2+ac2 - g2 - grad2'*(z-x) <= 0;
%             cvx_end
%             if z'*Q*z+z'*q<funv
%                 funv = z'*Q*z+z'*q;
%                 x_new = z;
%             end
            
            cvx_begin quiet
            variable z(n);
            minimize(z'*Q*z+z'*q)
            subject to
                z'*A1*z+z'*a1+ac1  +  z'*B21*z+z'*b21+bc21 <= g1 + grad1'*(z-x) + g2 + grad2'*(z-x);
                z'*A1*z+z'*a1+ac1  +  z'*B22*z+z'*b22+bc22 <= g1 + grad1'*(z-x) + g2 + grad2'*(z-x);
                z'*A2*z+z'*a2+ac2  +  z'*B11*z+z'*b11+bc11 <= g1 + grad1'*(z-x) + g2 + grad2'*(z-x);
                z'*A2*z+z'*a2+ac2  +  z'*B12*z+z'*b12+bc12 <= g1 + grad1'*(z-x) + g2 + grad2'*(z-x);
            cvx_end
            if z'*Q*z+z'*q<funv
                funv = z'*Q*z+z'*q;
                x_new = z;
            end
        end
    end
    
    if norm(x_new - x)<0.01
        x = x_new;
        break
    else
        x = x_new;
    end
    
end


end
