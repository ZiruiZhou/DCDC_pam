function x = PM2_newiter(Q,q,A,a,ac,B,b,bc,rho,eta,eps,x0)

x = x0;
n = length(x);
delta = eta;

for k=1:1000
%     fprintf('QPM: Inner iteration = %d\n',k);

    g11 = x'*B(:,:,1,1)*x + x'*b(:,1,1) + bc(1,1);
    g12 = x'*B(:,:,2,1)*x + x'*b(:,2,1) + bc(2,1);
    g21 = x'*B(:,:,1,2)*x + x'*b(:,1,2) + bc(1,2);
    g22 = x'*B(:,:,2,2)*x + x'*b(:,2,2) + bc(2,2);
    
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
    bc1 = zeros(1,nj1);
    for i=1:nj1
        B1(:,:,i) = B(:,:,j1(i),1);
        b1(:,i) = b(:,j1(i),1);
        bc1(i) = bc(j1(i),1);
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
    bc2 = zeros(1,nj2);
    for i=1:nj2
        B2(:,:,i) = B(:,:,j2(i),2);
        b2(:,i) = b(:,j2(i),2);
        bc2(i) = bc(j2(i),2);
    end
                
    funv = Inf;
    for J1=1:nj1
        for J2=1:nj2
            BB1 = B1(:,:,J1);
            bb1 = b1(:,J1);
            bbc1 = bc1(J1);
            BB2 = B2(:,:,J2);
            bb2 = b2(:,J2);
            bbc2 = bc2(J2);
            
            A0 = Q;
            a0 = q;
            A1 = A(:,:,1);
            a1 = a(:,1) - 2*BB1*x - bb1;
            aa1 = (2*BB1*x+bb1)'*x - x'*BB1*x - bb1'*x + ac(1) - bbc1;
            A2 = A(:,:,2);
            a2 = a(:,2) - 2*BB2*x - bb2;
            aa2 = (2*BB2*x+bb2)'*x - x'*BB2*x - bb2'*x + ac(2) - bbc2;
            
            
            [z,funv_n] = NPG(A0,a0,A1,a1,aa1,A2,a2,aa2,rho,delta,x);
            
            if funv_n<funv
                funv = funv_n;
                x_new = z;
            end
        end
    end
    
    if norm(x_new - x)<eta
        x = x_new;
        break
    else
        x = x_new;
        delta = max(0.5*delta,1e-6);
    end
    
end

end

