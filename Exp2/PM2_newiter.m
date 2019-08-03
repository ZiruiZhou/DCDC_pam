function x_new=PM2_newiter(A,b,rho,s,c,x0,eta)

epsilon=1e-4;
x=x0;
n=length(x);
delta = eta;

for iter=1:100
    
    I=zeros(n,3);
    dxt=max(x-s,0)+max(0,-x-s)-epsilon;
    I(x-s>=dxt,1) = 1;
    I(0>=dxt,2) = 1;
    I(-x-s>=dxt,3) = 1;
    Is = sum(I,2);
    if any(Is>2)
        error('found an element in Is greater than 2')
    end

    num = sum(Is==2);
    Num = 2^num;
    if Num>32
        fprintf('Num = %d in some iteration of QPM\n', Num);
    end
    Gx = [ones(n,1),zeros(n,1),-ones(n,1)]; 
    nzg = sum(I.*Gx,2);
    gx = zeros(n,1);
    gx(Is==1) = nzg(Is==1);

    funv = Inf;
    for i=1:Num
        i22 = bitget(i-1,num:-1:1)';
        gx(Is==2) = sign(x(Is==2)).*i22;
        hxs = sum(-x(gx<0)-s) + sum(x(gx>0)-s);
        gf = A'*(A*x-b);

%         cvx_begin quiet
%         cvx_precision low
%         variables z(n) t
%         minimize((z-x)'*(z-x)/2 + gf'*z + rho*t^2)
%         subject to
%             0 <= t;
%             norm(z,1) - hxs - gx'*(z-x) - c <= t;
%         cvx_end

        z=NPG(1,rho,delta,c,hxs,gf,gx,x);

        hzt = sum(max(z-s,0)+max(0,-z-s));
        F_rho = norm(A*z-b)^2/2 + rho*max(0,norm(z,1)-hzt-c)^2;
        if F_rho<funv
            x_new=z;
            funv=F_rho;
        end
    end
    
%     fprintf('QPM: Inner Iter = %d,\t F_rho = %1.8e\n',iter,funv);
    if norm(x_new-x)<eta
        break
    end

    x = x_new;
    delta = max(delta*0.1,1e-6);
end





    




