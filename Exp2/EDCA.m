function x_new = EDCA(A,b,xinit,s,c)

epsilon=1e-4;
n=size(A,2);

% Generate an initial point by solving convex relaxation
if isempty(xinit)
    cvx_begin quiet
    variable x(n)
    minimize((A*x-b)'*(A*x-b))
    subject to
        norm(x,1) - c <= 0;
    cvx_end
else
    x = xinit;
end

if norm(x,1) - sum(max(0,x-s)) - sum(max(0,-x-s)) - c > 0
    error('The input of EDCA must be feasible.')
end

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
        fprintf('Num = %d in some iteration of EDCA\n', Num);
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
        
        cvx_begin quiet
        variable z(n)
        minimize((A*z-b)'*(A*z-b))
        subject to
            norm(z,1) - hxs - gx'*(z-x) - c <= 0;
        cvx_end
        
        ff = norm(A*z-b)^2;
        if ff<funv
            x_new = z;
            funv = ff;
        end
    end
    
%     fprintf('EDCA: Iter = %d, \t FF = %1.3e\n',iter,funv);
    if norm(x_new-x)<1e-2
        break
    end
    
    x = x_new;
end