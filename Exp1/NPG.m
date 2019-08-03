function [x,val] = NPG(A0,a0,A1,a1,aa1,A2,a2,aa2,rho,delta,x)

% Apply gradient with BB-step size to solve the problem:
% min x'*A0*x + a0'*x + rho*[x'*A1*x + a1'*x + aa1]_+^2 + rho*[x'*A2*x +
% a2'*x + aa2]_+^2

M = 5;
iter = 0;
alpha = 1;
gamma = 0.5;
Fc = zeros(M,1);
[val,grad] = fun_val(A0,a0,A1,a1,aa1,A2,a2,aa2,rho,x);
res = norm(grad);

while (res>delta)&&(iter<=40000)
    
    i = mod(iter,M)+1;
    Fc(i) = val;
    while 1
        x_new = x - alpha*grad;
        [val_new,grad_new] = fun_val(A0,a0,A1,a1,aa1,A2,a2,aa2,rho,x_new);
        if val_new <= max(Fc) - gamma*alpha*res^2
            break
        end
        alpha = 0.1*alpha;
    end
    dx = x - x_new;
    dg = grad - grad_new;
    alpha = min(1e+10,max(1e-10,abs(dx'*dg)/(dg'*dg)));
    x = x_new;
    val = val_new;
    grad = grad_new;
    res = norm(grad);
    iter = iter+1;
end

%     fprintf('iter = %d \t res = %1.2e \t funv = %1.2e \t alpha = %1.2e\n',iter,res,val,alpha);


end



function [val,grad] = fun_val(A0,a0,A1,a1,aa1,A2,a2,aa2,rho,x)


p1 = max(0,x'*A1*x + a1'*x + aa1);
p2 = max(0,x'*A2*x + a2'*x + aa2);
val = x'*A0*x + a0'*x + rho*p1^2 + rho*p2^2;
grad = 2*A0*x + a0 + 2*rho*p1*(2*A1*x+a1) + 2*rho*p2*(2*A2*x+a2);

end
