function [obj,f1,g11,g12,f2,g21,g22] = funv_eval(Q,q,A,a,ac,B,b,bc,x)

g11 = x'*B(:,:,1,1)*x + x'*b(:,1,1) + bc(1,1);
g12 = x'*B(:,:,2,1)*x + x'*b(:,2,1) + bc(2,1);
g21 = x'*B(:,:,1,2)*x + x'*b(:,1,2) + bc(1,2);
g22 = x'*B(:,:,2,2)*x + x'*b(:,2,2) + bc(2,2);
f1 =  x'*A(:,:,1)*x   + x'*a(:,1)   + ac(1);
f2 =  x'*A(:,:,2)*x   + x'*a(:,2)   + ac(2);
obj = x'*Q*x + q'*x;


end

