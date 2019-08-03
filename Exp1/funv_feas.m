function [funv,feas] = funv_feas(Q,q,A,a,ac,B,b,bc,x)

funv = x'*Q*x+q'*x;
fg1 = x'*A(:,:,1)*x+x'*a(:,1)+ac(1) - max(x'*B(:,:,1,1)*x+x'*b(:,1,1)+bc(1,1), x'*B(:,:,2,1)*x+x'*b(:,2,1)+bc(2,1));
feas1 = max(0,fg1);
fg2 = x'*A(:,:,2)*x+x'*a(:,2)+ac(2) - max(x'*B(:,:,1,2)*x+x'*b(:,1,2)+bc(1,2), x'*B(:,:,2,2)*x+x'*b(:,2,2)+bc(2,2));
feas2 = max(0,fg2);
% feas = [feas1,feas2];
feas = max(feas1,feas2);

end

