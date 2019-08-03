function [funv,feas,err] = funv_feas_err(A,b,x,xs,s,c)

funv = norm(A*x-b)^2;
feas = max(0, norm(x,1) - sum(max(x-s,0)) - sum(max(-x-s,0)) - c);
err = norm(x - xs)/norm(xs);


end

