% Matlab codes for the methods PM, ALM, EPM, and EDCA for solving problem
% (5.4) in the paper.
% PM: our penalty method;
%   - PM1: choose p=1;
%   - PM2: choose p=2;
% ALM: our augmented Lagrangian method;
% EPM: exact penalty method in [Pang et al.'17];
% EDCA: enhanced DCA in [Pang et al.'17].

clear all
close all

filename = 'computational result.txt';
fid = fopen(filename,'w');

nr = [2^8,2^9,2^10,2^11];
mr = [2^6,2^7,2^8,2^9];
num_dim = length(nr);
Kr = [5,10];
num_K = 2;

num_test = 10;


for kk = 1:num_K
    for nn = 1:num_dim

        n = nr(nn);
        m = mr(nn);
        K = Kr(kk);
        s = 0.1; 
        c=K*s;
        
        x_cvx = zeros(n,num_test);
        funv_cvx = zeros(1,num_test);
        feas_cvx = zeros(1,num_test);
        err_cvx = zeros(1,num_test);
        time_cvx = zeros(1,num_test);
        
        x_pm1 = zeros(n,num_test);
        funv_pm1 = zeros(1,num_test);
        feas_pm1 = zeros(1,num_test);
        err_pm1 = zeros(1,num_test);
        time_pm1 = zeros(1,num_test);

        x_pm2 = zeros(n,num_test);
        funv_pm2 = zeros(1,num_test);
        feas_pm2 = zeros(1,num_test);
        err_pm2 = zeros(1,num_test);
        time_pm2 = zeros(1,num_test);

        x_alm = zeros(n,num_test);
        funv_alm = zeros(1,num_test);
        feas_alm = zeros(1,num_test);
        err_alm = zeros(1,num_test);
        time_alm = zeros(1,num_test);

        x_epm = zeros(n,num_test);
        funv_epm = zeros(1,num_test);
        feas_epm = zeros(1,num_test);
        err_epm = zeros(1,num_test);
        time_epm = zeros(1,num_test);

        x_edca = zeros(n,num_test);
        funv_edca = zeros(1,num_test);
        feas_edca = zeros(1,num_test);
        err_edca = zeros(1,num_test);
        time_edca = zeros(1,num_test);

        for nnt=1:num_test
            
            fprintf('\n\nExperiment on m = %d, n = %d, K = %d, \t No. test = %d.\n', m,n,K,nnt);
            fprintf('----------------------------------------------------------\n');
     
            randn('seed',nnt);
            rand('seed',nnt);
            
            % random +/- 1 signal
            xs = zeros(n,1);
            q = randperm(n);
            xs(q(1:K)) = sign(randn(K,1));

            A = randn(m,n);
            A = orth(A')';

            % noisy observations
            sigma = 0.001; 
            b = A*xs + sigma*randn(m,1);
            
            % Solving the convex relaxation to generate an initial point
            st_cvx = tic;
            cvx_begin quiet
            variable x(n)
            minimize((A*x-b)'*(A*x-b))
            subject to
                norm(x,1) - c <= 0;
            cvx_end
            
            % See the performance of the solution by convex relaxation
            time_cvx(nnt) = toc(st_cvx);
            x_cvx(:,nnt) = x;
            [funv_cvx(nnt),feas_cvx(nnt),err_cvx(nnt)] = funv_feas_err(A,b,x,xs,s,c);
            fprintf('CVX:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_cvx(nnt),feas_cvx(nnt),err_cvx(nnt),time_cvx(nnt));
            
            
            xinit = x;
            
            % Solve (5.4) by PM1
            st_pm1 = tic;
            xs_pm1 = PM1(A,b,xinit,s,c);
            time_pm1(nnt) = toc(st_pm1);
            x_pm1(:,nnt) = xs_pm1;
            [funv_pm1(nnt),feas_pm1(nnt),err_pm1(nnt)] = funv_feas_err(A,b,xs_pm1,xs,s,c);
            fprintf('PM1:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_pm1(nnt),feas_pm1(nnt),err_pm1(nnt),time_pm1(nnt));

            % Solve (5.4) by PM2          
            st_pm2 = tic;
            xs_pm2 = PM2(A,b,xinit,s,c);
            time_pm2(nnt) = toc(st_pm2);
            x_pm2(:,nnt) = xs_pm2;
            [funv_pm2(nnt),feas_pm2(nnt),err_pm2(nnt)] = funv_feas_err(A,b,xs_pm2,xs,s,c);
            fprintf('PM2:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_pm2(nnt),feas_pm2(nnt),err_pm2(nnt),time_pm2(nnt));

            % Solve (5.4) by ALM
            st_alm = tic;
            xs_alm = ALM(A,b,xinit,s,c);
            time_alm(nnt) = toc(st_alm);
            x_alm(:,nnt) = xs_alm;
            [funv_alm(nnt),feas_alm(nnt),err_alm(nnt)] = funv_feas_err(A,b,xs_alm,xs,s,c);
            fprintf('ALM:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_alm(nnt),feas_alm(nnt),err_alm(nnt),time_alm(nnt));

            % Solve (5.4) by EPM            
            st_epm = tic;
            xs_epm = EPM(A,b,xinit,s,c);
            time_epm(nnt) = toc(st_epm);
            x_epm(:,nnt) = xs_epm;
            [funv_epm(nnt),feas_epm(nnt),err_epm(nnt)] = funv_feas_err(A,b,xs_epm,xs,s,c);
            fprintf('EPM:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_epm(nnt),feas_epm(nnt),err_epm(nnt),time_epm(nnt));

            % Solve (5.4) by EDCA
            st_edca = tic;
            xs_edca = EDCA(A,b,xinit,s,c);
            time_edca(nnt) = toc(st_edca);
            x_edca(:,nnt) = xs_edca;
            [funv_edca(nnt),feas_edca(nnt),err_edca(nnt)] = funv_feas_err(A,b,xs_edca,xs,s,c);
            fprintf('EDCA: funv=%1.4e feas=%1.4e err=%1.4e time=%1.4f\n',funv_edca(nnt),feas_edca(nnt),err_edca(nnt),time_edca(nnt));
            
            fprintf('----------------------------------------------------------\n');

        end

        filenames = ['m=' num2str(m) 'n=' num2str(n) 'K=' num2str(K)];
        save(filenames);
        
        fprintf(fid,'\n\nAveraged result for m = %d, n = %d, K = %d\n', m,n,K);
        fprintf(fid,'----------------------------------------------------------\n');
        fprintf(fid,'CVX:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_cvx)/num_test,sum(feas_cvx)/num_test,sum(err_cvx)/num_test,sum(time_cvx)/num_test);
        fprintf(fid,'PM1:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_pm1)/num_test,sum(feas_pm1)/num_test,sum(err_pm1)/num_test,sum(time_pm1)/num_test);
        fprintf(fid,'PM2:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_pm2)/num_test,sum(feas_pm2)/num_test,sum(err_pm2)/num_test,sum(time_pm2)/num_test);
        fprintf(fid,'ALM:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_alm)/num_test,sum(feas_alm)/num_test,sum(err_alm)/num_test,sum(time_alm)/num_test);
        fprintf(fid,'EPM:  funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_epm)/num_test,sum(feas_epm)/num_test,sum(err_epm)/num_test,sum(time_epm)/num_test);
        fprintf(fid,'EDCA: funv=%1.4e feas=%1.4e err=%1.4e time=%1.4e\n',sum(funv_edca)/num_test,sum(feas_edca)/num_test,sum(err_edca)/num_test,sum(time_edca)/num_test);
        fprintf(fid,'----------------------------------------------------------\n');

    end
end

fclose(fid);
