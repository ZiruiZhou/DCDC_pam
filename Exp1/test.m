% Matlab codes for the methods PM, ALM, EPM, and EDCA for solving problem
% (5.1) in the paper.
% PM: our penalty method;
%   - PM1: choose p=1;
%   - PM2: choose p=2;
% ALM: our augmented Lagrangian method;
% EPM: exact penalty method in [Pang et al.'17];
% EDCA: enhanced DCA in [Pang et al.'17].

clear all

dim = [50,100,250,500];     % Problem dimensions

grp = 1;    % 0 or 1: different types of generation

filename = 'computational result.txt';
fid = fopen(filename,'w');


for nn=1:4
    
    randn('seed',1);
    rand('seed',1);
    
    % Generate a problem instance
    n = dim(nn);
    I = 2;
    A = zeros(n,n,I);
    a = zeros(n,I);
    ac = zeros(1,I);
    B = zeros(n,n,2,I);
    b = zeros(n,2,I);
    bc = zeros(2,I);
    
    U = randn(n,n);
    U = orth(U);
    diag_Q = (1+20*rand(n,1));
    Q = U*diag(diag_Q)*U';
    if grp==0
        q = zeros(n,1);    
    else
        q = randn(n,1);
    end

    for i=1:I
        U = randn(n,n);
        U = orth(U);
        diag_A = (1+20*rand(n,1));
        A(:,:,i) = U*diag(diag_A)*U';
        a(:,i) = randn(n,1);
        if grp==0
            ac(i) = 0;
        else
            ac(i) = randn;
        end
        for j=1:2
            U = randn(n,n);
            U = orth(U);
            diag_B = (1+20*rand(n,1));
            B(:,:,j,i) = U*diag(diag_B)*U';
            b(:,j,i) = randn(n,1);
            if grp==0
                bc(j,i) = 0;
            else
                bc(j,i) = randn;
            end
        end
    end
    
    % Set storage for num_test number of tests
    num_test = 10;

    x_init = zeros(n,num_test);
    funv_init = zeros(1,num_test);
    feas_init = zeros(1,num_test);

    x_epm = zeros(n,num_test);
    funv_epm = zeros(1,num_test);
    feas_epm = zeros(1,num_test);
    time_epm = zeros(1,num_test);

    x_pm1 = zeros(n,num_test);
    funv_pm1 = zeros(1,num_test);
    feas_pm1 = zeros(1,num_test);
    time_pm1 = zeros(1,num_test);

    x_pm2 = zeros(n,num_test);
    funv_pm2 = zeros(1,num_test);
    feas_pm2 = zeros(1,num_test);
    time_pm2 = zeros(1,num_test);

    x_edca = zeros(n,num_test);
    funv_edca = zeros(1,num_test);
    feas_edca = zeros(1,num_test);
    time_edca = zeros(1,num_test);

    x_alm = zeros(n,num_test);
    funv_alm = zeros(1,num_test);
    feas_alm = zeros(1,num_test);
    time_alm = zeros(1,num_test);

    for tt = 1:num_test
        
        fprintf('\n\n Running experiment for n = %d and tt = %d:\n', n,tt);
        fprintf('----------------------------------------------------------\n');


        % Generate an initial point
        x = randn(n,1);
        [~,feas_x] = funv_feas(Q,q,A,a,ac,B,b,bc,x);
        gen = 0;
        while feas_x>0
            x = randn(n,1);
            [~,feas_x] = funv_feas(Q,q,A,a,ac,B,b,bc,x);
            gen = gen + 1;
            if gen>100
                error('There are no feasible points found after 100 times random generation.')
            end
        end
        if feas_x>0
            error('The input of EDCA must be a feasible point.')
        end
        x_init(:,tt) = x;

        [funv_init(tt),feas_init(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x);

        eps = 0.01;
        rho = 10;
        

        % Solve by our penalty method with p=1
        estart = tic;
        x_pm1(:,tt) = PM1(Q,q,A,a,ac,B,b,bc,rho,eps,x,10);
        time_pm1(tt) = toc(estart);
        [funv_pm1(tt),feas_pm1(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x_pm1(:,tt));
        fprintf('PM1:  funv=%1.4f feas=%1.4e time=%1.4f\n',funv_pm1(tt),feas_pm1(tt),time_pm1(tt));


        % Solve by our penalty method with p=2.
        qstart = tic;
        x_pm2(:,tt) = PM2(Q,q,A,a,ac,B,b,bc,rho,eps,x,10);
        time_pm2(tt) = toc(qstart);
        [funv_pm2(tt),feas_pm2(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x_pm2(:,tt));
        fprintf('PM2:  funv=%1.4f feas=%1.4e time=%1.4f\n',funv_pm2(tt),feas_pm2(tt),time_pm2(tt));


        % Solve by our ALM.
        qstart = tic;
        x_alm(:,tt) = ALM(Q,q,A,a,ac,B,b,bc,rho,eps,x,10);
        time_alm(tt) = toc(qstart);
        [funv_alm(tt),feas_alm(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x_alm(:,tt));
        fprintf('ALM:  funv=%1.4f feas=%1.4e time=%1.4f\n',funv_alm(tt),feas_alm(tt),time_alm(tt));
        
        % Solve by the exact penalty method in Pang et al.'16.
        estart = tic;
        x_epm(:,tt) = EPM(Q,q,A,a,ac,B,b,bc,rho,eps,x,10,estart);
        time_epm(tt) = toc(estart);
        [funv_epm(tt),feas_epm(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x_epm(:,tt));
        fprintf('EPM:  funv=%1.4f feas=%1.4e time=%1.4f\n',funv_epm(tt),feas_epm(tt),time_epm(tt));

        % Solve by EDCA in Pang et al.'16.
        edstart = tic;
        x_edca(:,tt) = EDCA(Q,q,A,a,ac,B,b,bc,eps,x,100);
        time_edca(tt) = toc(edstart);
        [funv_edca(tt),feas_edca(tt)] = funv_feas(Q,q,A,a,ac,B,b,bc,x_edca(:,tt));
        fprintf('EDCA: funv=%1.4f feas=%1.4e time=%1.4f\n',funv_edca(tt),feas_edca(tt),time_edca(tt));
        
        fprintf('----------------------------------------------------------\n');

    end
        
    filenames = ['n=' num2str(n)];
    save(filenames);
    
    
    fprintf(fid,'\n\nAveraged result for n = %d:\n', n);
    fprintf(fid,'----------------------------------------------------------\n');
    fprintf(fid,'PM1:  funv=%1.4f feas=%1.4e time=%1.4f\n',sum(funv_pm1)/num_test,sum(feas_pm1)/num_test,sum(time_pm1)/num_test);
    fprintf(fid,'PM2:  funv=%1.4f feas=%1.4e time=%1.4f\n',sum(funv_pm2)/num_test,sum(feas_pm2)/num_test,sum(time_pm2)/num_test);
    fprintf(fid,'ALM:  funv=%1.4f feas=%1.4e time=%1.4f\n',sum(funv_alm)/num_test,sum(feas_alm)/num_test,sum(time_alm)/num_test);
    fprintf(fid,'EPM:  funv=%1.4f feas=%1.4e time=%1.4f\n',sum(funv_epm)/num_test,sum(feas_epm)/num_test,sum(time_epm)/num_test);
    fprintf(fid,'EDCA: funv=%1.4f feas=%1.4e time=%1.4f\n',sum(funv_edca)/num_test,sum(feas_edca)/num_test,sum(time_edca)/num_test);
    fprintf(fid,'----------------------------------------------------------\n');

end

fclose(fid);


