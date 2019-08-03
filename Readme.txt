Aug 1, 2019

This subdirectory contains the following MATLAB source codes and folders:


./Exp1: This subdirectory contains codes for the first experiment (solving problem (5.1)) in our paper.	

	ALM.m			main iteration of our augmented Lagrangian method

	ALM_newiter.m		solving the augmented Lagrangian subproblem by Algorithm 3.2 in our paper

	EDCA.m			enhanced DCA in [Pang et al.'17]

	EPM.m			main iteration of the exact penalty method in [Pang et al.'17]

	EPM_newiter.m		solving the exact penalty subproblem by Algorithm 1 of [Pang et al.'17]

	funv_eval.m		evaluate the component functions in problem (5.1)

	funv_feas.m		evaluate the objective value and constraint violation	

	NPG.m			nonmonotone gradient method to solve the convex subproblem arising in ALM_newiter.m and PM2_newiter.m

	PM1.m			main iteration of our penalty method with p=1

	PM1_newiter.m		solving the subproblem of our penalty method with p=1

	PM2.m			main iteration of our penalty method with p=2

	PM2_newiter.m		solving the subproblem of our penalty method with p=2

	test.m			run code; generate problem instances and record computational results.



./Exp2: This subdirectory contains codes for the second experiment (solving problem (5.4)) in our paper.	

	ALM.m			main iteration of our augmented Lagrangian method

	ALM_newiter.m		solving the augmented Lagrangian subproblem by Algorithm 3.2 in our paper

	EDCA.m			enhanced DCA in [Pang et al.'17]

	EPM.m			main iteration of the exact penalty method in [Pang et al.'17]

	EPM_newiter.m		solving the exact penalty subproblem by Algorithm 1 of [Pang et al.'17]

	funv_feas.m		evaluate the objective value and constraint violation	

	NPG.m			nonmonotone gradient method to solve the convex subproblem arising in ALM_newiter.m and PM2_newiter.m

	PM1.m			main iteration of our penalty method with p=1

	PM1_newiter.m		solving the subproblem of our penalty method with p=1

	PM2.m			main iteration of our penalty method with p=2

	PM2_newiter.m		solving the subproblem of our penalty method with p=2

	test.m			run code; generate problem instances and record computational results.



The runcodes require the package CVX, downloadable at http://cvxr.com/cvx/

Implementation and numerical experience with the above codes are described in Section 5 of the paper: 
    Zhaosong Lu, Zhe Sun, and Zirui Zhou
    "Penalty and Augmented Lagrangian Methods for a Class of Structured Nonsmooth DC Constrained DC Program",
    Submitted.
This code was last updated on Aug 1, 2019.

Questions/comments/suggestions about the codes are welcome.  

Zhaosong Lu, 	zhaosong@sfu.ca
Zhe Sun,	snzma@126.com
Zirui Zhou,	zirui-zhou@hkbu.edu.hk