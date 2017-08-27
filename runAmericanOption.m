% runAmericanOption
clc;close all; clear variables;
T = 1; K = 1; S0 = 1; r = 0.001;rng(1992);
Nsim = 3e5; Nstep = 200;N = 300; M = 200;
% the goodness of the approximation stringly depend on the tol
omega = 1.5; tol = 1e-5; maxiter = 110;
epsilon = 0.01;
NumMethod = 'implicit';
theta = 1/2;
model = 'Kou';
optionType = 'Put';
%% MC
switch model
	case 'BS'
		param = 0.6;
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
		% 2) PDE FD
		[PriceFD] = FDLogPriceAmerican1( S0,K,r,T,N,M,param,NumMethod,theta,...
			maxiter,tol,omega);
		% 3) PDE FEM
		[PriceFEM] = FEMLogPriceAmerican1( S0,K,r,T,N,M,param,NumMethod,theta,...
			maxiter,tol,omega);
		% 4) PDE FEM full implementation
		[PriceFEM2] = FEMLogPriceAmerican2( S0,K,r,T,N,M,param,NumMethod,theta,...
			maxiter,tol,omega);
	case 'Merton'
		param = [0.5;0.6;0.1;0.1];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
		% 2) PIDE logmoneyness
		[pricePIDElogmon] = FDLevyLogMoneynessAmerican1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
			maxiter,tol,omega);
	case 'Kou'
		param = [0.120381;0.330966;9.65997;3.13868;0.2761];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
		% 2) PIDE logmoneyness
		[pricePIDElogmon] = FDLevyLogMoneynessAmerican1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
			maxiter,tol,omega);
	case 'VG'
		param = [0.2;0.2];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
		% 2) PIDE
		[PricePIDElogmon] = FDLevyLogMoneyness1AIAmerican( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
			epsilon,maxiter,tol,omega);
	case 'NIG'
		param = [0.2; 0.8];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
		[price] = FDLevyLogMoneyness1AIAmerican( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
			epsilon,maxiter,tol,omega);
	otherwise
end
