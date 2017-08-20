% runAmericanOption
clc;close all; clear variables;
r=0.01; sigma=0.2; T=2; K=110; S0=100;
rng(1992);
Nsim = 1e5; Nstep = 40;N = 300; M = 100;
omega = 1.5; tol = 1e-3; maxiter = 20;
NumMethod = 'CN';
theta = 1/2;
model = 'BS';
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
	case 'Kou'
		param = [0.6;0.01;9;6;0.3];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
	case 'VG'
		param = [0.2;0.2];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
	case 'NIG'
		param = [0.2; 0.8];
		% 1) MC
		[PriceMC,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep);
	otherwise
end
