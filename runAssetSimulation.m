% run asset function
clc;close all; clear variables;

T = 1;S0 = 30;r = 0.1;
Nsim = 5e5; Nstep = 30;

model = 'BS';

switch model
	case 'BS'
		sigma = 0.6;
		[ S, check] = assetBS(S0,r,T,sigma,Nsim,Nstep);
	case 'Merton'
		% lambda, sigma, mu, delta
		param = [0.5; 0.1859; 0.5864; 0.0015];
		[ S, check ] = assetMerton(S0,r,T,param,Nsim,Nstep);
	case 'Kou'
		% sigma,lambda, lambda_p, lambda_m, p
		param = [0.1619; 0.0211; 9.65997; 3.13868; 0.3155];
		[ S, check ] = assetKou(S0,r,T,param,Nsim,Nstep);
	case 'VG'
		% sigma, kappa
		param = [0.12;0.2];
		[ S, check ] = assetVG( S0,r,T,param,Nsim,Nstep );
	case 'NIG'
		% sigma, kappa
		param = [0.12;0.2];
		[ S, check ] = assetNIG( S0,r,T,param,Nsim,Nstep );
	otherwise
end