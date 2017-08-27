% runBarrierOption
% pricing barrier option
% REMARK : the MC approximation strongly depends on the time horizon, if
% it's short (T<1) the approxiamtion is good if it's large the
% approximation is really poor
clc; close all; clear variables;
K=100; S0=100; r=0.3; T=1;
Nsim = 2e5; Nstep = 200;N = 500; M = 552;
rng(1992);
model = 'Kou';
optionType = 'Call';
barrierType = 'DO';
NumMethod = 'implicit';
theta = 1;
barrier = 70;
epsilon = 0.01;
switch optionType
	case 'Call'
		switch barrierType
			case 'DO'
				switch model
					case 'BS'
						param = 0.6;
						% 1) closed formula
						[ PriceExact] = barrierBS( S0, r, T, param,barrier, K, barrierType);
						% 2) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 3) PDE LogPrice implementation short
						[pricePDELogPrice] = FDLogPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						% 4) PDE LogPrice implementation full
						[pricePDELogPrice2] = FDLogPrice2( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						% 5) PDE Price short implementation
						[pricePDEPrice] = FDPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						% 6) PDE FEM LopgPrice
						[priceFEM] = FEMLogPrice2( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Merton'
						% sigma, lambda, mu, delta
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 2) PIDE logmoneyness --> don't use this with barrier
						% option since the barrier is time-dependent
% 						[pricePIDElogmon] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
% 							barrierType,barrier);
						% 3) PIDE logprice
						[pricePIDElogprice] = FDLevyLogPrice1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Kou'
						param = [0.6;0.5;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 3) PIDE logprice short implementation
						[pricePIDElogprice] = FDLevyLogPrice1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'VG'
						% theta_VG, sigma_VG, kappa
						param = [0.03;0.12;0.2;0];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 3) PIDE logpriceFwdPrice
						[pricePIDEprice] = FDLevyLogpriceFwdprice1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							epsilon,barrierType,barrier);
					case 'NIG'
						param = [0.2;0.8];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					otherwise
				end
			case 'UO'
				switch model
					case 'BS'
						param = 0.6;
						% 2) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 5) PDE Price short implementation
						[pricePDEPrice] = FDPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Merton'
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'VG'
						% theta_VG, sigma_VG, kappa
						param = [0.03;0.12;0.2;0];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'NIG'
						param = [0.2;0.8];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					otherwise
				end
			otherwise
		end
		
	case 'Put'
		switch barrierType
			case 'DO'
				switch model
					case 'BS'
						param = 0.6;
						% 2) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 5) PDE Price short implementation
						[pricePDEPrice] = FDPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						model = 'Merton'; param = [0.6;0.;0.1;0.1];
					case 'Merton'
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'VG'
						% theta_VG, sigma_VG, kappa
						param = [0.03;0.12;0.2;0];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'NIG'
						param = [0.2;0.8];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					otherwise
				end
			case 'UO'
				switch model
					case 'BS'
						param = 0.6;
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 2) PDE Price short implementation
						[pricePDEPrice] = FDPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						% 3) PDE logprice
						[pricePDElogprice] = FDLogPrice1( S0,K,r,T,N,M,param,optionType,NumMethod,theta,...
							barrierType,barrier);
						% check model merton
						model = 'Merton'; param = [0.6;0.;0.1;0.1];
					case 'Merton'
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					case 'VG'
						% theta_VG, sigma_VG, kappa
						param = [0.03;0.12;0.2;0.2];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 3) PIDE logpriceFwdPrice
						[pricePIDEprice] = FDLevyLogpriceFwdprice1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							epsilon,barrierType,barrier);
					case 'NIG'
						param = [0.2;0.8];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
					otherwise
				end
			otherwise
		end
end
