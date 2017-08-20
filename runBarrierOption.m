% runBarrierOption
% pricing barrier option
clc; close all; clear variables;rng(1192)
T = 1;S0 = 30;r = 0.1;K = 30;
Nsim = 2e5; Nstep = 80;N = 400; M = 400;
model = 'VG';
optionType = 'Call';
barrierType = 'DO';
NumMethod = 'implicit';
theta = 1;
barrier = 15;
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
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'VG'
						param = [0.2;0.2];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
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
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'VG'
						param = [0.2;0.2];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
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
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Merton'
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'VG'
						param = [0.2;0.2];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
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
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Merton'
						param = [0.5;0.6;0.1;0.1];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'Kou'
						param = [0.6;0.01;9;6;0.3];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
							barrierType,barrier);
					case 'VG'
						param = [0.2;0.2];
						% 1) MC
						[PriceMC, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
							model,optionType, barrierType, barrier);
						% 4) PIDE
						[pricePIDE] = FDLevyLogMoneyness1AI( S0,K,r,T,N,M,param,model,optionType,NumMethod,theta,...
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
