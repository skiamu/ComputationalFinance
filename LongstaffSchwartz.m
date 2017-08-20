function [Price,IC] = LongstaffSchwartz(S0,K,r,T,param,model,Nsim,Nstep)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
dt = T / Nstep;
switch model
	case 'BS'
		[ S, check] = assetBS( S0,r,T,param,Nsim,Nstep);
	case 'Merton'
		[ S, check ] = assetMerton(S0,r,T,param,Nsim,Nstep);
	case 'Kou'
		[ S, check ] = assetKou(S0,r,T,param,Nsim,Nstep);
	case 'VG'
		[ S, check ] = assetVG(S0,r,T,param,Nsim,Nstep);
	case 'NIG'
		[ S, check ] = assetNIG(S0,r,T,param,Nsim,Nstep);
end
%% 1) Initialization
ExerciseTime = Nstep * ones([Nsim 1]);
Put = subplus(K - S(:,end));
%% 2) backward procedure
for j = Nstep-1 : -1 : 1
	InMoney = find(S(:,j) < K);
	S_InMoney = S(InMoney,j);
	IV = K - S_InMoney;
	
	A = [ones(length(IV),1) S_InMoney S_InMoney.^2];
	b = exp(-r * (ExerciseTime(InMoney) - j) * dt) .* Put(InMoney);
	alpha = A \ b;
	CV = A * alpha;
	Index = find(IV > CV);
	EarlyExercise = InMoney(Index);
	Put(EarlyExercise) = IV(Index);
	ExerciseTime(EarlyExercise) = j;
end
[Price, ~, IC] = normfit(exp(-r * ExerciseTime * dt) .* Put);
end

