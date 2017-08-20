function [Price, check, IC] = PricingOptionMC( S0,r,T,K,param,Nsim,Nstep,...
	model,optionType, barrierType, barrier)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 9
	switch optionType
		case 'Call'
			Payoff = @(x) subplus(x(:,end) - K);
		case 'Put'
			Payoff = @(x) subplus(K - x(:,end));
		otherwise
	end
elseif nargin > 9
	[ Payoff ] = optionPayoff(K,optionType,barrierType, barrier);
end
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
DiscountedPayoff = exp(-r * T) * Payoff(S);
[Price, ~, IC] = normfit(DiscountedPayoff);

end

