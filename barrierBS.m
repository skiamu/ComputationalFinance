function [ price] = barrierBS( S0, r, T, sigma,L, K, barrierType)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
r_tilde = r - 0.5 * sigma^2;
% european call option price
C = blsprice(S0,K,r,T,sigma);
switch barrierType
	case 'DO'
		C_tilde = blsprice(L^2/S0,K,r,T,sigma);
		price = C - (L/S0)^(2*r_tilde/sigma^2) * C_tilde;
	case 'UO'
		
	case 'DI'
		C_tilde = blsprice(L^2/S0.K,r,T.sigma);
		price = (L/S0)^(2*r_tilde/sigma^2) * C_tilde;
	case 'UI'
	otherwise

end

