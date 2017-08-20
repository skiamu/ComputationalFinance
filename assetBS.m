function [ S, check] = assetBS( S0,r,T,sigma,Nsim,Nstep )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
dt = T / Nstep;
dX = zeros([Nsim Nstep+1]);
Z = randn([Nsim Nstep]);
dX(:,2:end) = (r - sigma^2/2) * dt + sigma * sqrt(dt) * Z;
X = cumsum(dX,2);
S = S0 .* exp(X);
check = mean(S(:,end) / exp(r * T));
end

