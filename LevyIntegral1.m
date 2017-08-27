function [alpha,lambda_num,Bl,Br] = LevyIntegral1(k,N)
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here
tol = 1e-10;
%% method 1
Br = 10;Bl = -10;
while k(Bl) < tol
	Bl = Bl + 0.02;
end
while k(Br) < tol
	Br = Br - 0.02;
end
%% method 2
% this methos works well for Merton and Kou
% Bl = -1;Br = 1;
% while (k(Bl)>tol)
%     Bl=Bl-2;
% end
% while (k(Br)>tol)
%     Br=Br+2;
% end
%% method 3 
% this methos works well for inifinite activity process (VG)
% Bl = -0.01; Br = 0.01;
% Bl = fsolve(@(y) k(y) - tol, Bl);
% Br = fsolve(@(y) k(y) - tol, Br);
qnodes = linspace(Bl, Br, 2*N);
lambda_num = trapz(qnodes,k(qnodes));
alpha = trapz(qnodes, (exp(qnodes) - 1) .* k(qnodes));
end % end function

