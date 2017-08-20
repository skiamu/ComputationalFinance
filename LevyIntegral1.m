function [alpha,lambda_num,Bl,Br] = LevyIntegral1(k,N)
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here
tol = 1e-10;
Br = 0.01;Bl = -0.01;
% while (k(Bl)>tol)
%     Bl=Bl-0.2;
% end
% while (k(Br)>tol)
%     Br=Br+0.2;
% end
Bl = fsolve(@(y) k(y) - tol, Bl);
Br = fsolve(@(y) k(y) - tol, Br);
qnodes = linspace(Bl, Br, 2*N);
lambda_num = trapz(qnodes,k(qnodes));
alpha = trapz(qnodes, (exp(qnodes) - 1) .* k(qnodes));
end % end function

