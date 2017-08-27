function I = LevyIntegral2(k,nodes,u,Br,Bl,N,K,S0,optionType,transf,barrierType)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
qnodes = linspace(Bl,Br,N+1);
dq = qnodes(2) - qnodes(1);
I = zeros([N-1 1]);
if nargin == 10 % european option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) + qnodes,nodes,u,K,S0,optionType,transf) .* k(qnodes));
	end
elseif nargin > 10 % barrier option
	for i = 1 : N-1
		I(i) = dq * trapz(f_u(nodes(i) + qnodes,nodes,u,K,S0,optionType,transf,barrierType) .* k(qnodes));
	end
end

end % function

