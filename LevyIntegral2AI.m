function I = LevyIntegral2AI(k,nodes,u,Br,Bl,N,K,S0,epsilon,optionType,transf,barrierType)
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here
qnodes1 = linspace(Bl,-epsilon,round(N/2));
qnodes2 = linspace(epsilon,Br,round(N/2));
dq1 = qnodes1(2) - qnodes1(1);
dq2 = qnodes2(2) - qnodes2(1);
I = zeros([N-1 1]);
if nargin == 11 % european option
	for i = 1 : N-1
		I(i) = dq1 * trapz(f_u(nodes(i) + qnodes1,nodes,u,K,S0,optionType,transf) .* k(qnodes1)) + ...
			dq2 * trapz(f_u(nodes(i) + qnodes2,nodes,u,K,S0,optionType,transf) .* k(qnodes2));
	end
elseif nargin > 11 % barrier option
	for i = 1 : N-1
		I(i) = dq1 * trapz(f_u(nodes(i) + qnodes1,nodes,u,K,S0,optionType,transf,barrierType) .* k(qnodes1)) + ...
			dq2 * trapz(f_u(nodes(i) + qnodes2,nodes,u,K,S0,optionType,transf,barrierType) .* k(qnodes2));
	end
end

end % function

