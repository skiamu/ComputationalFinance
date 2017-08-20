function [ V ] = PSORalgorithm2(M1,b,V,maxiter,tol,omega,K,N,x)
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here
for i = 1 : maxiter
	Vold = V;
	for ii = 1 : N-1
		if ii == 1
			y = (b(ii) - M1(ii,ii+1) * Vold(ii+1)) / M1(ii,ii);
		elseif ii == N+1
			y = (b(ii) - M1(ii,ii-1) * V(ii-1)) / M1(ii,ii);
		else
			y = (b(ii) - M1(ii,ii+1) * Vold(ii+1) - M1(ii,ii-1) * V(ii-1)) / M1(ii,ii);
		end
		V(ii) = Vold(ii) + omega * (y - Vold(ii));
		% x is the full space grid while ii doesn't consider the extreme
		% points, hence x(ii+1). This is valid if the implementation is the
		% short one
		V(ii) = max(V(ii), K - exp(x(ii)));
	end
	
	if norm(V - Vold,inf) < tol
		break
	end	
end

end % end function

