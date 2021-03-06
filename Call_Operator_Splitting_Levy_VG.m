function price=Call_Operator_Splitting_Levy_VG
%% Call_Operator_Splitting_Levy
% European Call, VG process

% Model parameters
T = 1;S0 = 30;r = 0.1;K = 30;
nu=0.16; etaN=0.0694; etaP=0.0166; %Levy density parameters
% param = [6.25;60.2;14.4];
% Discretization parameters
N = 200; M = 100;
% Grids
dt=T/M;
Smin=10; Smax=300;
xmin=log(Smin/K); xmax=log(Smax/K);
x=linspace(xmin,xmax,N+1); dx=x(2)-x(1);
nodes=x(2:end-1)'; % interior nodes

% Compute sigma_Eps, alpha, lambda
epsilon=0.01;
k=@(y) 1./(nu*abs(y)).*exp(-abs(y)/etaN).*(y<0)+...
    1./(nu*abs(y)).*exp(-abs(y)/etaP).*(y>0);
% model = 'VG2';
% k = LevyDensity(param, model);
[alpha,lambda,sigma,Bl,Br]=Levy_Integrals(k,N,epsilon)

% Matrix
Dd=-(sigma^2/(2*(dx^2))+(sigma^2/2+alpha)/(2*dx));
Dm=1/dt-(-sigma^2/(dx^2)-lambda);
Du=-(sigma^2/(2*(dx^2))-(sigma^2/2+alpha)/(2*dx));
e=ones(N-1,1);
M1=spdiags([Dd*e Dm*e Du*e],-1:1,N-1,N-1);

% Forward-in-time loop
u=max( exp(nodes)-1, 0); % initial solution
BC=zeros(N-1,1);
for j=1:M
    % Boundary condition
    BC(end)=-Du*(exp(xmax)-1);
    % Compute I
    I=Levy_Int(k,nodes,u,Bl,Br,N,epsilon);
    %Solve the linear system
    u=M1\(u/dt+BC+I);
end
% Plot and interpolation
S=K*exp(nodes-r*T);
C=K*u*exp(-r*T);
figure
plot(S,C);
price=interp1(S,C,S0,'spline')
if lambda==0
    price_ex=blsprice(S0,K,r,T,sigma)
end
end

function [alpha,lambda,sigma,Bl,Br]=Levy_Integrals(k,N,epsilon)
Tol=10^-12; Bl=-0.01; Br=0.01;
% while (k(Bl)>Tol)
%     Bl=Bl-0.002;
% end
% while (k(Br)>Tol)
%     Br=Br+0.002;
% end
Bl = fsolve(@(y) k(y) - Tol, Bl);
Br = fsolve(@(y) k(y) - Tol, Br);
qnodes=linspace(-epsilon,epsilon,2*N);
sigma=trapz(qnodes,(qnodes.^2).*k(qnodes));
qnodes1=linspace(Bl,-epsilon,N);
qnodes2=linspace(epsilon,Br,N);
lambda=trapz(qnodes1,k(qnodes1))+trapz(qnodes2,k(qnodes2)); 
alpha=trapz(qnodes1,(exp(qnodes1)-1).*k(qnodes1))+...
    trapz(qnodes2,(exp(qnodes2)-1).*k(qnodes2));
figure; hold on
plot(qnodes1,k(qnodes1),'b');
plot(qnodes2,k(qnodes2),'b');
plot(qnodes,k(qnodes),'r');


end

function I=Levy_Int(k,nodes,u,Bl,Br,N,epsilon)
qnodes1=linspace(Bl,-epsilon,round(N/2)); 
dq1=qnodes1(2)-qnodes1(1);
qnodes2=linspace(epsilon,Br,round(N/2)); 
dq2=qnodes2(2)-qnodes2(1);
I=zeros(N-1,1);
for i=1:N-1
  I(i)=trapz(f_u(nodes(i)+qnodes1,nodes,u).*k(qnodes1))*dq1+...
      trapz(f_u(nodes(i)+qnodes2,nodes,u).*k(qnodes2))*dq2;
end
end

function f=f_u(z,nodes,u)
f=zeros(size(z));
index=find(z>nodes(end));
f(index)=exp(z(index))-1;
index=find( (z>=nodes(1)).*(z<=nodes(end)) );
f(index)=interp1( nodes,u,z(index));









    


end




