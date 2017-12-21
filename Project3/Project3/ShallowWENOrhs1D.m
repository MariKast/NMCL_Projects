function [dq] = ShallowWENOrhs1D(x,q,h,time,k,m,Crec,dw,beta,lambda_max)
% function [dq] = ShallowWENOrhs1D(x,q,h,k,m,Crec,dw,beta,maxvel)
% Purpose: Evaluate right hand side for Shallow water equations using a WENO method
% with reconstruction on the conserved variables
N = length(x)-1; 
dq = zeros(N,2); 
ql = zeros(N,2); 
qr = zeros(N,2); 
qp = zeros(N,2); 
qm = zeros(N,2);

% Extend data and assign boundary conditions
[xe,he] = extendP3(x,q(:,1),h,m,'P',0,'P',0);
[xe,me] = extendP3(x,q(:,2),h,m,'P',0,'P',0);

% define cell left and right interface values
hl = zeros(N+2,1); 
hr = zeros(N+2,1);
ml = zeros(N+2,1); 
mr = zeros(N+2,1);

for i=1:N+2
   [hl(i),hr(i)] = WENO(xe(i:(i+2*(m-1))),he(i:(i+2*(m-1))),m,Crec,dw,beta);
   [ml(i),mr(i)] = WENO(xe(i:(i+2*(m-1))),me(i:(i+2*(m-1))),m,Crec,dw,beta);
end

%Compute source term here
source = calcSourceTerm_ex(x,time);
%source =0;
% Compute rhs - also change numerical flux here
ql = [hl(2:N+1) ml(2:N+1)]; 
qr = [hr(2:N+1) mr(2:N+1)]; 
qp = [hl(3:N+2) ml(3:N+2)]; 
qm = [hr(1:N)   mr(1:N)];
dq = - (ShallowRoe(qr,qp) - ShallowRoe(qm,ql))/h +source;
return