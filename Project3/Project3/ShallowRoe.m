function [numflux] = ShallowRoe(u,v)
% function [numflux] = EulerLF(u,v,gamma,lambda,maxvel);
% Purpose: Evaluate global Lax Friedrich numerical flux for 
% the Euler equations

% Compute flux for u
fu = calcFluxSW3(u);

% Compute flux for v
fv = calcFluxSW3(v);

% Transform into z coordinates
fact_u = 1./sqrt(u(:,1));
z_u = [fact_u fact_u].*u;
fact_v = 1./sqrt(v(:,1));
z_v = [fact_v fact_v].*v;

F = zeros( size(fu));
for i=1:size(u,1)
    A_hat = calcAhatuv(z_u(i,:),z_v(i,:));
    [S,lambdas] = eig(A_hat);
    lambdas = abs(lambdas);
    A_hat = S* lambdas * inv(S) ;
    F(i,:)= (0.5*A_hat* (u(i,:)-v(i,:))');
end

% Evaluate numerical flux
numflux = (fu+fv)/2 + F; 
return