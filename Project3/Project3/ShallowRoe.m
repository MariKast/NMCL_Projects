function [numflux] = ShallowRoe(u,v,q)
% function [numflux] = EulerLF(u,v,gamma,lambda,maxvel);
% Purpose: Evaluate global Lax Friedrich numerical flux for 
% the Euler equations

% Compute flux for u
fu = calcFluxSW(u);

% Compute flux for v
fv = calcFluxSW(v);

% Transform into z coordinates
trans_factor = 1./sqrt(q(1,:));

z = [trans_factor; trans_factor].*q;


% Evaluate numerical flux
numflux = (fu+fv)/2 - lambdamax/2*(v-u);
return