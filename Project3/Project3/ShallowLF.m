function [numflux] = ShallowLF(u,v,lambdamax)
% function [numflux] = EulerLF(u,v,gamma,lambda,maxvel);
% Purpose: Evaluate global Lax Friedrich numerical flux for 
% the Euler equations

% Compute flux for u
fu = calcFluxSW(u);

% Compute flux for v
fv = calcFluxSW(v);

% Evaluate numerical flux
numflux = (fu+fv)/2 - lambdamax/2*(v-u);
return