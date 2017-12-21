function [ Lterm ] = calcLterm(q, lambda_max, delta_x, x_mid, Roe, limiter, time )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Calculate the respective numerical fluxes between cells
if Roe
    num_flux = calcFluxRoe(q, periodic);
else
    num_flux = calcFluxLF(q,lambda_max, periodic);
end



end

