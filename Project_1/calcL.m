function [ L_term ] = calcL(q, lambda_max, delta_x, x_mid, periodic, Roe, limiter, source,  time )
%calcL Calculate the term L at T= time and with solution q with chosen
%options

% Calculate the respective numerical fluxes between cells
if Roe
    num_flux = calcFluxRoePr2(q, periodic, limiter);
else
    num_flux = calcFluxLFPr2(q,lambda_max, periodic, limiter);
end

% Calculate source term if necessary
if source
    source_term = calcSourceTerm(x_mid, time);
else
    source_term=0;
end

% Assemble the right hand side term
L_term = -1/delta_x*(num_flux(:,2:end)-num_flux(:,1:end-1)) +source_term;

end

