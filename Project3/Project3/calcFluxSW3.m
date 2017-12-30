function [ flux] = calcFluxSW3( q )
%calcFluxSW  Calculates the exact flux of the shallow water equation.
flux = zeros(size(q));
h = q(:,1);
m = q(:,2);
flux(:,1) = m;
flux(:,2) = m.^2./h + 0.5*h.^2; 
end

