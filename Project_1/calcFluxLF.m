function [F] = calcFluxLF(q, lambda_max, periodic)
%calcFluxLF Calculates the LF Flux
F= zeros(2, size(q,2)+1);
% Calculate exact flux
flux = calcFluxSW(q);
% LF flux
F(:,2:end-1)=0.5* (flux(:,1:end-1)+ flux(:,2:end)) - 0.5*lambda_max * (q(:,2:end)-q(:,1:end-1));

if periodic
    %Periodic BC
    F(:,1) = 0.5* (flux(:,end)+ flux(:,1)) - 0.5*lambda_max * (q(:,1)-q(:,end));
    F(:,end) = F(:,1);
else
    % Open BC
    F(:,1) = flux(:,1);
    F(:,end) = flux(:,end);
    
end