function [F] = calcFluxLFPr2(q, lambda_max, periodic, limiter)
%calcFluxLF Calculates the LF Flux

[ q_plus, q_minus ] = calcqPlusMinus( q, periodic, limiter);
% We need flux values from j=1: F_-0.5 to j=N: F_N+0.5
% The flux values needed for u_j= ..- (F_j+0.5 -F_j-0.5) are located for
% j-0.5 at  F(:,j) and j+0.5 at F(:,j+1)

% For the flux, we need the q_minus q_plus values of j+/- 0.5, we have
% computed them above, including BC's. However, we also need the actual
% flux of the problem, so we still impose BC's below

F= zeros(2, size(q,2)+1);
% Calculate exact flux
flux = calcFluxSW(q);
% LF flux for inner values
F(:,2:end-1)=0.5* (flux(:,1:end-1)+ flux(:,2:end)) - 0.5*lambda_max * (q_plus(:,2:end-1)-q_minus(:,2:end-1));

% BC's
if periodic
    %Periodic BC
    F(:,1) = 0.5* (flux(:,end)+ flux(:,1)) - 0.5*lambda_max * (q_plus(:,1)-q_minus(:,1));
    F(:,end) = 0.5* (flux(:,end)+ flux(:,1)) - 0.5*lambda_max * (q_plus(:,end)-q_minus(:,end));
else
    % Open BC
    F(:,1) = flux(:,1)  - 0.5*lambda_max * (q_plus(:,1)-q_minus(:,1));
    F(:,end) = flux(:,end) - 0.5*lambda_max * (q_plus(:,end)-q_minus(:,end));
    
end