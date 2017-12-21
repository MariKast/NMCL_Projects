function [F] = calcFluxRoePr2(q, periodic, limiter)
%calcFluxLF Calculates the Roe Flux
F = zeros(2, size(q,2)+1);

[ q_plus, q_minus ] = calcqPlusMinus( q, periodic, limiter);
% Calculate exact flux
flux = calcFluxSW(q);

% Transform into z coordinates
trans_factor = 1./sqrt(q(1,:));

z = [trans_factor; trans_factor].*q;

for i=2:size(q,2)
    A_hat = calcAhat(z,i-1,i);
    [S,lambdas] = eig(A_hat);
    lambdas = abs(lambdas);
    A_hat = S* lambdas * inv(S) ;
    F(:,i)=0.5* (flux(:,i-1)+ flux(:,i)) - 0.5*A_hat* (q_plus(:,i)-q_minus(:,i));
end


if periodic
    %Periodic BC
    A_hat = calcAhat(z,1,size(q,2));
    [S,lambdas] = eig(A_hat);
    lambdas = abs(lambdas);
    A_hat = S* lambdas* inv(S);
    F(:,1)= 0.5* (flux(:,1)+ flux(:,end)) - 0.5*A_hat* (q_plus(:,1)-q_minus(:,1));
    F(:,end) = 0.5* (flux(:,1)+ flux(:,end)) - 0.5*A_hat* (q_plus(:,end)-q_minus(:,end));
else
    % Open BC
    F(:,1) = flux(:,1)- 0.5*A_hat* (q_plus(:,1)-q_minus(:,1));
    F(:,end) = flux(:,end)- 0.5*A_hat* (q_plus(:,end)-q_minus(:,end));
    
end