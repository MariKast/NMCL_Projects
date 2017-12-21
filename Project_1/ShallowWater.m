function [h, m] = ShallowWater(x_mid,h_0,m_0,CFL,final_time, periodic, Roe, source)

% Define all parameters and functions in the right dimensions
time = 0;
q = [h_0; m_0];

u = m_0./h_0;
clear m_0 h_0;
delta_x = x_mid(2)-x_mid(1); 

% Time integration
while ( time < final_time)
    % Determine the time step k
    crit_term = abs(u)+ sqrt(q(1,:));
    lambda_max = max(crit_term);
    k = CFL*delta_x /lambda_max;
    
    % Calculate new time
    if (time + k< final_time)
        time = time+k;
    else %This makes sure we hit the exact final time.
        k = final_time-time;
        time = final_time;
    end
    
    % Calculate the respective numerical fluxes between cells
    if Roe
        num_flux = calcFluxRoe(q, periodic);
    else
        num_flux = calcFluxLF(q,lambda_max, periodic);
    end
    % Conservative scheme for next time step
    q_next = q - k/delta_x*(num_flux(:,2:end)-num_flux(:,1:end-1));
    
    % Evaluate and add source term
    if source
        source_term = calcSourceTerm(x_mid, time-k/2);
    else
        source_term=0;
    end
    q = q_next+source_term*k;
    
    % Update u
    h = q(1,:);
    m = q(2,:);
    u = m./h;
end
h= q(1,:);
m= q(2,:);

end