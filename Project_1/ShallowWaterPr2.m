function [h, m] = ShallowWaterPr2(x_mid,h_0,m_0,CFL,final_time, periodic, Roe, source, limiter)

% Define all parameters and functions in the right dimensions
time = 0;
q = [h_0; m_0];

u = m_0./h_0;
clear m_0 h_0;
delta_x = x_mid(2)-x_mid(1); 

% Time integration
while ( time < final_time)
    %% Determine the time step k
    crit_term = abs(u)+ sqrt(q(1,:));
    lambda_max = max(crit_term);
    k = CFL*delta_x /lambda_max;
    
    if (time + k> final_time)
    %This makes sure we hit the exact final time.
        k = final_time-time;
    end
    
    
    %% RK 3 written down explicitely for convenience
    % Stage 1
    L_term = calcL(q, lambda_max, delta_x, x_mid, periodic, Roe, limiter, source, time );
    q_1 = q + k*L_term;
    % Stage 2
    L_term = calcL(q_1, lambda_max, delta_x, x_mid,periodic, Roe, limiter, source, time+k );
    q_2 = 0.75*q +0.25*(q_1 + k*L_term);
    % Stage 3
    L_term = calcL(q_2, lambda_max, delta_x, x_mid, periodic, Roe, limiter, source, time+k/2);
    q_3 = q/3.0 + 2.0/3.0* (q_2 + k*L_term);
    
    % Final/ global stage
    q = q_3;
    time = time+k;
    
    %% Update u
    h = q(1,:);
    m = q(2,:);
    u = m./h;
%     plot(u)
%     hold on;
%     plot(m);
%     plot(h);
%     drawnow
%     hold off;
end
h= q(1,:);
m= q(2,:);

end