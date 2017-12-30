function [q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime, BC, flux, source)
% function [q] = ShallowWENODriver1D(x,q,h,m,CFL,gamma,FinalTime)
% Purpose  : Integrate 1D Shallow water equations until FinalTime using a WENO
%            scheme and a 3rd order SSP-RK
time = 0; tstep = 0;

% Initialize reconstruction weights
Crec = zeros(m+1,m);
for r=-1:m-1
    Crec(r+2,:) = ReconstructWeights(m,r);
end;

% Initialize linear weights
dw = LinearWeights(m,0);

% Compute smoothness indicator matrices
beta = zeros(m,m,m);
for r=0:m-1
    xl = -1/2 + [-r:1:m-r];
    beta(:,:,r+1) = betarcalc(xl,m);
end

% integrate scheme
while (time<FinalTime)  
%% Determine the time step k
    crit_term = abs(q(:,2)./q(:,1))+ sqrt(q(:,1));
    lambda_max = max(crit_term);
    k = CFL*h /lambda_max;
    
    if (time + k> FinalTime)
    %This makes sure we hit the exact final time.
        k = FinalTime-time;
    end

    
  %% Update solution
  rhsq  = ShallowWENOrhs1D(x,q,h,time,k,m,Crec,dw,beta, lambda_max, BC, flux, source); 
  q1 = q + k*rhsq;
  rhsq  = ShallowWENOrhs1D(x,q1,h,time+k,k,m,Crec,dw,beta,lambda_max, BC, flux, source); 
  q2 = (3*q + q1 + k*rhsq)/4;
  rhsq  = ShallowWENOrhs1D(x,q2,h,time+k/2,k,m,Crec,dw,beta,lambda_max, BC, flux, source); 
  q  = (q + 2*q2 + 2*k*rhsq)/3;
  time = time+k; 
  tstep = tstep+1;

end
return