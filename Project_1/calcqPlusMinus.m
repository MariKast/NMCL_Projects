function [ q_plus, q_minus ] = calcqPlusMinus( q, periodic, limiter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% For the BC's and domain we need slope values for u_0 to u_N+1,
% that means we need delta values from delta u_0 minus to delta u_N+1 plus

% calculate the delta triangles: u_j_minus is at delta(:,j+1), u_j_plus is at
% delta(:,j+1) ( because we start at j=0 for the BCs)
deltas= zeros(2, size(q,2)+3); % We need two ghost cells for the BC's
deltas(:,3:end-2)= q(:,2:end) -q(:, 1:end-1);
if periodic
    %Periodic BC
    deltas(:,1) = q(:,end) -q(:,end-1);  % delta u_0 minus
    deltas(:,2) = q(:,1) -q(:,end);      % delta u_1 minus
    deltas(:,end-1) = q(:,1)-q(:,end);   % delta u_N plus
    deltas(:,end) = q(:,2)-q(:,1);       % delta u_N+1 plus
    
else
    % Open BC
    deltas(:,1) = 0;
    deltas(:,2) = 0;
    deltas(:,end-1) = 0;
    deltas(:,end) = 0;
end

% We have slope of u_0 until slope N+1
slope= zeros(2, size(q,2)+2);
for j=1:size(slope,2)
    a = deltas(:,j);  % delta u_(j-1) minus
    b = deltas(:,j+1); % delta u_(j-1) plus
    if limiter==1
        slope(:,j) = minmod2(a, b);
    elseif limiter==2
        
        slope(:,j)= minmod3((a+b)/2, 2*a, 2*b);
    else
        slope(:,j)=0;
    end
    
end

% We need q_plus and q_minus for every value of j +/- 0.5
% For a given j, j-0.5 can be found in q_plus(:,j), j+0.5 can be found in
% q_plus(:,j+1), q_plus has the appropriate length for these indices, the
% same goes for indexing in q_minus

% q_minus uses the values from u_0 to u_N
q_minus =zeros(2, size(q,2)+1);
q_minus(:,2:end) = q + slope(:,2:end-1)/2.0;

% q_plus uses the values from u_1 to u_N+1
q_plus= zeros(2, size(q,2)+1);
q_plus(:,1:end-1) = q - slope(:,2:end-1)/2.0;


if periodic
    %Periodic BC
    q_minus(:,1)  = q(:,end) + slope(:,1)/2.0;  % values u_0
    q_plus(:,end) = q(:,1) - slope(:,end)/2.0;  % values u_N+1
    
    
else
    % Open BC
    q_minus(:,1)  = q(:,1) + slope(:,1)/2.0;      % values u_0
    q_plus(:,end) = q(:,end) - slope(:,end)/2.0;  % values u_N+1
end

end

