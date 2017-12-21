% Driver script for solving the 1D Shallow water equations using a WENO scheme 
clear all
close all


% Weno order and number of cells for loop
order =[1];
num_cells = [ 64, 128, 256, 512 1024];

% Set problem parameters
L = 2; 
FinalTime = 2;  
CFL = 0.50; 
source = true;
u_0 = 0.25;
flux ='LF';
BC= 'P';

error= zeros(length(order),length(num_cells));
for i=1:length(order)
    m= order(i);
    fprintf('Doing order %i. \n', m);
    for j= 1:length(num_cells)
        fprintf('Doing  %i cells. \n', num_cells(j));
        h = L/num_cells(j);
        % Initialize solution 
        x = [0:h:2]'; %values correspond to cell boundaries
        x_mid = x(1:end-1)+ 0.5*h;
        [h_init, m_init] = initial_1_ex(x);
        % Solve Problem
        q = [h_init m_init];
        [q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime,BC, flux, source);
        % Calculate the error
        u = q(:,2)./q(:,1);
        diff= u-u_0;
        error(i,j) = norm(diff,1)*h;
    end
end


%% Plot the error
% Plot error
figure;
loglog(L./num_cells, error(1,:),'-x', 'LineWidth', 1.5);
polyfit(log(L./num_cells), log(error(1,:)),1)
%% Show plot for really fine solution.
[h_an, m_an] = initial_1_ex(x-FinalTime);

figure;
plot(x_mid,q(:,1),'b');
hold on;
plot(x_mid,h_an,'--');

figure;
plot(x_mid,q(:,2),'b');
plot(x_mid,m_an,'--');

