% Driver script for solving the 1D Shallow water equations using a WENO scheme
clear all;
close all;
% Order of schem
m=2;

% Problem parameters
CFL = 0.50;
BC= 'P';
source= false;
L = 2;
FinalTime = 2;
% For the flux, put either 'R' for the Roe flux or 'LF' for the LF flux.
flux = 'LF';

% For fine reference solution
recalc_ref = false;
file_name = 'u_fine3_Weno.mat';
%% Obtain reference solution
if recalc_ref
    % Set fine resolution
    N = 1024;
    h = L/N;
    
    % Initialize solution
    x_fine = [0:h:2]'; %values correspond to cell boundaries
    [h_init, m_init] = initial_3_ex(x_fine);
    
    % Solve Problem
    q = [h_init m_init];
    [q] = ShallowWENO1D(x_fine,q,h,m,CFL,FinalTime,BC, flux, source);
    
    % Save reference solution
    u_fine = q(:,2)./q(:,1);
    save(file_name,'u_fine')
else
    %Load reference solution
    % Set fine resolution
    N = 1024;
    h = L/N;
    % Initialize solution
    x_fine = [0:h:2]';
    load(file_name);
end

%% Check convergence
num_cells = [ 64, 128, 256];
order = [2];
flux='R';
% Initialize solution
error= zeros(length(order),length(num_cells));
for i=1:length(order)
    m= order(i);
    fprintf('Doing order %i. \n', m);
    for j= 1:length(num_cells)
        fprintf('Doing  %i cells. \n', num_cells(j));
        h = L/num_cells(j);
        % Initialize solution
        x = [0:h:2]'; %values correspond to cell boundaries
        [h_init, m_init] = initial_3_ex(x);
        q_c = [h_init m_init];
        [q_c] = ShallowWENO1D(x,q_c,h,m,CFL,FinalTime,BC, flux, source);
        % Calculate the error
        u = q_c(:,2)./q_c(:,1);
        u_ref = interp1(x_fine, u_fine, x);
        diff= u-u_ref;
        error(i,j) = norm(diff,1)*h;
    end
end

%% Do plots
%% Plot the error
% Plot error
figure;
loglog(L./num_cells, error(1,:),'-x', 'LineWidth', 1.5);
hold on;
grid on;
polyfit(log(L./num_cells), log(error(1,:)),1)
%% Show plot for really fine solution.
[h_an, m_an] = initial_1_ex(x-FinalTime);

figure;
plot(x_fine,q(:,1),'b');
hold on;
plot(x_fine,h_an,'--');

figure;
plot(x_fine,q(:,2),'b');
plot(x_fine,m_an,'--');


