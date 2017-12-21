% Driver script for solving the 1D Burgers equation using a monotone scheme
clear all
close all;
plot_it = true;
%cell_sizes =[ 16, 32, 64, 128, 256  512, 1024]; %, 2048, 4096];
cell_sizes = [1024];
errors = zeros(size(cell_sizes));
% Set problem parameters
domain_size = 2;
final_time = 0.5;
CFL = 0.50;
periodic = false;
Roe = false;
source = false;

%% Ref sol

num_cells = 10000;
delta_x = domain_size/num_cells;

% Define domain and initial conditions
cell_bounds = (0:delta_x:domain_size);
x_mid_fine = cell_bounds(1:end-1)+0.5*delta_x;
load('u_fine4.mat');
% %Midpoint rule
% [h_0, m_0] = initial_4(x_mid_fine);
% %  [h_0, m_0] = initial_4(x_mid_fine);
% 
% %Solve Problem
%  [h_fine, m_fine] = ShallowWater(x_mid_fine,h_0,m_0,CFL,final_time, periodic, Roe, source);
% 
% % analytical solution
%  
% 
% u_fine= m_fine./h_fine;
%file_name = 'u_fine3.mat';
%save(file_name,'u_fine')
%% Conv plot

for iter=1:length(cell_sizes)
num_cells = cell_sizes(iter);
delta_x = domain_size/num_cells;

% Define domain and initial conditions
cell_bounds = (0:delta_x:domain_size);
x_mid = cell_bounds(1:end-1)+0.5*delta_x;
u_0 = 0.25;
% Midpoint rule
[h_0, m_0] = initial_4(x_mid);

% Solve Problem
[h, m] = ShallowWater(x_mid,h_0,m_0,CFL,final_time, periodic, Roe, source);

% analytical solution
u_ref = interp1(x_mid_fine, u_fine, x_mid);
u= m./h;

error = norm((u-u_ref),1)*delta_x;
errors(iter)= error;
end
%% Plot things
if false
%load('error3Roe.mat');

figure('DefaultAxesFontSize',14)
loglog(domain_size./cell_sizes, errors,'-x', 'LineWidth', 1.5);
hold on;
%loglog(domain_size./cell_sizes, errorsRoe,'-.x', 'LineWidth', 1.5);
xlabel('\Delta x');
ylabel( 'error ||u-u^*||_1');
polyfit(log(cell_sizes), log(errors),1)
grid on;
legend('LF','Roe')

end

if plot_it
figure('DefaultAxesFontSize',14)
%plot(x_mid, h_analyt,'-')
grid on;
set(gca, 'XMinorGrid','on', 'YMinorGrid','on')
hold on;

plot(x_mid, h,'.', 'LineWidth', 1.5)
plot(x_mid, m,'o', 'LineWidth', 1.5)
%plot(x_mid, m_analyt,'-')
plot(x_mid, u,'x', 'LineWidth', 1.5)
xlabel('x')
%axis([ 0 2 0.2 1.2])

legend('h(x)', 'm(x)','u(x)')

end