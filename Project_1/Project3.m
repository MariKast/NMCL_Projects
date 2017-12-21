% Driver script for solving the 1D Burgers equation using a monotone scheme
clear all
close all;
plot_it = true;
cell_sizes =[ 16, 32, 64, 128, 256  512]; %, 2048, 4096];
cell_sizes= [512];
errors = zeros(size(cell_sizes));
% Set problem parameters
domain_size = 2;
final_time = 0.5;
CFL = 0.50;
periodic = false;
Roe = true;
source = false;
limiter= 2;
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
[h, m] = ShallowWaterPr3(x_mid,h_0,m_0,CFL,final_time, periodic, Roe, source, limiter);

% analytical solution

h_analyt = func_h_0(x_mid-final_time);
m_analyt = u_0*h_analyt;
u_analyt = m_analyt./h_analyt; 

u= m./h;

error = norm((u-u_analyt),1)*delta_x;
errors(iter)= error;
end
%% Plot things
if true
%load('error1Roe.mat');
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

plot(x_mid, h,'-.', 'LineWidth', 1.5)
plot(x_mid, m,':', 'LineWidth', 1.5)
%plot(x_mid, m_analyt,'-')
plot(x_mid, u,'--', 'LineWidth', 1.5)
xlabel('x')


legend('h(x)', 'm(x)','u(x)')

end