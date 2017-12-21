% Driver script for solving the 1D Shallow water equations using a WENO scheme 
clear all



recalc_ref = true;
%% Obtain reference solution 
if recalc_ref
% Order of method
m=3;
file_name = 'u_fine_Weno.mat';
% Set problem parameters
L = 2; 
FinalTime = 2; 
N = 2048; 
CFL = 0.50; 
h = L/N; 

% Initialize solution 
x = [0:h:2]'; %values correspond to cell boundaries
[h_init, m_init] = initial_1_ex(x);

% Solve Problem
q = [h_init m_init];
BC= 'P';
% For the flux, put either 'R' for the Roe flux or 'LF' for the LF flux.
flux = 'LF';
source= true;
[q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime,BC, flux, source);
u_fine = q(:,2)/q(:,1);
save(file_name,'u_fine')
end
% Order of method
m=2;

% Set problem parameters
L = 2; 
FinalTime = 2; 
N = 256; 
CFL = 0.50; 
h = L/N; 
%%
% Initialize solution 
x = [0:h:2]'; %values correspond to cell boundaries
[h_init, m_init] = initial_2_ex(x);

% Solve Problem
q = [h_init m_init];
BC= 'P';
source= false;

%% 
% For the flux, put either 'R' for the Roe flux or 'LF' for the LF flux.
flux = 'R';

[q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime,BC, flux, source);