% Driver script for solving the 1D Shallow water equations using a WENO scheme 
clear all

% Order of method
m=2;

% Set problem parameters
L = 2; 
FinalTime = 0.5; 
N = 256; 
CFL = 0.50; 
h = L/N; 

%%
% Initialize solution 
x = [0:h:2]'; %values correspond to cell boundaries
[h_init, m_init] = initial_4_ex(x);

% Solve Problem
q = [h_init m_init];
BC= 'O';
source= false;

%% 
% For the flux, put either 'R' for the Roe flux or 'LF' for the LF flux.
flux = 'LF';

[q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime,BC, flux, source);