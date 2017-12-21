% Driver script for solving the 1D Shallow water equations using a WENO scheme 
clear all

% Order of method
m=2;

% Set problem parameters
L = 2; 
FinalTime = 2; 
N = 256; 
CFL = 0.50; 
h = L/N; 

% Define domain, materials and initial conditions
h_init = zeros(N+1,1); 
m_init = zeros(N+1,1);


% Initialize solution
x = [0:h:2];
[h_init; m_init] = initial_2_ex(x);

% Solve Problem
q = [h_init m_init];
[q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime);