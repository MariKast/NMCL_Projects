% Driver script for solving the 1D Shallow water equations using a WENO scheme 
clear all
close all
% Order of method
m=3;

% Set problem parameters
L = 2; 
FinalTime = 2; 
N = 512; 
CFL = 0.50; 
h = L/N; 

% Initialize solution 
x = [0:h:2]'; %values correspond to cell boundaries
[h_init, m_init] = initial_1_ex(x);

% Solve Problem
q = [h_init m_init];
[q] = ShallowWENO1D(x,q,h,m,CFL,FinalTime);