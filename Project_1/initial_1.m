function [ h_0, m_0 ] = initial_1(x)
%initial_1 Initial conditions of first exercise
u_0 = 0.25;
h_0 = func_h_0(x);
m_0 = h_0*u_0;
end

