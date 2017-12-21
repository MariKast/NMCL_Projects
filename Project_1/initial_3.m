function [ h_0, m_0 ] = initial_3(x)
%initial_3 Initial condition of second exercise b)
h_0 = 1- 0.2 * sin(2*pi*x);
m_0 = ones(size(h_0))*0.5;
end

