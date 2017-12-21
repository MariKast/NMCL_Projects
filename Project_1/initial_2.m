function [ h_0, m_0 ] = initial_2(x)
%initial_2 Initial condition of second exercise a)
h_0 = 1- 0.1 * sin(pi*x);
m_0 = zeros(size(h_0));
end

