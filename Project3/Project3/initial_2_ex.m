function [ h_0, m_0 ] = initial_2_ex(x)
%initial_1 Initial conditions of first exercise
x_left = x(1:end-1);
x_right = x(2:end);
h = x(2)-x(1);
h_0 = 1 + 0.1/(pi*h) *( cos(pi*x_right)- cos(pi*x_left)) ;
m_0 = zeros(length(x)-1,1);
end

