function [ h_0, m_0 ] = initial_3_ex(x)
%initial_1 Initial conditions of first exercise
h = x(2)-x(1);
x_left = x-0.5*h;
x_right = x+0.5*h;
h_0 = 1 + 0.1/(pi*h) *( cos(2*pi*x_right)- cos(2*pi*x_left)) ;
m_0 = 0.5* ones(length(x),1);
end

