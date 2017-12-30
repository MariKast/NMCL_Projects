function [ h_0, m_0 ] = initial_2_ex(x)
%initial_1 Initial conditions of first exercise
h = x(2)-x(1);
x_left= x-0.5*h;
x_right= x+0.5*h;
h_0 = 1 + 0.1/(pi*h) *( cos(pi*x_right)- cos(pi*x_left)) ;
m_0 = zeros(length(x),1);
end

