function [ h_0, m_0 ] = initial_1_ex(x)
%initial_1 Initial conditions of first exercise
h = x(2)-x(1);
x_left = x-0.5*h;
x_right = x+0.5*h;
u_0 = 0.25;
h_0 = 1- 0.5/(pi*h) *( cos(pi*x_right)- cos(pi*x_left)) ;
m_0 = h_0*u_0;
end

