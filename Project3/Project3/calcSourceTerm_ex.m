function [ source_term] = calcSourceTerm_ex(x,t)
%calcSourceTerm Calculates the source term in all the points of vector x
% at time t
h = x(2)-x(1);
x_left = x(1:end-1);
x_right = x(2:end);
u_0 = 0.25;
g = 1;

source_term = zeros(length(x)-1,1);
arg_right = pi* (x_right-t);
arg_left =  pi* (x_left-t);
source_term(:,1) = (u_0-1)/(h*2)* (sin(arg_right) -sin(arg_left) );
source_term(:,2) = (-u_0+u_0^2+g) /(h*2)*(sin(arg_right) -sin(arg_left) ) - ...
                    (cos(arg_right).^2 -cos(arg_left).^2)/ (8.0*h) ;
end

