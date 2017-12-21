function [ source_term] = calcSourceTerm(x,t)
%calcSourceTerm Calculates the source term in all the points of vector x
% at time t
u_0 = 0.25;
source_term = zeros(2, length(x));
common_factor = pi/2 * cos(pi*(x-t));
source_term(1,:) = (u_0-1)*common_factor;
source_term(2,:) = common_factor.*( -u_0 +u_0^2 + func_h_0(x-t));
end

