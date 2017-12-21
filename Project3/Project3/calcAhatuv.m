function [ A_hat] = calcAhatuv(z_u,z_v)
%calcA_hat Function that calculates the linearized Roe matrix for the
%interface between the "left" and "right" cell with the transformed
%coordinates z.

z1_bar = (z_u(1)+ z_v(1))/2;
z2_bar =  (z_u(2)+ z_v(2))/2;
z1_squared_bar = (z_u(1).^2+ z_v(1).^2)/2;

A_hat = [ 0, 1; ...
    (z1_squared_bar-z2_bar^2/z1_bar^2), (2*z2_bar/z1_bar)];
end