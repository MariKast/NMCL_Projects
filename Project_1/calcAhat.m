function [ A_hat] = calcAhat(z,left, right)
%calcA_hat Function that calculates the linearized Roe matrix for the
%interface between the "left" and "right" cell with the transformed
%coordinates z.

z1_bar = (z(1,left)+ z(1,right))/2;
z2_bar =  (z(2,left)+ z(2,right))/2;
z1_squared_bar = (z(1,left).^2+ z(1,right).^2)/2;

A_hat = [ 0, 1; ...
    (z1_squared_bar-z2_bar^2/z1_bar^2), (2*z2_bar/z1_bar)];
end