function [ h_0, m_0 ] = initial_4(x)
%initial_4 Initial condition of fourth exercise
h_0 = ones(size(x));
m_0 = zeros(size(x));
m_0(x<1)= -1.5;
end

