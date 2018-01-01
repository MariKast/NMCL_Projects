function [ h_0, m_0 ] = initial_4_ex(x)
%initial_1 Initial conditions of first exercise
h= x(2)-x(1);
x_left = x-0.5*h;
x_right = x+0.5*h;
h_0 = ones(length(x),1);
m_0 = ones(length(x),1)*(-1.5);
for i=1:length(x)
    if (x_left(i) >1)
        m_0(i)=0;
    elseif ( x_right(i) >1)
        m_0(i) = (x_right(i)-1)*(-1.5);
    end
end
end

