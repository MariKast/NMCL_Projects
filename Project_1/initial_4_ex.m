function [ h_0, m_0 ] = initial_4_ex(x)
%initial_1 Initial conditions of first exercise
x_left = x(1:end-1);
x_right = x(2:end);
h_0 = ones(length(x)-1,1);
m_0 = ones(length(x)-1,1)*(-1.5);
for i=1:length(x)-1
    if (x_left(i) >1)
        m_0(i)=0;
    elseif ( x_right(i) >1)
        m_0(i) = (1-x_right(i))*(-1.5);
    end
end
end
