function [result] = inner_product(a,b,n)
result = 0;
for i = 1:n
    result = result+a(i)*b(i);
end
end
