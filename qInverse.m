function [inverse] = qInverse(q)
inverse = [-q(1);-q(2);-q(3);q(4)]./norm(q)^2;
end