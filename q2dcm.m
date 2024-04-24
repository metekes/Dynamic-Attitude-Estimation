function [A] = q2dcm(q)
% Transfer from quaternion to DCM. The convention for the quaterion is the 4th term is the scaler one.

A(1,1) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
A(1,2) = 2*q(1)*q(2) + 2*q(3)*q(4);
A(1,3) = 2*q(1)*q(3) - 2*q(2)*q(4);
A(2,1) = 2*q(1)*q(2) - 2*q(3)*q(4);
A(2,2) = -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2;
A(2,3) = 2*q(3)*q(2) + 2*q(1)*q(4);
A(3,1) = 2*q(3)*q(1) + 2*q(2)*q(4);
A(3,2) = 2*q(3)*q(2) - 2*q(1)*q(4);
A(3,3) = -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;
end