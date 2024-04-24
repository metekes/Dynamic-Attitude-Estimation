function [w_dot] = w_dot_func(J,N,w)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
w_dot=pinv(J)*(N-cross(w,(J*w)));
end