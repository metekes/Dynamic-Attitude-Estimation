function [output] = skew (input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = [0, -input(3), input(2); input(3), 0, -input(1); -input(2), input(1), 0];
end