function [q_out] = quat_multipication(q1,q2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
q_out = [ q1(4),  q1(3), -q1(2), q1(1);
		  q1(3),  q1(4),  q1(1), q1(2);
		 -q1(2), -q1(1),  q1(4), q1(3);
 	 	 -q1(1), -q1(2), -q1(3), q1(4)] * q2;

q_out = q_out/norm(q_out);
end