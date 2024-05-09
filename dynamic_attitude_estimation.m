% Load data
% VARIABLES:
%	sun_mes_hat, sun_eci_hat: Sun vectors
%	mag_mes_hat, mag_eci_hat: Magnetic filed vectors
% 	Wx, Wy, Wz: Angular velocities in x, y, z
%   q_real: Generated attitude

run("data_generate_v1.m")
clc; clear; close all;
load("my_data.mat")

% Define constants
t_max = 7200;       % max time
dt = 1;		        % time step
bias = [0; 0; 0];   % initial bias
omega = [0; 0; 0];  % initial omega
P = blkdiag(eye(3) * 0.01, eye(3) * 10^-4); 	                   % initial covarince matrix
initial_q = [1*sin(pi/6); 2*sin(pi/6); 3*sin(pi/6); cos(pi/6)];    % inital quaternion
initial_q = initial_q / norm(initial_q);                           % normalize quaternion

% quaternion input of filter for the first itter
q = initial_q;

% Normilze the measurement vectors
for i = 1:t_max
    sun_eci_hat(:,i) = sun_eci(:,i)/norm(sun_eci(:,i));
    sun_mes_hat(:,i) = sun_mes(:,i)/norm(sun_mes(:,i));
    mag_eci_hat(:,i) = mag_eci(:,i)/norm(mag_eci(:,i));
    mag_mes_hat(:,i) = mag_mes(:,i)/norm(mag_mes(:,i));
end

% Start the loop
for i = 1:t_max

	% Call the UKF filter
	[delta_euler, delta_bias, P_k] = EKF_attitude_quest(q, P, omega, sun_mes_hat(:,i), mag_mes_hat(:,i), sun_eci_hat(:,i), mag_eci_hat(:,i), sigma_u, sigma_v);

	% Update P matrix
	P = P_k;

    % Store covariance history
    covar_history(:,:,i) = P;

	% Update quaternion using delta euler angles
	q = q + 0.5 * [ q(4), -q(3),  q(2);
		            q(3),  q(4), -q(1);
		           -q(2),  q(1),  q(4);
   	 	           -q(1), -q(2), -q(3)] * delta_euler;

	% Normilze the updated quaternion
   	q = q / norm(q);

	% Update the bias
	bias = bias + delta_bias;

    % Store bias history
    delta_bias_history(:,i) = bias - bias_real(:,i);

	% Update the omegas obtained form gyro using biases
	omega = [Wx(:,i); Wy(:,i); Wz(:,i)];
	omega = omega - bias;

	% Update quaternion in time ( t --> t+1 )
% 	q_dot = 0.5 * [ q(4), -q(3),  q(2);
% 			        q(3),  q(4), -q(1);
% 		           -q(2),  q(1),  q(4);
%  	 	           -q(1), -q(2), -q(3)] * omega;

    q_dot = 0.5 * [ 0,         omega(3), -omega(2), omega(1);
                   -omega(3),  0,         omega(1), omega(2);
                    omega(2), -omega(1),  0,        omega(3);
                   -omega(1), -omega(2), -omega(3), 0] * q;

	q = q + q_dot * dt;
    q = q / norm(q);

    % Save the euler difference
    delta_q_history(:,i) = cat(2, [q(4),  q(3), -q(2);
                                  -q(3),  q(4),  q(1);
                                   q(2), -q(1),  q(4);
                                  -q(1), -q(2), -q(3)], q) * qInverse(q_real(:,i));
end

figure(1)
subplot(3,1,1)
hold on
grid on
plot(2 * delta_q_history(1,:));
plot(reshape(6 * sqrt(covar_history(1,1,:)), 1, t_max))
plot(reshape(-6 * sqrt(covar_history(1,1,:)), 1, t_max))
title("Roll Error")
xlabel("t [s]")
ylabel("roll error [rad]")
legend("Error", "3\sigma", "-3\sigma")
% ylim([-0.15, 0.15])

subplot(3,1,2)
hold on
grid on
plot(2 * delta_q_history(2,:));
plot(reshape(6 * sqrt(covar_history(2,2,:)), 1, t_max))
plot(reshape(-6 * sqrt(covar_history(2,2,:)), 1, t_max))
title("Pitch Error")
xlabel("t [s]")
ylabel("pitch error [rad]")
legend("Error", "3\sigma", "-3\sigma")
% ylim([-0.15, 0.15])

subplot(3,1,3)
hold on
grid on
plot(2 * delta_q_history(3,:));
plot(reshape(6 * sqrt(covar_history(3,3,:)), 1, t_max))
plot(reshape(-6 * sqrt(covar_history(3,3,:)), 1, t_max))
title("Yaw Error")
xlabel("t [s]")
ylabel("yaw error [rad]")
legend("Error", "3\sigma", "-3\sigma")
% ylim([-0.15, 0.15])
exportgraphics(gcf,'attitude_error_quest.jpg', BackgroundColor='none',ContentType='image')

figure(2)
subplot(3,1,1)
hold on
grid on
plot(delta_bias_history(1,:));
plot(reshape(3 * sqrt(covar_history(4,4,:)), 1, t_max))
plot(reshape(-3 * sqrt(covar_history(4,4,:)), 1, t_max))
title("b_x Error")
xlabel("t [s]")
ylabel("b_x error [rad/s]")
legend("Error", "3\sigma", "-3\sigma")
ylim([-2E-4, 2E-4])

subplot(3,1,2)
hold on
grid on
plot(delta_bias_history(2,:));
plot(reshape(3 * sqrt(covar_history(5,5,:)), 1, t_max))
plot(reshape(-3 * sqrt(covar_history(5,5,:)), 1, t_max))
title("b_y Error")
xlabel("t [s]")
ylabel("b_y error [rad/s]")
legend("Error", "3\sigma", "-3\sigma")
ylim([-2E-4, 2E-4])

subplot(3,1,3)
hold on
grid on
plot(delta_bias_history(3,:));
plot(reshape(3 * sqrt(covar_history(6,6,:)), 1, t_max))
plot(reshape(-3 * sqrt(covar_history(6,6,:)), 1, t_max))
title("b_z Error")
xlabel("t [s]")
ylabel("b_z error [rad/s]")
legend("Error", "3\sigma", "-3\sigma")
ylim([-2E-4, 2E-4])
exportgraphics(gcf,'bias_error_quest.jpg', BackgroundColor='none',ContentType='image')