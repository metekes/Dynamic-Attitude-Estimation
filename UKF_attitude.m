function [delta_euler, delta_bias, P_k] = UKF_attitude(q, P, omega, sun_mes_hat, mag_mes_hat, sun_eci_hat, mag_eci_hat, sigma_u, sigma_v)

    % INPUTS:
    %	 q: quaternion
    %	 P: covariance matrix
    %	 sun_mes_hat, sun_eci_hat: sun vectors
    %	 mag_mes_hat, mag_eci_hat: magnetic filed vectors
    
    % OUTPUTS:
    %	 delta_euler: error in euler angles
    %	 delta_bias: error in biases
    
    % Dynamic model marix
    F = eye(6) + cat(1, cat(2,-skew(omega), -eye(3)), zeros(3,6));
    
    % Dynamic model noise covarinace
    Q = [ (sigma_v^2 + 1/3 * sigma_u^2) * eye(3), -1/2 * sigma_u^2 * eye(3);
         -(1/2 * sigma_u^2) * eye(3), sigma_u^2 * eye(3)];
    
    % Observation model noise covariance
    R = blkdiag(eye(3) * 0.002^2, eye(3) * 0.01^2);
    
    % Define variables to find sigma points
    n      = 6;     % dimension of state space
    kapa   = 5;
    % alpha  = 0.4;  % should be between 0 and 1
    % beta   = 2;    % generally optimum value is 2
    % lambda = alpha^2 * (n + kapa) - n;
    
    % Select sigma points
    X(:,1) = zeros(6,1);  % mean is all 0 vector
    
    % Defining weights
    w_m(1) = kapa/(n+kapa);
    w_c(1) = kapa/(n+kapa);
    
    for i = 1:n
	    w_m(end + 1) = 1 / (2 * (n + kapa));
	    w_m(end + 1) = 1 / (2 * (n + kapa));
	    w_c(end + 1) = 1 / (2 * (n + kapa));
	    w_c(end + 1) = 1 / (2 * (n + kapa));
    end
    
    % Selecting sigma points
    X(:,end+1:end+n) =  chol(((n+kapa) * P), "lower");
    X(:,end+1:end+n) = -chol(((n+kapa) * P), "lower");
    
    % Dynamic update for sigma points
    X = F * X;
    
    % Calculating mean of the sigma points
    X_mean(6,1) = 0;
    
    for i = 1:2*n+1
	    X_mean = X_mean + w_m(i) * X(:,i);
    end
    
    % Calculate dynamic model covariance matrix
    P_ = zeros(6,6);
    
    for i = 1:2*n+1
	    P_ = P_ + w_c(i) * (X(:,i) - X_mean) * (X(:,i) - X_mean)';
    end
    
    P_ = P_ + Q;
    
    % Calculate observations for sigma points
    for i = 1:2*n+1
 	    % Update quaternion using delta euler angles
	    q_ = q + 0.5 * [ q(4), -q(3),  q(2);
		                 q(3),  q(4), -q(1);
		                -q(2),  q(1),  q(4);
   	 	                -q(1), -q(2), -q(3)] * X(1:3,i);
    
	    % Normalize the updated quaternion
   	    q_ = q_ / norm(q_);
    
        % Find observation for sigma points
	    Y(:,i) = cat(1, q2dcm(q_)*sun_eci_hat, q2dcm(q_)*mag_eci_hat);
    end
    
    % Calculate y mean
    Y_mean(6,1) = 0;
    
    for i = 1:2*n+1
	    Y_mean = Y_mean + w_m(i) * Y(:,i);
    end
    
    % Calculate innovation covariance
    P_innov = zeros(6,6);
    
    for i = 1:2*n+1
	    P_innov = P_innov + w_c(i) * (Y(:,i) - Y_mean) * (Y(:,i) - Y_mean)';
    end
    
    P_innov = P_innov + R;
    
    % Calculate cross covarince
    P_cc = zeros(6,6);
    
    for i = 1:2*n+1
	    P_cc = P_cc + w_c(i) * (X(:,i) - X_mean) * (Y(:,i) - Y_mean)';
    end
    
    K = P_cc * inv(P_innov);
    
    % Update quaternion using delta euler angles, which is used to find
    % observation for mean of the sgima points
    q_ = q + 0.5 * [ q(4), -q(3),  q(2);
	                 q(3),  q(4), -q(1);
	                -q(2),  q(1),  q(4);
   	                -q(1), -q(2), -q(3)] * X_mean(1:3);
    
    % Normilze the updated quaternion
    q_ = q_ / norm(q_);
    
    % Calculate the observation of the mean of the sigma points
    Y_x_mean = cat(1, q2dcm(q_)*sun_eci_hat, q2dcm(q_)*mag_eci_hat);
    
    % Calculate state variables
    delta_x = K * ([sun_mes_hat; mag_mes_hat] - Y_mean);
    
    delta_euler = delta_x(1:3);
    delta_bias = delta_x(4:6);
    
    % Propagte the covariance matrix
    P_k = P_ - K*P_innov*K.';

end
