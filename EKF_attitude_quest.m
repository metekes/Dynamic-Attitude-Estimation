function [delta_theta, delta_bias, P] = EKF_attitude_quest(q0, P_prev, omega, sun_mes_hat, mag_mes_hat, sun_eci_hat, mag_eci_hat, sigma_u, sigma_v)
    dt = 1;
    A = eye(6) + cat(1, cat(2,-skew(omega), -eye(3)), zeros(3,6)) * dt + 1/2*[skew(omega)^2, skew(omega); zeros(3), zeros(3)] * dt^2;
    B = 0;
    Q = [(sigma_v^2 + 1/3 * sigma_u^2) * eye(3), -1/2 * sigma_u^2 * eye(3);
         -(1/2 * sigma_u^2) * eye(3), sigma_u^2 * eye(3)];
    
    % q = cat(1, 1/2*x_previous(1:3), 1);
    H = [eye(3), zeros(3)];
    R = eye(3)*0.02^2;
    
    % prediction
    % G = blkdiag(-eye(3), eye(3));
    P = A*P_prev*A.'+ Q; %%%%%% G*Q*G.'
    
    % correction
    K = P*H.'* inv(H*P*H.'+ R);
    % K = (H*P*H.'+ R)/(P*H.');
    P = (eye(6)-K*H)*P;

    q_quest = quest([sun_eci_hat, mag_eci_hat], [sun_mes_hat, mag_mes_hat]);
    delta_q_quest = quat_multipication(q_quest, qInverse(q0));
%     [r, p, y] = dcm2angle(q2dcm(delta_q_quest));
%     delta_euler_quest = [r; p; y];
    delta_euler_quest = [2*delta_q_quest(1); 2*delta_q_quest(2); 2*delta_q_quest(3)]/delta_q_quest(4);

    e = delta_euler_quest;
    delta_x = K*e;
    
    delta_theta = delta_x(1:3);
    delta_bias = delta_x(4:6);
end