function [delta_theta, delta_bias, P] = EKF_attitude(q0, P_prev, omega, sun_mes_hat, mag_mes_hat, sun_eci_hat, mag_eci_hat, sigma_u, sigma_v)
    dt = 1;
    A = eye(6) + cat(1, cat(2,-skew(omega), -eye(3)), zeros(3,6)) * dt + 1/2*[skew(omega)^2, skew(omega); zeros(3), zeros(3)] * dt;
    B = 0;
    Q = [(sigma_v^2 + 1/3 * sigma_u^2) * eye(3), -1/2 * sigma_u^2 * eye(3);
         -(1/2 * sigma_u^2) * eye(3), sigma_u^2 * eye(3)];
    
    % q = cat(1, 1/2*x_previous(1:3), 1);
    H = cat(1, cat(2, skew(q2dcm(q0)*sun_eci_hat), zeros(3,3)), cat(2, skew(q2dcm(q0)*mag_eci_hat), zeros(3,3)));%+ rand(6,6)*0.0001;
    R = blkdiag(eye(3)*0.002^2, eye(3)*0.015^2);
    
    % prediction
    % G = blkdiag(-eye(3), eye(3));
    P = A*P_prev*A.'+ Q; %%%%%% G*Q*G.'
    
    % correction
    K = P*H.'* inv(H*P*H.'+ R);
    % K = (H*P*H.'+ R)/(P*H.');
    P = (eye(6)-K*H)*P;
    e = (cat(1, sun_mes_hat, mag_mes_hat) - cat(1, q2dcm(q0)*sun_eci_hat, q2dcm(q0)*mag_eci_hat));
    delta_x = K*e;
    
    delta_theta = delta_x(1:3);
    delta_bias = delta_x(4:6);
end