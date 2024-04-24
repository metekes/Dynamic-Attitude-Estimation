clc; clearvars -except delta_q_quest q;
load("mes_nanosat_15082018.mat","sun_eci","mag_eci");

mag_eci=mag_eci.';
J=[0.037, 0, 0; 0, 0.037, 0; 0, 0, 0.01];
q0=[0; 0; 0; 1];
w0=[0.001; -0.006; 0.002] * pi/180 * 10;
N=[10^-10; 5*10^-7; 0]*0;
b0=[0; 0; 0];

dt = 1;
t_max = 7200;

sigma_u = 6.65E-6;
sigma_v = 2.32E-5;

for i=1:dt:t_max
    q1_dot=q_dot_func(w0,q0);
    w1_dot=w_dot_func(J,N,w0);

    q2_dot=q_dot_func(w0+dt/2*w1_dot,q0+dt/2*q1_dot);
    w2_dot=w_dot_func(J,N,w0+dt/2*w1_dot);

    q3_dot=q_dot_func(w0+dt/2*w2_dot,q0+dt/2*q2_dot);
    w3_dot=w_dot_func(J,N,w0+dt/2*w2_dot);

    q4_dot=q_dot_func(w0+dt*w3_dot,q0+dt*q3_dot);
    w4_dot=w_dot_func(J,N,w0+dt*w3_dot);

    q_dot=1/6*(q1_dot+2*q2_dot+2*q3_dot+q4_dot);
    w_dot=1/6*(w1_dot+2*w2_dot+2*w3_dot+w4_dot);

    q0=q0+dt*q_dot;
    w0=w0+dt*w_dot;

    q0=q0/norm(q0);
    
    Wx_ac(:,i)=w0(1);
    Wy_ac(:,i)=w0(2);
    Wz_ac(:,i)=w0(3);
    sun_mes(:,i)=q2dcm(q0)*sun_eci(:,i)+randn(3,1)*0.002;
    mag_mes(:,i)=q2dcm(q0)*mag_eci(:,i)+randn(3,1)*300;
    
    b_dot=randn(3,1)* sigma_u;
    b0 = b0 + dt * b_dot;
    gauss_noise = randn(3,1) * sigma_v;

    Wx(:,i) = Wx_ac(:,i) + b0(1) + gauss_noise(1);
    Wy(:,i) = Wy_ac(:,i) + b0(2) + gauss_noise(2);
    Wz(:,i) = Wz_ac(:,i) + b0(3) + gauss_noise(3);

    % omega(:,i) = cat(1,Wx(:,i),Wy(:,i),Wz(:,i));

    q_real(:,i) = q0;

end

clear("q0")
save("my_data.mat")


