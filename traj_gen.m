
%% Time stuff
samp_rate = 100;
T_slew = 10;
t_final = 9.99;
t_vec = 0:(1/samp_rate):t_final;

%% initial conditions 
init_quat = phi_2_quat(rand(3,1));
init_DCM = quat_2_dcm(init_quat);
des_quat = phi_2_quat(rand(3,1));
des_DCM = quat_2_dcm(des_quat);
init_omega = deg2rad([0 0 0]);
init_tau = [0 0 0];
init = [init_quat',init_omega,init_tau];

%spacecraft properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

%% trajectory generator
B0_Q_B = (init_DCM')*des_DCM;      %B0_Q_B
phi = unhat(logm(B0_Q_B)); %ADDED A NEGATIVE HERE. WHY??
r = phi/norm(phi);
theta_max = norm(phi);

%versine function and 1st and 2nd derrivatives
theta_fx = @(t,theta_max,T_slew) theta_max*.5*(1-cos((pi/T_slew)*t));
Dtheta_fx = @(t,theta_max,T_slew) theta_max*.5*(pi/T_slew)*sin((pi/T_slew)*t);
DDtheta_fx = @(t,theta_max,T_slew) theta_max*.5*(pi/T_slew)^2*cos((pi/T_slew)*t);

%time vec
Traj_t_vec = 0:(1/samp_rate):T_slew;

%trajectory sim 
for i = 1:length(Traj_t_vec)
    
    %get time 
    t = Traj_t_vec(i);
    
    %calculate theta and 1st and 2nd derrivatives
    theta = theta_fx(t,theta_max,T_slew);
    Dtheta = Dtheta_fx(t,theta_max,T_slew);
    DDtheta = DDtheta_fx(t,theta_max,T_slew);
    
    %rotation matrix from B initial to current body axes
    Bt_Q_B = expm(hat(r*theta));
    
    %planned trajectory and velocity/accel/tau 
    quat_traj{i} = qdot(init_quat,phi_2_quat(r*theta));%REMOVED A NEGATIVE HERE SO THAT THE TRAJECTORY IS STILL RIGHT
    omega_traj(:,i) = r*Dtheta;
    alpha_traj(:,i) = r*DDtheta;
    tau_ff(:,i) = J*alpha_traj(:,i) + cross(omega_traj(:,i),J*omega_traj(:,i));
    
    
end

%% simulation 

%main loop 
for i = 1:length(t_vec)
    
    %store for graphing 
    quat_hist(i,:) = init(1:4);
    omega_hist(i,:) = init(5:7);
    
    %add feed forward tau
    init(8:10) = tau_ff(:,i)';
    
    %propagate 1/samp_rate
    [~,y] = ode45(@trajODE,[0,1/samp_rate],init);
   
    %reset initial conditions 
    init = y(end,:)';
      
end

%% Plotting 
figure
hold on 
title('Quaternion')
plot(t_vec,quat_hist(:,1));
plot(t_vec,quat_hist(:,2));
plot(t_vec,quat_hist(:,3));
plot(t_vec,quat_hist(:,4));
legend('q_1','q_2','q_3','q_4')
xlabel('Time (s)')
hold off

figure
hold on 
title('Angular Velocity')
plot(t_vec,rad2deg(omega_hist(:,1)));
plot(t_vec,rad2deg(omega_hist(:,2)));
plot(t_vec,rad2deg(omega_hist(:,3)));
legend('\omega_x','\omega_y','\omega_z')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
hold off

quat_hist(end,:)' - des_quat
quat_hist(end,:)' - quat_traj{end}


%% ODE

function [X_dot] = trajODE(t,X)

%get spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

%unpack state 
X = X(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);

%dynamics
X_dot = zeros(size(X));
X_dot(1:4) = .5*qdot(quat,[omega;0]);
X_dot(5:7) = invJ*(tau - cross(omega,J*omega));


end