%% Spacecraft Sim
% the below sim uses the Kane/Levinson quaternion convention. Two vector
% basis are used, N which is fixed to the inertial (newtonian) reference 
% frame, and B which is fixed to the spacecraft reference frame. Attitude 
% is stored as parametrizations of the ^N R ^B direction cosine matrix.

% the script 'Slew_manuever_script.m' has the following features that can
% be migrated over to this one easily:
%
%     - Spline trajectory generation between Modified Rodrigues Parameters
%     - Eigen Axis Slew
%     - Time Variant LQR
%     - LQR Station Keeping 
%     - 3D Gausian noise on attitude/angular velocity 

clear

% spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

% timing stuff
samp_rate = 100;                                  % hz
t_f = 20;                                         % s
t_vec = 0:1/samp_rate:(t_f-(1/samp_rate));

% initial conditions 
q_init = randq;
q_init = ([1 2 3 4]/norm([1 2 3 4]))';
w_initial = deg2rad([30 15 25]);                   % rad/s
init = [q_init', w_initial, 0 0 0,0 0 0,0 0 0];

% allocate arrays 
quat_hist = zeros(4,length(t_vec));
omega_hist = zeros(3,length(t_vec));

%main loop 
for i = 1:length(t_vec)
    
    %store for graphing 
    quat_hist(:,i) = init(1:4);
    omega_hist(:,i) = init(5:7);
    
    % controller
    [init] = controller(init,'detumble');
    
    %propagate 1/samp_rate
    [~,y] = ode45(@trajODE,[0,1/samp_rate],init);
   
    %reset initial conditions 
    init = y(end,:)';
      
end

%% Post Processing 

% allocate arrays 
H_B = zeros(3,length(t_vec));
H_N = zeros(3,length(t_vec));

for i = 1:length(t_vec)
    
    % DCM [^N R ^B]
    N_R_B = dcm_from_q(quat_hist(:,i));
    
    % angular momentum expressed in basis b (body fixed)
    H_B(:,i) = J*omega_hist(:,i);
    
    % angular momentum expressed in basis n (inertial)
    H_N(:,i) = N_R_B*H_B(:,i);
    
end


%% Plotting 
figure
hold on 
title('Quaternion')
plot(t_vec,quat_hist(1,:));
plot(t_vec,quat_hist(2,:));
plot(t_vec,quat_hist(3,:));
plot(t_vec,quat_hist(4,:));
legend('q_1','q_2','q_3','q_4')
xlabel('Time (s)')
hold off

figure
hold on 
title('Angular Velocity')
plot(t_vec,rad2deg(omega_hist(1,:)));
plot(t_vec,rad2deg(omega_hist(2,:)));
plot(t_vec,rad2deg(omega_hist(3,:)));
legend('\omega_x','\omega_y','\omega_z')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
hold off

figure
hold on 
title('Angular Momentum (Expressed in Inertial Basis N)')
plot(t_vec,rad2deg(H_N(1,:)));
plot(t_vec,rad2deg(H_N(2,:)));
plot(t_vec,rad2deg(H_N(3,:)));
legend('H nx','H ny','H nz')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
hold off


%% Supporting functions 

function [X_dot] = trajODE(t,X)

% spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

% wheel actuator jacobian
B_w = eye(3);                       

% unpack state 
X = X(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);
rotor = X(11:13);
rotor_dot = X(14:16);

% wheel momentum
rho = B_w*rotor;
rho_dot = B_w*rotor_dot;

% dynamics
X_dot = zeros(size(X));
X_dot(1:4) = .5*qdot(quat,[omega;0]);
X_dot(5:7) = -invJ*(rho_dot + cross(omega,J*omega+rho) - tau);
X_dot(11:13) = rotor_dot;

end

function [init] = controller(init,mode)

% unpack state:
X = init(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);
rotor = X(11:13);
rotor_dot = X(14:16);

% spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

if mode == 'detumble'
    
    % get attitude 
    N_R_B = dcm_from_q(quat);
    
    % angular momentums 
    H_B = J*omega;
    H_N = N_R_B*H_B;
    
    % thrust against the inertial angular momentum
    max_thrust = 10;
    tau_out_N = -max_thrust*(H_N/norm(H_N)); % this is a smoother run down
    %tau_out_N = -max_thrust*sign(H_N);      % this is a more aggressive one
    tau_out_B = (N_R_B')*tau_out_N;
    init(8:10) = tau_out_B';
    
end

end

