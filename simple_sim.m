%% Spacecraft Sim
clear

% timing stuff
samp_rate = 100;                                  % hz
t_f = 10;                                         % s
t_vec = 0:1/samp_rate:(t_f-(1/samp_rate));

% initial conditions 
q_init = randq;
w_initial = deg2rad([5 10 15]);                   % rad/s
init = [q_init', w_initial, 0 0 0,0 0 0,0 0 0];

% allocate arrays 
quat_hist = zeros(length(t_vec),4);
omega_hist = zeros(length(t_vec),3);

%main loop 
for i = 1:length(t_vec)
    
    %store for graphing 
    quat_hist(i,:) = init(1:4);
    omega_hist(i,:) = init(5:7);
    
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

%% supporting functions 

function [X_dot] = trajODE(t,X)

%get spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

% wheel actuator jacobian
B_w = eye(3);                       

%unpack state 
X = X(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);
rotor = X(11:13);
rotor_dot = X(14:16);

% wheel momentum
rho = B_w*rotor;
rho_dot = B_w*rotor_dot;

%dynamics
X_dot = zeros(size(X));
X_dot(1:4) = .5*qdot(quat,[omega;0]);
X_dot(5:7) = -invJ*(rho_dot + cross(omega,J*omega+rho) - tau);
X_dot(11:13) = rotor_dot;

end


