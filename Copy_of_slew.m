%% Slew manuever with non-0 angular velocities at start and stop 
% this script performs a slew between two attitudes, and has non zero
% angular velocities at both the start and the finish. The trajectory is
% created using Modified Rodrigues Parameters (MRP's) and a cubic spline
% connecting the initial and final conditions. The inverse dynamics are 
% performed on the trajectory, and a feed forward momentum wheel control 
% strategy is planned. 3DOF Gaussian noise is added to both the attitude
% and angular velocity, so the system is linearized at each time step
% before the slew, and Time Variant LQR (TVLQR) (Also known as discrete
% finite horizon LQR) tracks the nonlinear dynamics along the trajectory.
% The Kane/Levinson attitude/quaternion conventions are used. The initial
% and final attitudes are generated randomly each time the script is run.

%% Trajectory Generation
clear
% inputs that describe slew
t_f = 10;                            % s
w_initial = deg2rad([.2 1 .6]);      % rad/s
w_final = deg2rad([1 4 3]);          % rad/s
q_init = randq;                      % N_q_(B_initial)
q_final = randq;                     % N_q_(B_final)
samp_rate = 1000;                    % hz

% convert the error quaternion into Modified Rodrigues Parameters (MRP)
mrp_initial = [0 0 0];
q_init_conj = [-q_init(1:3);q_init(4)];
q_slew = qdot(q_init_conj,q_final);
mrp_final = p_from_q(q_slew)';

% get the MRP dot at start and finish
mrp_d1 = pdot_from_w(mrp_initial,w_initial)';
mrp_d2 = pdot_from_w(mrp_final,w_final)';

% create a cubic spline between the two MRP's
[A,B,C,D,MRP_FX,MRP_DFX,omega_FX,MRP_traj,omega_traj,alpha_traj,...
    tau_ff,rho,rho_dot,rotor,rotor_dot,A_k,B_k,quat_traj] = ...
    cubicspline(mrp_initial,mrp_d1,mrp_final,mrp_d2,t_f,samp_rate,q_init);

%% TVLQR (Finite Horizon Discrete Linear Quadratic Regulator)

% Q and R matrices for TVLQR
Q = .001*eye(9);
Q(1:6,1:6) = 100000*eye(6);
Q(4:6,4:6) = 10*Q(4:6,4:6);
R = .0000000000001*eye(3);
 
% Loop for calculating K gain matrix at each time step during trajectory
N = length(A_k);
S{N} = Q;
K{N} = zeros(3,9);
for k = (N-1):-1:1
   K{k} = inv(R+B_k{k}'*S{k+1}*B_k{k})*B_k{k}'*S{k+1}*A_k{k};
   S{k} = Q + K{k}'*R*K{k} + (A_k{k}-B_k{k}*K{k})'*S{k+1}*(A_k{k} - B_k{k}*K{k});  
end

disp('Generated TVLQR Gains')

%% simulation 

% time vector
t_vec = 0:1/samp_rate:(t_f-(1/samp_rate));
init = [q_init', w_initial, 0 0 0,0 0 0,0 0 0];

% noise covariances
V_quat = .0000001*eye(3);
V_omega = .0000001*eye(3);

disp('begin sim')
%main loop 
for i = 1:length(t_vec)
    
    %store for graphing 
    quat_hist(i,:) = init(1:4);
    omega_hist(i,:) = init(5:7);
    
    % TVLQR
    [init] = controller(init,K{i},quat_traj(:,i),omega_traj(:,i),rotor(:,i),rotor_dot(:,i));
    
    %propagate 1/samp_rate
    [~,y] = ode45(@trajODE,[0,1/samp_rate],init);
   
    %reset initial conditions 
    init = y(end,:)';
      
    %add noise
    [init] = addnoise(init,V_quat,V_omega);
end

disp('done with sim')

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
B = eye(3);


%unpack state 
X = X(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);
rotor = X(11:13);
rotor_dot = X(14:16);

% wheel momentum
rho = B*rotor;
rho_dot = B*rotor_dot;

%dynamics
X_dot = zeros(size(X));
X_dot(1:4) = .5*qdot(quat,[omega;0]);
X_dot(5:7) = -invJ*(rho_dot + cross(omega,J*omega+rho) - tau);
X_dot(11:13) = rotor_dot;

end

function [A,B,C,D,MRP_FX,MRP_DFX,omega_FX,MRP_traj,omega_traj,alpha_traj,tau_ff,rho,rho_dot,rotor,rotor_dot,A_k,B_k,quat_traj] = cubicspline(mrp1,mrp1_dot,mrp2,mrp2_dot,t,samp_rate,q_init)
% Takes MRP and initial and final conditions and creates cubic spline

% x1  - starting MRP
% xd1 - starting MRP dot
% samp_rate - hz


%% get spacecraft inertia properties 

% wheel actuator jacobian
B_w = eye(3);

% moore penrose psuedo inverse of B_w
invB_w = inv(B_w'*B_w)*B_w';

% Spacecraft Inertia
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);

MATA = [3*t^2, 2*t, 1, 0;...
        t^3,  t^2, t, 1 ;...
        0 0 1 0;...
        0 0 0 1];
    
MATB = [mrp2_dot;mrp2;mrp1_dot;mrp1];

ansVec = inv(MATA)*MATB;
A = ansVec(1,:);
B = ansVec(2,:);
C = ansVec(3,:);
D = ansVec(4,:);

MRP_FX= @(t) (A*t^3 + B*t^2 + C*t + D)';
MRP_DFX = @(t) (3*A*t^2 + 2*B*t + C)';
omega_FX = @(t) w_from_pdot(MRP_FX(t),MRP_DFX(t));


w1 = w_from_pdot(mrp1,mrp1_dot);


MRP_traj(:,1) = mrp1(:);
omega_traj(:,1) = w1(:);

rho(:,1) = [0 0 0]';

for i = 1:samp_rate*t
    t_val = (1/samp_rate)*i;
   MRP_traj(:,i+1) = MRP_FX(t_val);
   
   B1_quat_Bt = q_from_p(MRP_traj(:,i+1));
   quat_traj(:,i) = qdot(q_init,B1_quat_Bt);
   
   
   omega_traj(:,i+1) = w_from_pdot(MRP_traj(:,i+1),MRP_DFX(t_val));
   alpha_traj(:,i) = (omega_traj(:,i+1) - omega_traj(:,i))*samp_rate;
   tau_ff(:,i) = J*alpha_traj(:,i) + cross(omega_traj(:,i),J*omega_traj(:,i));
   rho_dot(:,i) = -J*alpha_traj(:,i) - cross(omega_traj(:,i),(J*omega_traj(:,i) + rho(:,i)));
   rho(:,i+1) = rho(:,i) + (1/samp_rate)*rho_dot(:,i);
   
   rotor(:,i) = invB_w*rho(:,i);
   rotor_dot(:,i) = invB_w*rho_dot(:,i);
   
   
   %linearize system about trajectory points 
   w = omega_traj(:,i);
   r = rotor(:,i);
   dt = 1/samp_rate;
   
   A_cont = [-hat(w),      eye(3),                           zeros(3,3);...
             zeros(3,3),  -invJ*(hat(w)*J - hat(J*w + B_w*r)), -invJ*hat(w)*B_w;...
             zeros(3,9)];
   B_cont = [zeros(3,3);...
              -invJ*B_w;...
              eye(3)];
          
   A_k{i} = expm(A_cont*dt);
   B_k{i} = B_cont*dt;
end

%add final point
alpha_traj(:,samp_rate*t+1) = alpha_traj(:,samp_rate*t);
tau_ff(:,samp_rate*t+1) = tau_ff(:,samp_rate*t);
rho_dot(:,samp_rate*t+1) = rho_dot(:,samp_rate*t);
disp('created spline')

end

function [init_out] = controller(init,K,quat_traj,omega_traj,rotor_traj,rotor_dot_traj)

%get spacecraft inertia properties 
J = diag([100 200 300]);
invJ = diag([1/100 1/200 1/300]);
B = eye(3);


%unpack state 
X = init(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);
tau = X(8:10);
rotor = X(11:13);
rotor_dot = X(14:16);

%columnize inputs
quat_traj = quat_traj(:);
omega_traj = omega_traj(:);
rotor_traj = rotor_traj(:);
rotor_dot_traj = rotor_dot_traj(:);

%error phi
planned_dcm = dcm_from_q(quat_traj);
actual_dcm = dcm_from_q(quat);
phi = unhat(logm(planned_dcm'*actual_dcm));
d_omega = omega - omega_traj;
d_rotor = rotor - rotor_traj;

x = [phi;d_omega;d_rotor];

% U = feed-forward term + feedback TVLQR term
U = rotor_dot_traj -K*x;

init_out = init;
init_out(14:16) = U';


end

function [init] = addnoise(init,V_quat,V_omega)
%Add 3D gaussian noise to both quaternion and angular velocity

init = init(:);

X = init(:);
quat = X(1:4)/norm(X(1:4));
omega = X(5:7);

%noise 
phi_noise = mvnrnd([0 0 0],V_quat)';
q_noise = q_from_phi(phi_noise);
w_noise = mvnrnd([0 0 0],V_omega)';

init(1:4) = qdot(q_noise,quat);
init(5:7) = omega + w_noise;

init = init';

end

function [p] = p_from_q(q)
% convert a quaternion, scalar last, to modified rodgrigues parameter

q = q(:);
v = q(1:3);
s = q(4);

p = v/(1+s);

if norm(p)>1
    p = -p/(norm(p)^2);
end
end

function [vechat] = hat(vec)
v1 = vec(1);
v2 = vec(2);
v3 = vec(3);

vechat = [0 -v3  v2    ;...
          v3  0 -v1    ;...
          -v2 v1   0   ];
      
end

function [pdot] = pdot_from_w(p,w)
% this is the kinematics of the modified rodrigues parameter assuming that
% attitude is being denoted as N_R_B using the kane/levinson convention
p = p(:);
w = w(:);

pdot = ((1+norm(p)^2)/4)*(eye(3) + 2*(hat(p)^2 + hat(p))/(1+norm(p)^2))*w;

end

function [q] = q_from_p(p)
% convert modified rodrigues parameter to quaternion scalar last

p = p(:);

q = (1/(1+norm(p)^2))*[2*p ;(1-norm(p)^2)];

q = q(:)/norm(q);

end

function [product] = qdot(q1,q2)
%scalar last quaternion multiplication, hamilton product, qdot

%make column vectors
q1=q1(:);
q2 = q2(:);

%break quats into s and v 
v1 = q1(1:3);
s1 = q1(4);
v2 = q2(1:3);
s2 = q2(4);

%multiply 
product = zeros(4,1);
product(1:3) = s1*v2 + s2*v1 + cross(v1,v2);
product(4) = s1*s2 - v1'*v2;
product = product(:);
end

function [q] = randq

q = rand(4,1);
q = q/norm(q);

end

