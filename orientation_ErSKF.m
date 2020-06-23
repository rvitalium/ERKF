%{
Error-State Kalman Filter

Nominal state kinematics: deterministic, nonlinear, corrected with error states

Nominal States (10 total)
q: orietation (quaternion formulation) (4x1)
wb: angular rate gyroscope bias (3x1)

Error state kinematics: stochastic, linear

Error States (6 total)
dth: quaternion error state (theta) (3x1)
dwb: angular rate gyroscope bias error state (3x1)

Process Noise Covariance Matrix (Q)
Ti: covariance matrix for random impulse applied to quaternion error state
Oi: covariance matrix for random impulse applied to angular rate gyroscope bias error state

Measurement Noise Covariance Matrix (R)
sig_an: covariance for accelerometer measurements
sig_mn: covariance for magnetometer measurements

%}

clc;
clear;
close all;
warning('off');

%% load file

switch 1
    case 1 % hardcode file name
        load('validation1.mat');
    case 2 % chose file from folder
        home = cd;
        [filename,pathname]=uigetfile('*.mat','file');
        cd(pathname);
        load(filename);
        cd(home);
end

seg = fieldnames(IMU);
seg_ind = find(~strcmp(seg,'time')&~strcmp(seg,'button'));

%% initialize IMU data

time = IMU.time;
dt = mean(diff(time)); % average sampling period in [s]
g = -[0; 0; 9.81]; % gravity vector in world frame

accel1 = IMU.(seg{seg_ind}).a; % IMU acceleration
gyro1 = IMU.(seg{seg_ind}).w; % IMU angular velocity
mag1 = IMU.(seg{seg_ind}).m; % IMU magnetic field
u1 = [accel1'; gyro1'; mag1']; % input to the system

% indices
u_ind.a = 1:3; % indices for acceleration
u_ind.w = 4:6; % indices for angular velocity
u_ind.m = 7:9; % indices for magnetic field

%% define noise matrices

% sensor noise
sig_an = (120e-6)*9.81*sqrt(50); % accelerometer noise in [m/s^2]
sig_wn = 0.025*(pi/180)*sqrt(50); % angular rate gyro noise in [rad/s]
sig_mn = (1.8e-3)*sqrt(32.5); % magnetometer noise in [Gauss]

% sensor bias noise
sig_ww = 0.025*(pi/180)*50; % angular rate gyro bias noise in [rad/s * sqrt(s)]

% random impulse noise for error state kinematics
Ti = (sig_wn.^2)*(dt^2).*eye(3); % [rad^2]
Oi = (2.25e-4.^2)*dt.*eye(3); % [rad^2/s^2]

% define process noise covariance matrix (Q)
Q = blkdiag(Ti,Oi); % process noise covariance matrix

% define measurement noise covariance matrix (R)
R = [(sig_an^2).*eye(3), zeros(3); zeros(3), (sig_mn^2).*eye(3)]; % accelerometer and magnetometer noise

%% initialize nominal states (x)

% initial orientation using first accelerometer and magnetometer measurements
Z = accel1(1,:)./norm(accel1(1,:)); % direction of gravity (world Z-axis)
X = mag1(1,:)./norm(mag1(1,:)); % direction of magnetic north (world X-axis)
Y = cross(Z,X); % cross product to get world Y-axis
Y = Y./norm(Y); % normalized the Y-axis
X = cross(Y,Z); % make sure axes are orthonomal (magnetic field has an inclination angle)
dcm = [X;Y;Z]; % rows are world frame axes resolved in senor frame
q = rotm2quat(dcm)'; % convert to quaternion

% initial angular rate gyro bias using static period for entire trial
I = sqrt(sum(accel1.^2,2)) - norm(g) < 0.25 & sqrt(sum(gyro1.^2,2)) < 0.05;
wb = mean(gyro1(I,:))'; % averages for angular rate gyroscope bias
sig_wb = std(sqrt(sum(gyro1(I,:).^2,2)));

x1 = [q; wb]; % initial nominal states

% inclination angle
I = find(I);
for i = 1:length(I)
    M(i) = acos((accel1(I(i),:)*(mag1(I(i),:)'))/(norm(accel1(I(i),:))*norm(mag1(I(i),:))));
end
psi1 = mean(M);
sigma_m = std(M);

% magnetic field strength
M_mag = sqrt(sum(mag1(1:128*3,:).^2,2));
mu_mag = mean(M_mag);
sigma_mag = std(M_mag);

% indices
x_ind.q = 1:4; % index for quaternion
x_ind.wb = 5:7; % index for angular rate gyro bias

%% initialize error states (dx)

dth = zeros(3,1); % dth(eta) instead of dq
dwb = zeros(3,1); % angular rate gyro bias

dx1 = [dth; dwb]; % initial error states

% indices
dx_ind.dth = 1:3; % index for quaternion error state
dx_ind.dwb = 4:6; % index for angular rate gyro bias error state

%% initialize uncertainty covariance matrix (P)

G = [0;0;1];
M = [1;0;0];
w = q(1);
v = q(2:4);
h1 = 2*[w.*G + cross(G,v), v'*G*eye(3) + v*G' - G*v' + w.*skew(G)];
h2 = 2*[w.*M + cross(M,v), v'*M*eye(3) + v*M' - M*v' + w.*skew(M)];
Hx = [h1; h2];
Hx = [Hx, zeros(6,3)];
Q_dth = 0.5*[-q(2), -q(3), -q(4);
              q(1), -q(4),  q(3);
              q(4),  q(1), -q(2);
             -q(3),  q(2),  q(1)];
X_dth = blkdiag(Q_dth,eye(3));
H = Hx*X_dth;

P1 = Q + (H')*R*H; % process noise plus measurement noise
P_hat1 = P1;

%% main loop

ka = round(0.025/dt);
km = round(0.075/dt);

for i = 1:length(time)
    
    % step 1: predict nominal state kinematics
    x1_hat(:,i+1) = nomStateKin(x1(:,i),u1(:,i),dt,x_ind,u_ind);

    % step 2: predict error state covariance
    P1_hat(:,:,i+1) = errStateKin(x1_hat(:,i+1),u1(:,i),dt,P1(:,:,i),Q,x_ind,u_ind);

    % step 3: determine if IMU is static -> accelerometer data usable for an update
    if i <= ka
        static(i) = motionDynamics(x1_hat(:,2:i+1),u1(:,1:i),g,sig_wb,x_ind,u_ind);
    else
        static(i) = motionDynamics(x1_hat(:,i-ka+1:i+1),u1(:,i-ka:i),g,sig_wb,x_ind,u_ind);
    end
    K.static = static(i);

    % step 4: determine if IMU magnetic field is normal -> magnetometer data usable for an update
    if i <= km && i >= 2
        magnet(i) = magneticField(static(i),psi1,sigma_m,mu_mag,sigma_mag,x1(:,1:i),x1_hat(:,2:i+1),u1(:,1:i),x_ind,u_ind,dt);
    elseif i <= km && i == 1
        magnet(i) = 0;
    else
        magnet(i) = magneticField(static(i),psi1,sigma_m,mu_mag,sigma_mag,x1(:,i-km:i),x1_hat(:,i-km+1:i+1),u1(:,i-km:i),x_ind,u_ind,dt);
    end
    K.magnet = magnet(i);
        
    if static(i) == 0 || magnet(i) == 0 % if there is an opportunity to update estimates
        % step 6: update error state kinematics
        if i <= km
            [dx1(:,i),P1(:,:,i+1)] = errStateUpd(x1_hat(:,i+1),u1(:,i),R,P1_hat(:,:,i+1),g,x_ind,u_ind,K);
        else
            [dx1(:,i),P1(:,:,i+1)] = errStateUpd(x1_hat(:,i:i+1),u1(:,i-km:i),R,P1_hat(:,:,i+1),g,x_ind,u_ind,K);
        end
        
        % step 7: update nominal state kinematics
        x1(:,i+1) = nomStateUpd(x1_hat(:,i+1),dx1(:,i),x_ind,dx_ind);
    
    else % if there is not an opportunity to update estimates
        dx1(:,i+1) = zeros(size(dx1,1),1);
        x1(:,i+1) = x1_hat(:,i+1);
        P1(:,:,i+1) = P1_hat(:,:,i+1);
    end
    
end

%% plotting results

E11=quat2eul([x1(x_ind.q(1),2:i+1)',x1(x_ind.q(2:end),2:i+1)'],'ZYX');

figure;
hold on;
plot(time(1:i),x1(x_ind.q(1),2:i+1),'k');
plot(time(1:i),x1(x_ind.q(2),2:i+1),'r');
plot(time(1:i),x1(x_ind.q(3),2:i+1),'g');
plot(time(1:i),x1(x_ind.q(4),2:i+1),'b');
title('quaternions');
legend('q_{ER,0}','q_{ER,1}','q_{ER,2}','q_{ER,3}')
xlabel('time (s)');
ylabel('q');

figure; 
hold on;
plot(time(1:i),unwrap(E11(:,3)),'r');
plot(time(1:i),unwrap(E11(:,2)),'g');
plot(time(1:i),unwrap(E11(:,1)),'b');
title('euler angles');
legend('\phi_{ER}','\theta_{ER}','\psi_{ER}')
xlabel('time (s)');
ylabel('angle (deg)');