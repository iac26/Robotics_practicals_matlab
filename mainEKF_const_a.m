% Main code for the EKF implementation in the Robotics Practivals course
clear all;  clc;

%% Constants
g = 9.81;        % gravity [m/s^2] (could come from a gravity model)

% Values required to calculate initial value for k in Baro observation
% equation:
M = 0.02897;  % Air molar mass [kg/mol]
R = 8.3145;  % Universal gas constant [J/mol/K]
T = 298.15;  % Air temperature (carefull about the units)

%% Loading observations
load('Quad_test_data.mat')
gps_measurements  = observations.gps.z;
time_gps          = observations.gps.time;

% gps_measurements(round(end/4):round(end/2))  = [];
% time_gps(round(end/4):round(end/2))  = [];

baro_measurements = observations.baro.z;
time_baro         = observations.baro.time;


%% KALMAN filter uncertainty parameters
%Standard deviation of the GPS position white noise
sigma_z_gps  = 5;       % [m]

%Standard deviation of the BARO readings white noise
sigma_z_baro = 5;       % [Pa]

% Process model noise
sigma_a  = 1;           % [m/s^2 /s]
sigma_p0 = 0.01;           % [Pa /s]
sigma_k  = 1e-9;           % [s^2/m^2 /s]

%% Generate the corresponding coviariance matrices
%% INITIAL VALUES
%States
X_0 = [gps_measurements(1), 0, 0, baro_measurements(1), M/(R*T), gps_measurements(1)]';

%Covariances
P_0       = diag([25, 0.25, 0.25, 25, 1e-13, 25]);

% Initialize the A PRIORI estimate
x_tilde = X_0;
P_tilde = P_0;


% State vector initial A POSTERIORI covariance matrix P_0
x_hat = x_tilde;
P_hat = P_tilde;

% Observations covariance matrices
R_gps  = sigma_z_gps^2;
R_baro = sigma_z_baro^2;

% Calculation of Kalman filter matrices
% Process model noise [PSD]
Q = diag([sigma_a^2, sigma_p0^2, sigma_k^2]);

% Dynamic matrix
F = [ 0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0 ];

% Noise shaping matrix
G = [ 0 0 0;
      0 0 0;
      1 0 0;
      0 1 0;
      0 0 1;
      0 0 0];

% MISC
i = 1;
j = 1;
k = 1;
time_last = min(time_baro(i),time_gps(j));

%% EKF Solution computation
% Preallocation of variables
kk     = length(union(time_gps,time_baro));
X      = nan(6,kk);
P      = nan(6,6,kk);
P_diag = nan(6,kk); % Only diagonal elements of covariance matrix
t      = nan(1,kk);

while true
    try
        X(:,k) = x_hat;
        P(:,:,k) = P_hat;
        P_diag(:,k) = diag(P_hat);
        
        if time_baro(i)<time_gps(j)
            % Kalman filter prediction step:
            [time_last,x_tilde,P_tilde] = prediction(time_baro(i),time_last,x_hat,P_hat,Q,F,G);
            
            % Kalman filter update step for Baro measurements:
            [x_hat,P_hat] = BARO_update_const_a(baro_measurements(i),x_tilde,P_tilde,R_baro,g);
            t(k) = time_baro(i);
            i = i + 1;
        else
            % Kalman filter prediction step:
            [time_last,x_tilde,P_tilde] = prediction(time_gps(j),time_last,x_hat,P_hat,Q,F,G);
            
            % Kalman filter update step for GPS measurements:
            [x_hat,P_hat] = GPS_update_const_a(gps_measurements(j),x_tilde,P_tilde,R_gps);
            t(k) = time_gps(j);
            j = j + 1;
        end
        
        k = k+1;
    
    catch
        disp(k-1)

        break
    end
    
end

% Removing non-allocated values
X = X(:,1:k-1);
P = P(:,:,1:k-1);
P_diag = P_diag(:,1:k-1);
t = t(:,1:k-1);





%% Plotting


subplot(221);
hold on

top1 = X(1,:)+3*sqrt(P_diag(1,:));
bot1 = X(1,:)-3*sqrt(P_diag(1,:));

plot(X(1,:), 'b');
plot(top1, 'r-');
plot(bot1, 'r-');
title("h");


% subplot(232);
% hold on
% top2 = X(2,:)+3*sqrt(P_diag(2,:));
% bot2 = X(2,:)-3*sqrt(P_diag(2,:));
% 
% plot(X(2,:), 'b');
% plot(top2, 'r-');
% plot(bot2, 'r-');
% title("v");
% 
% subplot(233);
% hold on
% top3 = X(3,:)+3*sqrt(P_diag(3,:));
% bot3 = X(3,:)-3*sqrt(P_diag(3,:));
% 
% plot(X(3,:), 'b');
% plot(top3, 'r-');
% plot(bot3, 'r-');
% title("a");


subplot(222);
hold on
top4 = X(4,:)+3*sqrt(P_diag(4,:));
bot4 = X(4,:)-3*sqrt(P_diag(4,:));

plot(X(4,:), 'b');
plot(top4, 'r-');
plot(bot4, 'r-');
title("p_0");

subplot(223);
hold on
top5 = X(5,:)+3*sqrt(P_diag(5,:));
bot5 = X(5,:)-3*sqrt(P_diag(5,:));

plot(X(5,:), 'b');
plot(top5, 'r-');
plot(bot5, 'r-');
title("K");

subplot(224);
hold on
top6 = X(6,:)+3*sqrt(P_diag(6,:));
bot6 = X(6,:)-3*sqrt(P_diag(6,:));

plot(X(6,:), 'b');
plot(top6, 'r-');
plot(bot6, 'r-');
title("h_0");









