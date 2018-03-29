clear; clc; close all;

% Time (s)
dt = 0.01;      % Step size
tf = 1000;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 1;
wy0 = 0;
wz0 = 0;
psi0 = 0;
theta0 = 0;
phi0 = 0;

% Torques (N-m). You can change these torque vectors to validate your
% model and try out the different cases required in the project. When
% you send me your function, I will try out some torques to see if
% your model accurately predicts the response.
% omega = wx0*pi/180;
% Mx = zeros(size(t));
% My = 3179*omega^2*ones(size(t));
% Mz = 1538*omega^2*ones(size(t));


% Part C.b
% Mx = 176*cos(0.2*t);
% My = 54*ones(size(t));
% Mz = 98*sin(0.3*t);
% Part C.c.1
Mx = zeros(size(t));
My = -0.968471643039707*ones(size(t));
Mz = -0.468444068275835*ones(size(t));
% Part C.c.2
% Mx = zeros(size(t));
% My = zeros(size(t));
% Mz = zeros(size(t));

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname1_lastname2() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=hurst(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

% Plot results
figure (1)
subplot(2,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z');
subplot(2,1,2);
plot(t,psi,t,theta,t,phi);
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi');
