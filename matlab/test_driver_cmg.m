clear; clc;

% Time (s)
dt = 0.01;      % Step size
tf = 30;        % Final time
t  = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0    = 0;
wy0    = 0;
wz0    = 0;
psi0   = 0;
theta0 = 0;
phi0   = 0;

% CMG gimbal angles (deg). You can change these angle vectors to
% validate your model and try out the different cases required in
% the project. When you send me your function, I will try out some
% angle vectors to see if your model accurately predicts the response.
% th1 = zeros(size(t));
% th2 = zeros(size(t));
% th3 = zeros(size(t));
% th1 = 15*(1./(1+exp(-0.3*t))-0.5);
% th2 = 5*sin(2*pi/30*t);
% th3 = -5*sin(2*pi/30*t);
% th1 = 15*sin(2*pi/30*t);
% th2 = zeros(size(t));
% th3 = zeros(size(t));
th1 =  10.0*cos(t);
th2 =  2.50*cos(t);
th3 =  4.00*sin(t);

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname1_lastname2_cmg() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=hurst_cmg(wx0,wy0,wz0,psi0,theta0,phi0,t,th1,th2,th3);

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


[wx,wy,wz,psi,theta,phi]=cmg_verbose(wx0,wy0,wz0,psi0,theta0,phi0,t);
% Plot results
figure (2)
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
