function [wx, wy, wz, psi, theta, phi] = hurst(wx0, wy0, wz0, psi0, theta0, phi0, t, Mx, My, Mz)
%% Apollo Spacecraft Attitude Dynamics
% Austin Hurst

%% General Notes
% CM      - Command Module
% SM      - Service Module
% CSM     - Command and Service Module
% 3-2-1 Euler Angle sequernces (yaw (psi), pitch (theta), roll(phi))
% B Frame - Body fixed frame at the center of mass
%           x axis points out of the tip of the CM, (roll)
%           y axis points out of the side, (pitch)
%           z axis points out of the belly, (yaw)
% A Frame - Body fixed frame 1000 inches in -x from base of the CM
%           All axes parallel to the B frame

%% Constants
in2m          = 0.0254;     % inches to meters conversion (note that this is exact)
g_english     = 32.174;     % gravity in ft/s^2, used to find mass of Spacecraft from launch weight
slugs2kg      = 14.5939;    % slugs to kilogram conversion, used to find mass of Spacecraft
slugftft2kgmm = 1.35581795; % slug*ft^2 to kg*m^2
slugs2lbm     = 32.174;     % conversion from slugs to pound mass

%% Task A - Inertial Properties
% a. General Case
%    i.  The location of the center of mass of the CSM, expressed in meters, as a vector relative to the A frame.
%        Include the propellant in your calculation.
%    ii. The total inertia matrix of the CSM, including propellant, about the body-fixed B axes at the center of mass
%        of the CSM. Express your answer as a matrix with units of kg-m^2
% b. Simplified Case, assuming the centers of mass of all components are located on the x-axis
%     i.  The location of the center of mass of the CSM, expressed in meters, as a vector relative to the A frame.
%         Include the propellant in your calculation.
%     ii. The total inertia matrix of the CSM, including propellant, about the body-fixed B axes at the center of mass
%         of the CSM. Express your answer as a matrix with units of kg-m^2
simplified_case = false;
% A frame properties - English Units
w_cm    =  9730.0; % weight of the CM, pounds
w_sm    =  9690.0; % weight of the SM, pounds
w_pr    = 37295.0; % weight of the propellant, pounds
cg_x_cm =  1043.1; % x positions of center of gravity, CM, inches
cg_y_cm =    -0.1; % y positions of center of gravity, CM, inches
cg_z_cm =     7.8; % z positions of center of gravity, CM, inches
cg_x_sm =   908.2; % x positions of center of gravity, SM, inches
cg_y_sm =     0.7; % y positions of center of gravity, SM, inches
cg_z_sm =    -0.6; % z positions of center of gravity, SM, inches
cg_x_pr =   905.9; % x positions of center of gravity, propellant, inches
cg_y_pr =     5.6; % y positions of center of gravity, propellant, inches
cg_z_pr =    -2.4; % z positions of center of gravity, propellant, inches
Ixx_cmA =  4474.0; % moment of inertia for the CM about x, slug-ft^2
Iyy_cmA =  3919.0; % moment of inertia for the CM about y, slug-ft^2
Izz_cmA =  3684.0; % moment of inertia for the CM about z, slug-ft^2
Ixx_smA =  6222.0; % moment of inertia for the SM about x, slug-ft^2
Iyy_smA = 10321.0; % moment of inertia for the SM about y, slug-ft^2
Izz_smA = 10136.0; % moment of inertia for the SM about z, slug-ft^2
Ixx_prA = 19162.0; % moment of inertia for the propellant about x, slug-ft^2
Iyy_prA = 19872.0; % moment of inertia for the propellant about y, slug-ft^2
Izz_prA = 26398.0; % moment of inertia for the propellant about z, slug-ft^2

if (simplified_case == true)
  cg_y_cm = 0.0; % y positions of center of gravity, CM, inches
  cg_z_cm = 0.0; % z positions of center of gravity, CM, inches
  cg_y_sm = 0.0; % y positions of center of gravity, SM, inches
  cg_z_sm = 0.0; % z positions of center of gravity, SM, inches
  cg_y_pr = 0.0; % y positions of center of gravity, propellant, inches
  cg_z_pr = 0.0; % z positions of center of gravity, propellant, inches
end

w    = w_cm + w_sm + w_pr;
cg_xA = ((w_cm*cg_x_cm) + (w_sm*cg_x_sm) + (w_pr*cg_x_pr))/w*in2m; % x position of the center of mass, A frame, meters
cg_yA = ((w_cm*cg_y_cm) + (w_sm*cg_y_sm) + (w_pr*cg_y_pr))/w*in2m; % y position of the center of mass, A frame, meters
cg_zA = ((w_cm*cg_z_cm) + (w_sm*cg_z_sm) + (w_pr*cg_z_pr))/w*in2m; % z position of the center of mass, A frame, meters
cg_A  = [cg_xA; cg_yA; cg_zA];

Icm = [Ixx_cmA, 0.0, 0.0; 0.0, Iyy_cmA, 0.0; 0.0, 0.0, Izz_cmA];
Ism = [Ixx_smA, 0.0, 0.0; 0.0, Iyy_smA, 0.0; 0.0, 0.0, Izz_smA];
Ipr = [Ixx_prA, 0.0, 0.0; 0.0, Iyy_prA, 0.0; 0.0, 0.0, Izz_prA];

% Parallel Axis Theorem
% I = Ig + m*(d.'*d*eye() - d*d.');
m_cm = w_cm/g_english*slugs2lbm;                              % mass of CM in pounds mass
cg_cmA = [cg_x_cm; cg_y_cm; cg_z_cm];
IcmA = Icm + m_cm*(cg_cmA.'*cg_cmA*eye(3) - cg_cmA*cg_cmA.');

m_sm = w_sm/g_english*slugs2lbm;                              % mass of CM in pounds mass
cg_smA = [cg_x_sm; cg_y_sm; cg_z_sm];
IsmA = Ism + m_sm*(cg_smA.'*cg_smA*eye(3) - cg_smA*cg_smA.');

m_pr = w_pr/g_english*slugs2lbm;                              % mass of CM in pounds mass
cg_prA = [cg_x_pr; cg_y_pr; cg_z_pr];
IprA = Ipr + m_pr*(cg_prA.'*cg_prA*eye(3) - cg_prA*cg_prA.');

IA  = IcmA + IsmA + IprA;
m = w/g_english*slugs2kg;                                     % mass of Spacecraft in kg
I_english = IA - m*(cg_A.'*cg_A*eye(3) - cg_A*cg_A.');        % Inertia Matrix in the B frame (english)
I = I_english*slugftft2kgmm;                                  % Inertia Matrix in the B frame (metric)

%% Task B - Equations of Motion
%% Task C -
%% Task D -
%% Task E -
%% Task F -
end
