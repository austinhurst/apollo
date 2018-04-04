function [wx, wy, wz, psi, theta, phi] = cmg_verbose(wx0, wy0, wz0, psi0, theta0, phi0, t)
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

% Compute Spacecraft Inertia
simplified_case = false;
get_I(simplified_case);

%% Task E - Simulation and Calculations
t   = makeRow(t);
Ts  = t.';

x0    = [phi0, theta0, psi0, wx0, wy0, wz0]*pi/180.0;
[t,x] = ode45(@apollo_eom_cmg_verbose, t, x0);
x     = x*180.0/pi;
phi   = x(:,1);
theta = x(:,2);
psi   = x(:,3);
wx    = x(:,4);
wy    = x(:,5);
wz    = x(:,6);

% minimum and maximum values
phi_max   = max(phi);
phi_min   = min(phi);
theta_max = max(theta);
theta_min = min(theta);
psi_max   = max(psi);
psi_min   = min(psi);
wx_max    = max(wx);
wx_min    = min(wx);
wy_max    = max(wy);
wy_min    = min(wy);
wz_max    = max(wz);
wz_min    = min(wz);

end
function [I] = compute_I(simplified_case)
  %% Task A - Inertial Properties

  % Constants
  in2ft         = 1.0/12.0;   % inches to feet conversion (note that this is exact)
  g_english     = 32.174;     % gravity in ft/s^2, used to find mass of Spacecraft from launch weight
  slugftft2kgmm = 1.35581795; % slug*ft^2 to kg*m^2

  % A frame properties - English Units
  w_cm    =  9730.0;    % weight of the CM, pounds
  w_sm    =  9690.0;    % weight of the SM, pounds
  w_pr    = 37295.0;    % weight of the propellant, pounds
  cg_x_cm =  1043.1;    % x positions of center of gravity, CM, inches
  cg_y_cm =    -0.1;    % y positions of center of gravity, CM, inches
  cg_z_cm =     7.8;    % z positions of center of gravity, CM, inches
  cg_x_sm =   908.2;    % x positions of center of gravity, SM, inches
  cg_y_sm =     0.7;    % y positions of center of gravity, SM, inches
  cg_z_sm =    -0.6;    % z positions of center of gravity, SM, inches
  cg_x_pr =   905.9;    % x positions of center of gravity, propellant, inches
  cg_y_pr =     5.6;    % y positions of center of gravity, propellant, inches
  cg_z_pr =    -2.4;    % z positions of center of gravity, propellant, inches
  Ixx_cmA =  4474.0;    % moment of inertia for the CM about x, slug-ft^2
  Iyy_cmA =  3919.0;    % moment of inertia for the CM about y, slug-ft^2
  Izz_cmA =  3684.0;    % moment of inertia for the CM about z, slug-ft^2
  Ixx_smA =  6222.0;    % moment of inertia for the SM about x, slug-ft^2
  Iyy_smA = 10321.0;    % moment of inertia for the SM about y, slug-ft^2
  Izz_smA = 10136.0;    % moment of inertia for the SM about z, slug-ft^2
  Ixx_prA = 19162.0;    % moment of inertia for the propellant about x, slug-ft^2
  Iyy_prA = 19872.0;    % moment of inertia for the propellant about y, slug-ft^2
  Izz_prA = 26398.0;    % moment of inertia for the propellant about z, slug-ft^2

  if (simplified_case == true)
    cg_y_cm = 0.0;      % y positions of center of gravity, CM, inches
    cg_z_cm = 0.0;      % z positions of center of gravity, CM, inches
    cg_y_sm = 0.0;      % y positions of center of gravity, SM, inches
    cg_z_sm = 0.0;      % z positions of center of gravity, SM, inches
    cg_y_pr = 0.0;      % y positions of center of gravity, propellant, inches
    cg_z_pr = 0.0;      % z positions of center of gravity, propellant, inches
  end

  w    = w_cm + w_sm + w_pr;
  cg_xA = ((w_cm*cg_x_cm) + (w_sm*cg_x_sm) + (w_pr*cg_x_pr))/w*in2ft; % x position of the center of mass, A frame, ft
  cg_yA = ((w_cm*cg_y_cm) + (w_sm*cg_y_sm) + (w_pr*cg_y_pr))/w*in2ft; % y position of the center of mass, A frame, ft
  cg_zA = ((w_cm*cg_z_cm) + (w_sm*cg_z_sm) + (w_pr*cg_z_pr))/w*in2ft; % z position of the center of mass, A frame, ft
  cg_A  = [cg_xA; cg_yA; cg_zA];

  Icm = [Ixx_cmA, 0.0, 0.0; 0.0, Iyy_cmA, 0.0; 0.0, 0.0, Izz_cmA];
  Ism = [Ixx_smA, 0.0, 0.0; 0.0, Iyy_smA, 0.0; 0.0, 0.0, Izz_smA];
  Ipr = [Ixx_prA, 0.0, 0.0; 0.0, Iyy_prA, 0.0; 0.0, 0.0, Izz_prA];

  % Parallel Axis Theorem
  % I = Ig + m*(d.'*d*eye() - d*d.');
  m_cm = w_cm/g_english;                                        % mass of CM in slugs
  cg_cmA = [cg_x_cm; cg_y_cm; cg_z_cm]*in2ft;
  IcmA = Icm + m_cm*(cg_cmA.'*cg_cmA*eye(3) - cg_cmA*cg_cmA.');

  m_sm = w_sm/g_english;                                        % mass of CM in slugs
  cg_smA = [cg_x_sm; cg_y_sm; cg_z_sm]*in2ft;
  IsmA = Ism + m_sm*(cg_smA.'*cg_smA*eye(3) - cg_smA*cg_smA.');

  m_pr = w_pr/g_english;                                        % mass of CM in slugs
  cg_prA = [cg_x_pr; cg_y_pr; cg_z_pr]*in2ft;
  IprA = Ipr + m_pr*(cg_prA.'*cg_prA*eye(3) - cg_prA*cg_prA.');

  IA  = IcmA + IsmA + IprA;
  m = w/g_english;                                              % mass of Spacecraft in slugs
  I_english = IA - m*(cg_A.'*cg_A*eye(3) - cg_A*cg_A.');        % Inertia Matrix in the B frame (slug*ft^2)
  I = I_english*slugftft2kgmm;                                  % Inertia Matrix in the B frame (kg*m^2)
end
function [xdot] = apollo_eom_cmg_verbose(t,x)
  %% Task B - Equations of Motion
  persistent I
  persistent invI
  persistent p_dot
  persistent q_dot
  persistent r_dot
  if isempty(I) || t == 0.0
    I        = get_I();
    invI     = inv(I);
    p_dot    = 0.0;
    q_dot    = 0.0;
    r_dot    = 0.0;
  end

  phi       = x(1);
  theta     = x(2);
  psi       = x(3);
  p         = x(4);
  q         = x(5);
  r         = x(6);
  omega     = [p; q; r];

  M         = compute_M_verbose(t, p, q, r, p_dot, q_dot, r_dot);

  s_phi     = sin(phi);
  c_phi     = cos(phi);
  s_theta   = sin(theta);
  c_theta   = cos(theta);

  phi_dot   = 1.0/c_theta*(q*s_theta*s_phi + r*s_theta*c_phi)+p;
  theta_dot = q*c_phi - r*s_phi;
  psi_dot   = 1.0/c_theta*(q*s_phi + r*c_phi);
  omega_dot = invI*(M - cross(omega,(I*omega)));

  p_dot     = omega_dot(1);
  q_dot     = omega_dot(2);
  r_dot     = omega_dot(3);
  xdot      = [phi_dot; theta_dot; psi_dot; p_dot; q_dot; r_dot];
end
function [M] = compute_M_verbose(t, p, q, r, p_dot, q_dot, r_dot);
  %% Task D - Control Moment Gyroscope Model
  persistent IR
  persistent Omega
  if (t == 0 || isempty(IR))
    Omega        = 6600.0*2.0*pi/60.0;
    IR           = get_IR();
  end

  theta1         =  10.0*cos(t)*pi/180.0;
  theta2         =  2.50*cos(t)*pi/180.0;
  theta3         =  4.00*sin(t)*pi/180.0;
  theta1_dot     = -10.0*sin(t)*pi/180.0;
  theta2_dot     = -2.50*sin(t)*pi/180.0;
  theta3_dot     =  4.00*cos(t)*pi/180.0;
  theta1_dot_dot = -10.0*cos(t)*pi/180.0;
  theta2_dot_dot = -2.50*cos(t)*pi/180.0;
  theta3_dot_dot = -4.00*sin(t)*pi/180.0;

  % CMG-1
  I       = zeros(3);
  I(1,1)  = IR;
  R1      = [cos(theta1), sin(theta1), 0.0; -sin(theta1), cos(theta1), 0.0; 0.0, 0.0, 1.0];
  w       = [Omega; 0.0; theta1_dot] + R1*[p; q; r];
  w_dot   = [
             -theta1_dot*sin(theta1)*p + p_dot*cos(theta1) + theta1_dot*cos(theta1)*q + q_dot*sin(theta1);
             -theta1_dot*cos(theta1)*p - p_dot*sin(theta1) - theta1_dot*sin(theta1)*q + q_dot*cos(theta1);
              theta1_dot_dot + r_dot;
            ];
  w_frame = [0.0; 0.0; theta1_dot] + R1*[p; q; r];
  M1      = I*w_dot + cross(w_frame,I*w);
  Mcsm1   = -R1.'*M1;

  % CMG-2
  I       = zeros(3);
  I(2,2)  = IR;
  R2      = [1.0, 0.0, 0.0; 0.0, cos(theta2), sin(theta2); 0.0, -sin(theta2), cos(theta2)];
  w       = [theta2_dot; Omega; 0.0] + R2*[p; q; r];
  w_dot   = [
              theta2_dot_dot + p_dot;
             -theta2_dot*sin(theta2)*q + q_dot*cos(theta2) + theta2_dot*cos(theta2)*r + r_dot*sin(theta2);
             -theta2_dot*cos(theta2)*q - q_dot*sin(theta2) - theta2_dot*sin(theta2)*r + r_dot*cos(theta2);
            ];
  w_frame = [theta2_dot; 0.0; 0.0] + R2*[p; q; r];
  M2      = I*w_dot + cross(w_frame,I*w);
  Mcsm2   = -R2.'*M2;

  % CMG-3
  I       = zeros(3);
  I(3,3)  = IR;
  R3      = [cos(theta3), 0.0, -sin(theta3); 0.0, 1.0, 0.0; sin(theta3), 0.0, cos(theta3)];
  w       = [0.0; theta3_dot; Omega] + R3*[p; q; r];
  w_dot   = [
             -theta3_dot*sin(theta3)*p + p_dot*cos(theta3) - theta3_dot*cos(theta3)*r - r_dot*sin(theta3);
              theta3_dot_dot + q_dot;
              theta3_dot*cos(theta3)*p + p_dot*sin(theta3) - theta3_dot*sin(theta3)*r + r_dot*cos(theta3);
            ];
  w_frame = [0.0; theta3_dot; 0.0] + R3*[p; q; r];
  M3      = I*w_dot + cross(w_frame,I*w);
  Mcsm3   = -R3.'*M3;

  M = Mcsm1 + Mcsm2 + Mcsm3;
end

%% Storage Functions
function [I_out] = get_I(varargin)
  persistent I
  if nargin == 1
    I = compute_I(varargin{1});
  end
  if isempty(I)
    error('set I first');
  end
  I_out = I;
end
function [IR] = get_IR()
  persistent IR_s
  if (isempty(IR_s))
    max_th_dot = 30*pi/180.0;                % maximum thetaN gimbal rate
    I          = get_I();
    Omega      = 6600.0*2.0*pi/60.0;

    % x axis req
    req_x_dd   = 5.0*pi/180.0;               % required x axis acceleration rad/s^2
    M_req      = I*[req_x_dd; 0.0; 0.0];
    IR_min_x   = M_req(1)/(max_th_dot*Omega);
    % y axis req
    req_y_dd   = 2.5*pi/180.0;               % required y axis acceleration rad/s^2
    M_req      = I*[0.0; req_y_dd; 0.0];
    IR_min_y   = M_req(2)/(max_th_dot*Omega);
    % z axis req
    req_z_dd   = 2.5*pi/180.0;               % required z axis acceleration rad/s^2
    M_req      = I*[0.0; 0.0; req_z_dd];
    IR_min_z   = M_req(3)/(max_th_dot*Omega);

    IR_s       = max([IR_min_x, IR_min_y, IR_min_z]);
  end
  IR = IR_s;
end

%% Miscellaneous Functions
function [col] = makeRow(input)
  if iscolumn(input)
    input = input.';
  end
  col = input;
end
