function [wx, wy, wz, psi, theta, phi] = hurst_cmg(wx0, wy0, wz0, psi0, theta0, phi0, t, th1, th2, th3)
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
th1 = makeRow(th1*pi/180.0);
th2 = makeRow(th2*pi/180.0);
th3 = makeRow(th3*pi/180.0);
t   = makeRow(t);

ths = [th1.', th2.', th3.'];
Ts  = t.';
getTsths(Ts, ths);

x0    = [phi0, theta0, psi0, wx0, wy0, wz0]*pi/180.0;
[t,x] = ode45(@apollo_eom_cmg, t, x0);
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
function [xdot] = apollo_eom_cmg(t,x)
  %% Task B - Equations of Motion
  persistent I
  persistent invI
  persistent Ts
  persistent ths
  if isempty(I) || t == 0.0
    I        = get_I();
    invI     = inv(I);
    [Ts, ths] = getTsths();
  end

  th        = interp1(Ts, ths, t).';
  phi       = x(1);
  theta     = x(2);
  psi       = x(3);
  p         = x(4);
  q         = x(5);
  r         = x(6);
  omega     = [p; q; r];

  M         = compute_M(t, th);

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
function [M] = compute_M(t, th)
  %% Task D - Control Moment Gyroscope Model
  persistent IR
  persistent Omega
  persistent t_last
  persistent tau
  persistent th_der_last
  persistent th_last
  if (t == 0 || isempty(IR))
    Omega       = 6600.0*2.0*pi/60.0;
    IR          = get_IR();
    t_last      = 0.0;
    tau         = 0.15;               % tunning parameter to get the derivative right
    th_der_last = [0.0; 0.0; 0.0];
    th_last     = [0.0; 0.0; 0.0];
  end

  ts         = t - t_last;
  theta1     = th(1);
  theta2     = th(2);
  theta3     = th(3);
  theta1_dot = time_derivative(theta1, th_last(1), th_der_last(1), ts, tau);
  theta2_dot = time_derivative(theta2, th_last(2), th_der_last(2), ts, tau);
  theta3_dot = time_derivative(theta3, th_last(3), th_der_last(3), ts, tau);

  % CMG-1
  Mcsm1 = [
            sin(theta1)*IR*Omega*theta1_dot;
           -cos(theta1)*IR*Omega*theta1_dot;
            0.0;
          ];
  % CMG-2
  Mcsm2 = [
            0.0;
           -sin(theta2)*IR*Omega*theta2_dot;
           -cos(theta2)*IR*Omega*theta2_dot;
          ];
  % CMG-3
  Mcsm3 = [
           -cos(theta3)*IR*Omega*theta3_dot;
            0.0;
           -sin(theta3)*IR*Omega*theta3_dot;
          ];
  M = Mcsm1 + Mcsm2 + Mcsm3;

  % age data
  t_last      = t;
  th_last     = [theta1; theta2; theta3];
  th_der_last = [ theta1_dot; theta2_dot; theta3_dot];
end

%% Storage Functions
function [Ts, ths] = getTsths(varargin)
  persistent th_in
  persistent T_in
  if nargin == 2
    T_in = varargin{1};
    th_in = varargin{2};
  end
  ths = th_in;
  Ts = T_in;
end
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
function der = time_derivative(val, val_last, der_last, ts, tau)
  der = (2.0*tau - ts)/(2.0*tau + ts)*der_last + 2.0/(2.0*tau + ts)*(val - val_last);
end
