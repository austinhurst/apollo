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

simplified_case = false;
get_I(simplified_case);

%% Task C - Simulation and Calculations
Mx = makeRow(Mx);
My = makeRow(My);
Mz = makeRow(Mz);
t  = makeRow(t);

Ms = [Mx.', My.', Mz.'];
Ts = t.';
getTsMs(Ts, Ms);

x0    = [phi0, theta0, psi0, wx0, wy0, wz0]*pi/180.0;
[t,x] = ode45(@apollo_eom, t, x0);
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

% Passive Thermal Control
% find moments needed to maintain a desired omega
omega_d = [1.0; 0.0; 0.0]*pi/180;
I = get_I();
M = cross((I*omega_d),omega_d);

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
function [xdot] = apollo_eom(t,x)
  %% Task B - Equations of Motion
  persistent I
  persistent invI
  persistent Ts
  persistent Ms
  if isempty(I) || t == 0.0
    I        = get_I();
    invI     = inv(I);
    [Ts, Ms] = getTsMs();
  end

  M = interp1(Ts, Ms, t).';
  phi       = x(1);
  theta     = x(2);
  psi       = x(3);
  p         = x(4);
  q         = x(5);
  r         = x(6);
  omega     = [p; q; r];

  s_phi     = sin(phi);
  c_phi     = cos(phi);
  s_theta   = sin(theta);
  c_theta   = cos(theta);

  phi_dot   = 1.0/c_theta*(q*s_theta*s_phi + r*s_theta*c_phi)+p;
  theta_dot = q*c_phi - r*s_phi;
  psi_dot   = 1.0/c_theta*(q*s_phi + r*c_phi);
  omega_dot = invI*(M - cross((I*omega),omega));

  p_dot     = omega_dot(1);
  q_dot     = omega_dot(2);
  r_dot     = omega_dot(3);
  xdot      = [phi_dot; theta_dot; psi_dot; p_dot; q_dot; r_dot];
end

%% Storage Functions
function [Ts, Ms] = getTsMs(varargin)
  persistent M_in
  persistent T_in
  if nargin == 2
    T_in = varargin{1};
    M_in = varargin{2};
  end
  Ms = M_in;
  Ts = T_in;
end
function [I_out] = get_I(varargin)
  persistent I
  if nargin == 1
    I = compute_I(varargin{1});
  end
  I_out = I;
end

%% Miscellaneous Functions
function [col] = makeRow(input)
  if iscolumn(input)
    input = input.';
  end
  col = input;
end
