#include <apollo/pid.h>

namespace apollo
{
PID::PID()
{
  // SETUP THE CONTROLLER HERE
  float Ixx, Iyy, Izz, tr, zeta;
  getRosParam("vehicle_description/Ixx",  Ixx);
  getRosParam("vehicle_description/Iyy",  Iyy);
  getRosParam("vehicle_description/Izz",  Izz);
  getRosParam("rise_time",                tr);
  getRosParam("zeta",                     zeta);

  float wn = 2.2f/tr;
  // PHI
  kP_phi_ = wn*wn*Ixx;
  kD_phi_ = 2.0f*zeta*wn*Ixx;

  // THETA
  kP_theta_ = wn*wn*Iyy;
  kD_theta_ = 2.0f*zeta*wn*Iyy;

  // PSI
  kD_psi_ = 0.0f;
  kP_psi_ = 5.0f*(Izz + kD_psi_)/tr;

  sigma_         = 0.05;
  r_last_        = 0.0;
  rd_            = 0.0;
}
void PID::control(const ros::TimerEvent& event)
{
  ros::Time new_time = ros::Time::now();
  mapControlChannels();
  ros::Duration time_step = new_time - last_time_;
  ts_ = time_step.toSec();

  // PHI - PD CONTROL
  float e_phi     = roll_desired_ - state_.phi;
  float tau_phi   = kP_phi_*e_phi  - kD_phi_*state_.p;

  // THETA - PD CONTROL
  float e_theta   = pitch_desired_ - state_.theta;
  float tau_theta = kP_theta_*e_theta  - kD_theta_*state_.q;

  // PSI - PD Control
  rd_             = ((2.0f*sigma_ - ts_)/(2.0f*sigma_ + ts_))*rd_ + (2.0f/(2.0f*sigma_ + ts_))*(state_.r - r_last_);
  float e_psi     = yaw_rate_desired_ - state_.r;
  float tau_psi   = kP_psi_*e_psi  - kD_psi_*rd_;

  // OUTPUT
  Mx_ = tau_phi;
  My_ = tau_theta;
  Mz_ = tau_psi;

  publishMomentsCommand();
  publishDesiredCommand();

  // Age Data
  last_time_ =  new_time;
  r_last_    =  state_.r;
}
float PID::saturate(float value_i, float min, float max)
{
  float  value_o  = value_i < min ? min : value_i;
  return value_o  = value_o > max ? max : value_o;
}
} // end namespace apollo
