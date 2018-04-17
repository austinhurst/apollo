#include <apollo/equations_of_motion.h>

namespace apollo
{
EquationsOfMotion::EquationsOfMotion() :
  nh_(ros::NodeHandle())
{
  //********************** PARAMETERS **********************//
  // Simulation Parameters
  float propogate_rate, update_viz_rate;
  if (!(ros::param::get("sim/propogate_rate",propogate_rate)))
    ROS_WARN("No param named 'propogate_rate'");
  if (!(ros::param::get("/apollo/ground_station/update_viz_rate",update_viz_rate)))
    ROS_WARN("No param named 'update_viz_rate'");
  if (!(ros::param::get("/apollo/vehicle_description/Ixx",Ixx_)))
    ROS_WARN("No param named 'Ixx");
  if (!(ros::param::get("/apollo/vehicle_description/Iyy",Iyy_)))
    ROS_WARN("No param named 'Iyy");
  if (!(ros::param::get("/apollo/vehicle_description/Izz",Izz_)))
    ROS_WARN("No param named 'Izz");
  if (!(ros::param::get("/apollo/vehicle_description/Ixy",Ixy_)))
    ROS_WARN("No param named 'Ixy");
  if (!(ros::param::get("/apollo/vehicle_description/Ixz",Ixz_)))
    ROS_WARN("No param named 'num_motors");
  if (!(ros::param::get("/apollo/vehicle_description/Iyz",Iyz_)))
    ROS_WARN("No param named 'Iyz");

  if (!(ros::param::get("/apollo/initial/phi",state_.phi)))
    ROS_WARN("No param named 'phi'");
  if (!(ros::param::get("/apollo/initial/theta",state_.theta)))
    ROS_WARN("No param named 'theta'");
  if (!(ros::param::get("/apollo/initial/psi",state_.psi)))
    ROS_WARN("No param named 'psi'");
  if (!(ros::param::get("/apollo/initial/p",state_.p)))
    ROS_WARN("No param named 'p'");
  if (!(ros::param::get("/apollo/initial/q",state_.q)))
    ROS_WARN("No param named 'q'");
  if (!(ros::param::get("/apollo/initial/r",state_.r)))
    ROS_WARN("No param named 'r'");

  //************** SUBSCRIBERS AND PUBLISHERS **************//
  moment_subscriber_ = nh_.subscribe("moment_command",1,&EquationsOfMotion::momentCallback, this);
  truth_publisher_ = nh_.advertise<apollo::VehicleState>("truth",1);

  //******************** CLASS VARIABLES *******************//

  odom_trans_.header.frame_id = "odom";
  odom_trans_.child_frame_id = "base_link";
  last_time_ = ros::Time::now();

  //***************** CALLBACKS AND TIMERS *****************//
  propogate_timer_ = nh_.createTimer(ros::Duration(1.0/propogate_rate), &EquationsOfMotion::propogate, this);
  update_viz_timer_ = nh_.createWallTimer(ros::WallDuration(1.0/update_viz_rate), &EquationsOfMotion::updateViz, this);

  //********************** FUNCTIONS ***********************//
  float det_I = Ixx_*(Iyy_*Izz_ - Iyz_*Iyz_) - Ixy_*(Ixy_*Izz_ + Ixz_*Iyz_) - Ixz_*(Ixy_*Iyz_ + Ixz_*Iyy_);
  invI11_     = (Iyy_*Izz_ - Iyz_*Iyz_)/det_I;
  invI12_     = (Ixy_*Izz_ + Ixz_*Iyz_)/det_I;
  invI13_     = (Ixy_*Iyz_ + Ixz_*Iyy_)/det_I;
  invI21_     = invI12_;
  invI22_     = (Ixx_*Izz_ - Ixz_*Ixz_)/det_I;
  invI23_     = (Ixx_*Iyz_ + Ixy_*Ixz_)/det_I;
  invI31_     = invI13_;
  invI32_     = invI23_;
  invI33_     = (Ixx_*Iyy_ - Ixy_*Ixy_)/det_I;
  Mx_ = 0.0f;
  My_ = 0.0f;
  Mz_ = 0.0f;
}
EquationsOfMotion::~EquationsOfMotion()
{

}
void EquationsOfMotion::propogate(const ros::TimerEvent&)
{
  // Runge-Kutta 4th Order - Compute Truth
  ros::Time new_time = ros::Time::now();
  float h = (new_time - last_time_).toSec();
  k1_ = derivative(state_               );
  k2_ = derivative(state_ + k1_*(h/2.0f));
  k3_ = derivative(state_ + k2_*(h/2.0f));
  k4_ = derivative(state_ + k3_*h       );
  state_ = state_ + (k1_ + k2_*2.0f + k3_*2.0f + k4_)*(h/6.0f);
  truth_publisher_.publish(state_.msg());
  last_time_ = new_time;
}
void EquationsOfMotion::momentCallback(const geometry_msgs::Vector3ConstPtr &msg)
{
  Mx_ = msg->x;
  My_ = msg->y;
  Mz_ = msg->z;
}
state_t EquationsOfMotion::derivative(state_t s)
{
  float c_phi   = cosf(s.phi);
  float s_phi   = sinf(s.phi);

  //********** EQUATIONS OF MOTION ************//
  state_t derivative_of;

  // KINEMATICS
  // solve for (phi_dot, theta_dot, psi_dot)^T -->  (p, q, r)^T = (phi_dot, 0, 0)^T + Rv2_b*(0, theta_dot, 0)^T + Rv1_b*(0, 0, psi_dot)^T
  float t_theta       = tanf(s.theta);
  float sec_theta     = 1.0f/cosf(s.theta);
  derivative_of.phi   = s.p +   s_phi*t_theta*s.q +   c_phi*t_theta*s.r;
  derivative_of.theta =                 c_phi*s.q -           s_phi*s.r;
  derivative_of.psi   =       s_phi*sec_theta*s.q + c_phi*sec_theta*s.r;

  // DYNAMICS
  // (p_dot, q_dot, r_dot)^T = I^-1*(-omega_b/i^b x (I*omega_b/i^b) + m^b)          omega_b/i^b = (p,q,r)^T
  float negwJwm_i = -s.q*(-Ixz_*s.p - Iyz_*s.q + Izz_*s.r) + s.r*(-Ixy_*s.p + Iyy_*s.q - Iyz_*s.r) + Mx_;
  float negwJwm_j =  s.p*(-Ixz_*s.p - Iyz_*s.q + Izz_*s.r) - s.r*( Ixx_*s.p - Ixy_*s.q - Ixz_*s.r) + My_;
  float negwJwm_k = -s.p*(-Ixy_*s.p + Iyy_*s.q - Iyz_*s.r) + s.q*( Ixx_*s.p - Ixy_*s.q - Ixz_*s.r) + Mz_;
  derivative_of.p = invI11_*negwJwm_i + invI12_*negwJwm_j + invI13_*negwJwm_k;
  derivative_of.q = invI21_*negwJwm_i + invI22_*negwJwm_j + invI23_*negwJwm_k;
  derivative_of.r = invI31_*negwJwm_i + invI32_*negwJwm_j + invI33_*negwJwm_k;

  return derivative_of;
}
void EquationsOfMotion::updateViz(const ros::WallTimerEvent&)
{
  // Pull in State (truth), translate to quaternion, broadcast tf
  odom_trans_.header.stamp = ros::Time::now();
  odom_trans_.transform.translation.x =  0.0;
  odom_trans_.transform.translation.y =  0.0;
  odom_trans_.transform.translation.z =  0.0;
  // Euler Angles in NED to Quaternion in NED to Quaternion in XYZ
  float qx,qy,qz,qw;
  qx =   cosf(state_.psi/2.0f)*sinf(state_.theta/2.0f)*cosf(state_.phi/2.0f)\
       + sinf(state_.psi/2.0f)*cosf(state_.theta/2.0f)*sinf(state_.phi/2.0f);
  qy =   cosf(state_.psi/2.0f)*cosf(state_.theta/2.0f)*sinf(state_.phi/2.0f)\
       - sinf(state_.psi/2.0f)*sinf(state_.theta/2.0f)*cosf(state_.phi/2.0f);
  qz = - sinf(state_.psi/2.0f)*cosf(state_.theta/2.0f)*cosf(state_.phi/2.0f)\
       + cosf(state_.psi/2.0f)*sinf(state_.theta/2.0f)*sinf(state_.phi/2.0f);
  qw =   cosf(state_.psi/2.0f)*cosf(state_.theta/2.0f)*cosf(state_.phi/2.0f)\
       + sinf(state_.psi/2.0f)*sinf(state_.theta/2.0f)*sinf(state_.phi/2.0f);
  tf::Quaternion q(qx,qy,qz,qw);
  tf::quaternionTFToMsg(q,odom_trans_.transform.rotation);
  pose_broadcaster_.sendTransform(odom_trans_);
}
} // end namespace apollo

//********************************************************//
//************************ MAIN **************************//
//********************************************************//
int main(int argc, char** argv)
{
  ros::init(argc, argv, "equations_of_motion");
  ros::NodeHandle nh("apollo");

  apollo::EquationsOfMotion eom_obj;
  ros::spin();

  return 0;
} // end main
