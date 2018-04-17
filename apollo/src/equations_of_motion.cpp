#include <apollo/equations_of_motion_base.h>

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
  if (!(ros::param::get("/apollo/vehicle_description/mass",mass_)))
    ROS_WARN("No param named 'mass");
  if (!(ros::param::get("/apollo/vehicle_description/Ixx",Ixx_)))
    ROS_WARN("No param named 'Jx");
  if (!(ros::param::get("/apollo/vehicle_description/Iyy",Iyy_)))
    ROS_WARN("No param named 'Jy");
  if (!(ros::param::get("/apollo/vehicle_description/Izz",Izz_)))
    ROS_WARN("No param named 'Jz");
  if (!(ros::param::get("/apollo/vehicle_description/Ixy",Ixy_)))
    ROS_WARN("No param named 'Jxy");
  if (!(ros::param::get("/apollo/vehicle_description/Ixz",Ixz_)))
    ROS_WARN("No param named 'num_motors");
  if (!(ros::param::get("/apollo/vehicle_description/Iyz",Iyz_)))
    ROS_WARN("No param named 'Jyz");

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
  truth_publisher_ = nh_.advertise<VehicleState>("truth",1);

  //******************** CLASS VARIABLES *******************//

  odom_trans_.header.frame_id = "odom";
  odom_trans_.child_frame_id = "base_link";
  srand(time(0));
  last_time_ = ros::Time::now();

  //***************** CALLBACKS AND TIMERS *****************//
  propogate_timer_ = nh_.createTimer(ros::Duration(1.0/propogate_rate), &EquationsOfMotion::propogate, this);
  update_viz_timer_ = nh_.createWallTimer(ros::WallDuration(1.0/update_viz_rate), &EquationsOfMotion::updateViz, this);

  //********************** FUNCTIONS ***********************//

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
apollo::state_struct EquationsOfMotion::derivative(apollo::state_struct)
{
  ROS_ERROR("CHILD CLASS FUNCTION 'derivative' WAS NOT CALLED.");
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

  apollo_sim::EquationsOfMotion eom_obj;
  ros::spin();
  
  return 0;
} // end main
