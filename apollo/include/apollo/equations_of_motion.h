#ifndef EQUATIONS_OF_MOTION_H
#define EQUATIONS_OF_MOTION_H

#include <cmath>
#include <string>
#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <tf/LinearMath/Quaternion.h>
#include <apollo/VehicleState.h>
#include <apollo/state_t.h>

#include <geometry_msgs/Vector3.h>

namespace apollo
{
class EquationsOfMotion
{
public:
  EquationsOfMotion();
  ~EquationsOfMotion();
private:
  //********************* NODE HANDLES *********************//
  ros::NodeHandle nh_;         // public node handle for publishing, subscribing

  //************** SUBSCRIBERS AND PUBLISHERS **************//
  ros::Subscriber moment_subscriber_;
  ros::Publisher truth_publisher_;
  tf::TransformBroadcaster pose_broadcaster_;

  //******************** CLASS VARIABLES *******************//
protected:
  float Ixx_;                 // moment  of inertia about i^b in Kg*m^2
  float Iyy_;                 // moment  of inertia about j^b in Kg*m^2
  float Izz_;                 // moment  of inertia about k^b in Kg*m^2
  float Ixy_;                 // product of inertia about i^b in Kg*m^2
  float Ixz_;                 // product of inertia about j^b in Kg*m^2
  float Iyz_;                 // product of inertia about k^b in Kg*m^2
  float invI11_;
  float invI12_;
  float invI13_;
  float invI21_;
  float invI22_;
  float invI23_;
  float invI31_;
  float invI32_;
  float invI33_;
  float piD30_;

  // Truth Variables
  state_t state_;

  float Mx_;
  float My_;
  float Mz_;

private:
  // RK4 Variables
  state_t k1_;
  state_t k2_;
  state_t k3_;
  state_t k4_;
  ros::Time last_time_;

private:
  VehicleState truth_msg_;
  geometry_msgs::TransformStamped odom_trans_;

  //***************** CALLBACKS AND TIMERS *****************//
  void propogate(const ros::TimerEvent& event);
  ros::Timer propogate_timer_;
  void updateViz(const ros::WallTimerEvent& event);
  ros::WallTimer update_viz_timer_;
  void momentCallback(const geometry_msgs::Vector3ConstPtr &ms);

  //********************** FUNCTIONS ***********************//
protected:
  state_t derivative(state_t state);

};// end class EquationsOfMotion
} // end namespace apollo

#endif // EQUATIONS_OF_MOTION_H
