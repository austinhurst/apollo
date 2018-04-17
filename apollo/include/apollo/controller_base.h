#ifndef CONTROLLER_BASE_H
#define CONTROLLER_BASE_H

#include <cmath>
#include <map>
#include <string>
#include <ros/ros.h>
#include <pegasus/VehicleState.h>
#include <pegasus/state_t.h>
#include <pegasus/DesiredControl.h>

#include <rosflight_msgs/RCRaw.h>

namespace apollo
{
class Controller
{
public:
  Controller();
  ~Controller();
private:
  //********************* NODE HANDLES *********************//
  ros::NodeHandle nh_;         // public node handle for publishing, subscribing

  //************** SUBSCRIBERS AND PUBLISHERS **************//
  ros::Subscriber vehicle_state_subscriber_;
  ros::Subscriber rx_subscriber_;

  ros::Publisher desired_command_publisher_;

  //******************** CLASS VARIABLES *******************//
protected:
  state_t state_;
  ros::Time last_time_;

  float roll_desired_;
  float pitch_desired_;
  float yaw_rate_desired_;
  float piD180_;

private:
  // rx Channel Variables
  int rx_[8];
  int aileron_stick_;
  int elevator_stick_;
  int throttle_stick_;
  int rudder_stick_;
  int aux1_stick_;
  int aux2_stick_;
  int aux3_stick_;
  int aux4_stick_;
  std::map<int, float> A_angle_map_;
  std::map<int, float> E_angle_map_;
  std::map<int, float> R_rate_map_;
  std::map<int, float> T_map_;
  int A_mid_us_;
  int E_mid_us_;

  int A_channel_;
  int E_channel_;
  int T_channel_;
  int R_channel_;

  //***************** CALLBACKS AND TIMERS *****************//
  void vehicleStateCallback(const VehicleStateConstPtr &msg);
  void rxCallback(const rosflight_msgs::RCRaw &msg);

protected:
  virtual void control(const ros::TimerEvent& event);
  ros::Timer control_timer_;
  //********************** FUNCTIONS ***********************//
  void pullParameters();
  void publishMomentCommands();
  void publishDesiredCommand();
  void mapControlChannels();
  void buildStickMap();
  void setChannels(std::string channel_map);
  void getRosParam(std::string parameter_name,      int    &param);
  void getRosParam(std::string parameter_name,      float  &param);
  void getRosParam(std::string parameter_name,      double &param);
  void getRosParam(std::string parameter_name, std::string &param);

};// end class Controller
} // end namespace apollo

#endif // CONTROLLER_BASE_H
