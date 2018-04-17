#include <apollo/controller_base.h>

#include <apollo/pid.h>

namespace apollo
{
Controller::Controller() :
  nh_(ros::NodeHandle())
{
  piD180_ = M_PI/180.0f;
  //********************** PARAMETERS **********************//
  float control_rate;
  getRosParam("control_rate", control_rate);
  pullParameters();

  //************** SUBSCRIBERS AND PUBLISHERS **************//
  vehicle_state_subscriber_ = nh_.subscribe("truth",1,&Controller::vehicleStateCallback, this);
  rx_subscriber_ = nh_.subscribe("/rosflight/rc_raw",1,&Controller::rxCallback, this);

  moment_command_publisher_   = nh_.advertise<TODO>("motor_command",1);
  desired_command_publisher_  = nh_.advertise<DesiredControl>("desired_command",1);

  //******************** CLASS VARIABLES *******************//
  last_time_ = ros::Time::now();
  yaw_rate_desired_   = 0.0f;
  roll_desired_       = 0.0f;
  pitch_desired_      = 0.0f;

  //***************** CALLBACKS AND TIMERS *****************//
  control_timer_ = nh_.createTimer(ros::Duration(1.0/control_rate), &Controller::control, this);

  //********************** FUNCTIONS ***********************//
  buildStickMap();
}
Controller::~Controller()
{

}
void Controller::control(const ros::TimerEvent& event)
{
  ROS_ERROR("CHILD CLASS FUNCTION 'control' WAS NOT CALLED.");
}
void Controller::vehicleStateCallback(const VehicleStateConstPtr &msg)
{
  state_.msg2struct(msg);
}
void Controller::rxCallback(const rosflight_msgs::RCRaw &msg)
{
  rx_[0] = msg.values[0];
  rx_[1] = msg.values[1];
  rx_[2] = msg.values[2];
  rx_[3] = msg.values[3];
  rx_[4] = msg.values[4];
  rx_[5] = msg.values[5];
  rx_[6] = msg.values[6];
  rx_[7] = msg.values[7];
}
void Controller::buildStickMap()
{
  float max_angle, angle_expo, rate_expo, thrust_expo, rc, trc,  a_super_rate, r_super_rate, t_super_rate, max_height;
  int A_min_us, A_max_us, E_min_us, E_max_us, T_min_us, T_mid_us, T_max_us;
  int R_min_us, R_mid_us, R_max_us;
  getRosParam("rx/aileron/min_us", A_min_us);
  getRosParam("rx/aileron/mid_us",A_mid_us_);
  getRosParam("rx/aileron/max_us",A_max_us);
  getRosParam("rx/elevator/min_us",E_min_us);
  getRosParam("rx/elevator/mid_us",E_mid_us_);
  getRosParam("rx/elevator/max_us",E_max_us);
  getRosParam("rx/thrust/min_us",T_min_us);
  getRosParam("rx/thrust/mid_us",T_mid_us);
  getRosParam("rx/thrust/max_us",T_max_us);
  getRosParam("rx/rudder/min_us",R_min_us);
  getRosParam("rx/rudder/mid_us",R_mid_us);
  getRosParam("rx/rudder/max_us",R_max_us);
  getRosParam("rx/curve/angle_mode/max_angle",max_angle);
  max_angle = max_angle*piD180_;
  getRosParam("rx/curve/angle_mode/expo",angle_expo);
  getRosParam("rx/curve/rate_mode/expo",rate_expo);
  getRosParam("rx/curve/thrust/expo",thrust_expo);
  getRosParam("rx/curve/rate_mode/rc",rc);
  getRosParam("rx/curve/thrust/rc",trc);
  getRosParam("rx/curve/angle_mode/super_rate",a_super_rate);
  getRosParam("rx/curve/rate_mode/super_rate",r_super_rate);
  getRosParam("rx/curve/thrust/super_rate",t_super_rate);
  float shorter;
  if (A_max_us - A_mid_us_ < A_mid_us_ - A_min_us && A_max_us - A_mid_us_ < 500.0f)
    shorter = A_max_us - A_mid_us_;
  else if (A_mid_us_ - A_min_us < 500.0f)
    shorter = A_mid_us_ - A_min_us;
  else
    shorter = 500.0f;
  for (int i = A_min_us; i <= A_max_us; i++)
  {
    if (i == A_mid_us_)
    {
      A_angle_map_[i]  = 0.0f;
    }
    else
    {
    float mid_us_ = A_mid_us_;
    A_angle_map_[i]  = ( pow((abs(i-mid_us_)/shorter),pow(10.0,angle_expo))*max_angle\
                       + pow((abs(i-mid_us_)/shorter),pow(10.0, .8))*max_angle*0.4f*a_super_rate)\
                       *(i-mid_us_)/abs(i-mid_us_);
    }
  }
  if (E_max_us - E_mid_us_ < E_mid_us_ - E_min_us && E_max_us - E_mid_us_ < 500.0f)
    shorter = E_max_us - E_mid_us_;
  else if (E_mid_us_ - E_min_us < 500.0f)
    shorter = E_mid_us_ - E_min_us;
  else
    shorter = 500.0f;
  for (int i = E_min_us; i <= E_max_us; i++)
  {
    if (i == E_mid_us_)
    {
      E_angle_map_[i]  = 0.0f;
    }
    else
    {
    float mid_us_ = E_mid_us_;
    E_angle_map_[i]  = ( pow((abs(i-mid_us_)/shorter),pow(10.0,angle_expo))*max_angle\
                       + pow((abs(i-mid_us_)/shorter),pow(10.0, .8))*max_angle*0.4f*a_super_rate)\
                       *(mid_us_ - i)/abs(i-mid_us_);
    }
  }
  for (int i = T_min_us; i <= T_max_us; i++)
  {
    float mid_us_ = T_mid_us;
    T_map_[i] = ( pow((i - T_min_us)/((float) T_max_us - (float) T_min_us),pow(10.0,thrust_expo))*trc
                + pow((i - T_min_us)/((float) T_max_us - (float) T_min_us),pow(10.0, 0.8))*trc*0.4f*t_super_rate);
  }
  if (R_max_us - R_mid_us < R_mid_us - R_min_us && R_max_us - R_mid_us < 500.0f)
    shorter = R_max_us - R_mid_us;
  else if (R_mid_us - R_min_us < 500.0f)
    shorter = R_mid_us - R_min_us;
  else
    shorter = 500.0f;
  for (int i = R_min_us; i <= R_max_us; i++)
  {
    if (i == R_mid_us)
    {
      R_rate_map_[i]   = 0.0f;
    }
    else
    {
    float mid_us_ = R_mid_us;
    R_rate_map_[i]   = ( pow((abs(i-mid_us_)/shorter),pow(10.0,rate_expo))*400.0f*piD180_*rc\
                       + pow((abs(i-mid_us_)/shorter),pow(10.0, .8))*400.0f*piD180_*rc*0.4f*r_super_rate)\
                       *(i-mid_us_)/abs(i-mid_us_);
    }
  }
}
void Controller::mapControlChannels()
{
  aileron_stick_      = rx_[A_channel_];
  elevator_stick_     = rx_[E_channel_];
  throttle_stick_     = rx_[T_channel_];
  rudder_stick_       = rx_[R_channel_];

  roll_desired_       = A_angle_map_[aileron_stick_];
  pitch_desired_      = E_angle_map_[elevator_stick_];
  thrust_desired_     = T_map_[throttle_stick_];
  yaw_rate_desired_   = R_rate_map_[rudder_stick_];
}
void Controller::publishDesiredCommand()
{
  DesiredControl des_msg;
  des_msg.thrust_desired     = thrust_desired_;
  des_msg.yaw_rate_desired   = yaw_rate_desired_;
  des_msg.roll_desired       = roll_desired_;
  des_msg.pitch_desired      = pitch_desired_;
  desired_command_publisher_.publish(des_msg);
}
void Controller::publishMomentCommands()
{

}
void Controller::setChannels(std::string channel_map)
{
  for (unsigned int i = 0; i < 8; i++)
  {
    switch (channel_map[i])
    {
    case 'A':
      A_channel_ = i;
      break;
    case 'E':
      E_channel_ = i;
      break;
    case 'T':
      T_channel_ = i;
      break;
    case 'R':
      R_channel_ = i;
      break;
    case '1':
      aux1_channel_ = i;
      break;
    case '2':
      aux2_channel_ = i;
      break;
    case '3':
      aux3_channel_ = i;
      break;
    case '4':
      aux4_channel_ = i;
      break;
    }
  }
}
void Controller::pullParameters()
{
  std::string channel_map;
  getRosParam("rx/channel_map",channel_map);
  setChannels(channel_map);
}
void Controller::getRosParam(std::string parameter_name, int &param)
{
  if (!(ros::param::get(parameter_name ,param)))
    ROS_WARN("%s",("No param named '" + parameter_name + "'").c_str());
}
void Controller::getRosParam(std::string parameter_name, float &param)
{
  if (!(ros::param::get(parameter_name ,param)))
    ROS_WARN("%s",("No param named '" + parameter_name + "'").c_str());
}
void Controller::getRosParam(std::string parameter_name, double &param)
{
  if (!(ros::param::get(parameter_name ,param)))
    ROS_WARN("%s",("No param named '" + parameter_name + "'").c_str());
}
void Controller::getRosParam(std::string parameter_name, std::string &param)
{
  if (!(ros::param::get(parameter_name ,param)))
    ROS_WARN("%s",("No param named '" + parameter_name + "'").c_str());
}
} // end namespace apollo

//********************************************************//
//************************ MAIN **************************//
//********************************************************//
int main(int argc, char** argv)
{
  ros::init(argc, argv, "controller");
  ros::NodeHandle nh("apollo");

  pegasus::Controller *controller_obj;
  std::string control_type;
  if (!(ros::param::get("control_type",control_type)))
    ROS_WARN("No param named 'control_type'");
  if (control_type == "PID")
    controller_obj = new apollo::PID;
  else
    ROS_ERROR("NO CONTROLLER INITIALIZED");

  ros::spin();

  delete controller_obj;
  return 0;
} // end main
