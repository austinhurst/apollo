#ifndef PID_H
#define PID_H

#include <apollo/controller_base.h>

namespace apollo
{
class PID : public Controller
{
public:
  PID();

private:
  virtual void control(const ros::TimerEvent& event);
  double ts_;
  float mass_;

  // PHI
  float kP_phi_;
  float kD_phi_;

  // THETA
  float kP_theta_;
  float kD_theta_;

  // PSI
  float kP_psi_;
  float kD_psi_;
  float sigma_;
  float r_last_;
  float rd_;

  float saturate(float value_i, float min, float max);
};// end class PID
} // end namespace apollo

#endif // PID_H
