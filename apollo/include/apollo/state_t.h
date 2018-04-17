#ifndef STATE_T_H
#define STATE_T_H

#include <apollo/VehicleState.h>


namespace apollo
{
  struct state_t
  {
    //*********************** VARIABLE ***********************//
    // EDIT EVERY FUNCTION IF YOU ADD A VARIABLE TO THE STRUCT OR MESSAGE
    float phi;       // Roll, angle rotated around i^(vehicle-2)
    float theta;     // Pitch, angle rotated around j^(vehicle-1)
    float psi;       // Yaw, angle rotated around k^(vehicle)
    float p;         // Roll rate around i^(body)
    float q;         // Pitch rate around j^(body)
    float r;         // Yaw rate around k^(body)

  private:
    VehicleState _msg_;
  public:
    state_t()
    {
      phi   = 0.0f;
      theta = 0.0f;
      psi   = 0.0f;
      p     = 0.0f;
      q     = 0.0f;
      r     = 0.0f;
    }
    //****************** OPERATOR FUNCTIONS ******************//
    state_t operator+(const state_t s)
    {
      state_t n;
      n.phi   = s.phi   + phi;
      n.theta = s.theta + theta;
      n.psi   = s.psi   + psi;
      n.p     = s.p     + p;
      n.q     = s.q     + q;
      n.r     = s.r     + r;
      return n;
    }
    state_t operator*(const float num)
    {
      state_t n;
      n.phi   = phi   * num;
      n.theta = theta * num;
      n.psi   = psi   * num;
      n.p     = p     * num;
      n.q     = q     * num;
      n.r     = r     * num;
      return n;
    }
    //********************** FUNCTIONS ***********************//
    VehicleState msg()
    {
      _msg_.phi   = phi;
      _msg_.theta = theta;
      _msg_.psi   = psi;
      _msg_.p     = p;
      _msg_.q     = q;
      _msg_.r     = r;
      return _msg_;
    }
    void msg2struct(const VehicleState msg_in)
    {
      phi   = msg_in.phi;
      theta = msg_in.theta;
      psi   = msg_in.psi;
      p     = msg_in.p;
      q     = msg_in.q;
      r     = msg_in.r;
    }
    void msg2struct(const VehicleStateConstPtr &msg_in)
    {
      phi   = msg_in->phi;
      theta = msg_in->theta;
      psi   = msg_in->psi;
      p     = msg_in->p;
      q     = msg_in->q;
      r     = msg_in->r;
    }
  };
} // end namespace apollo

#endif // STATE_T_H
