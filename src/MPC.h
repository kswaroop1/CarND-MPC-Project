#ifndef MPC_H
#define MPC_H

#define HAVE_CSTDDEF
#include <vector>
#include <string>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
constexpr double Lf = 2.67;

struct Offsets
{
  Offsets(size_t N_) {
    N = N_;
    x_start = 0;
    y_start = x_start + N;
    psi_start = y_start + N;
    v_start = psi_start + N;
    cte_start = v_start + N;
    epsi_start = cte_start + N;
    delta_start = epsi_start + N;
    a_start = delta_start + N - 1;
  }
  size_t N, x_start, y_start, psi_start, v_start, cte_start, epsi_start, delta_start, a_start;
};

struct Params
{
  Params(int MPC_pts_offsets_, double dt_, double ref_v_, double delay_) : MPC_pts_offsets(MPC_pts_offsets_), dt(dt_), ref_v(ref_v_), delay(delay_), ref_cte(0), ref_epsi(0) {  
    projected = int(delay_ / dt_); // for 100 ms delay => 0.1 sec; dt=0.05 => 0.1/0.05=2 steps in future
    w_v = 1; w_cte = w_epsi = 100; w_delta = 5000; w_a = 10; w_delta2 = 100; w_a2 = 10; w_exp = 0.9;
  }
  int MPC_pts_offsets, projected;
  double dt, ref_v, ref_cte, ref_epsi, delay;
  double w_v, w_cte, w_epsi, w_delta, w_a, w_delta2, w_a2, w_exp;
};


class MPC {
 public:
  MPC(int N, const Params& p);

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double detla, double a);
  string solver_options;
  vector<double> mpc_x_vals;
  vector<double> mpc_y_vals;
private:
  double dt, ref_v, delay;
  Offsets o;
  Params p;
};

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
inline double deg2rad(double x) { return x * pi() / 180; }
inline double rad2deg(double x) { return x * 180 / pi(); }

// Evaluate a polynomial.
inline double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Evaluate derivative of polynomial.
inline double polyevalderiv(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * pow(x, i-1);
  }
  return result;
}

#endif /* MPC_H */
