#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  const Eigen::VectorXd& coeffs;
  const Offsets& o;
  const Params& p;
  FG_eval(const Eigen::VectorXd& coeffs_, const Offsets& o_, const Params& p_) : coeffs(coeffs_), o(o_), p(p_) {}

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    // fg a vector of constraints, x is a vector of constraints.
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    auto exp_wt = 1.0;
    for (int i = 0; i < o.N; i++) {
      fg[0] += CppAD::pow(vars[o.cte_start + i] - p.ref_cte, 2)   * p.w_cte   * exp_wt;
      fg[0] += CppAD::pow(vars[o.epsi_start + i] - p.ref_epsi, 2) * p.w_epsi  * exp_wt;
      fg[0] += CppAD::pow(vars[o.v_start + i] - p.ref_v, 2)       * p.w_v     * exp_wt;
      exp_wt *= p.w_exp;
    }
    // Minimize the use of actuators.
    exp_wt = 1.0;
    for (int i = p.projected; i < o.N - 1; i++) {
      fg[0] += CppAD::pow(vars[o.delta_start + i], 2) * p.w_delta * exp_wt;
      fg[0] += CppAD::pow(vars[o.a_start + i], 2)     * p.w_a     * exp_wt;
      exp_wt *= p.w_exp;
    }
    // Minimize the value gap between sequential actuations.
    exp_wt = 1.0;
    for (int i = (p.projected >0 ? p.projected-1 : 0); i < o.N - 2; i++) {
      fg[0] += CppAD::pow(vars[o.delta_start + i + 1] - vars[o.delta_start + i], 2) * p.w_delta2  * exp_wt; 
      fg[0] += CppAD::pow(vars[o.a_start + i + 1] - vars[o.a_start + i], 2)         * p.w_a2      * exp_wt; 
      exp_wt *= p.w_exp;
    }

    // Setup Model Constraints

    // Initial constraints
    // We add 1 to each of the starting indices due to cost being located at index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + o.x_start] = vars[o.x_start];
    fg[1 + o.y_start] = vars[o.y_start];
    fg[1 + o.psi_start] = vars[o.psi_start];
    fg[1 + o.v_start] = vars[o.v_start];
    fg[1 + o.cte_start] = vars[o.cte_start];
    fg[1 + o.epsi_start] = vars[o.epsi_start];

    // The rest of the constraints // i.e. the state update equations for the model
    for (int i = 0; i < o.N - 1; i++) {
      // The state at time t+1 .
      AD<double> x1 = vars[o.x_start + i + 1];
      AD<double> y1 = vars[o.y_start + i + 1];
      AD<double> psi1 = vars[o.psi_start + i + 1];
      AD<double> v1 = vars[o.v_start + i + 1];
      AD<double> cte1 = vars[o.cte_start + i + 1];
      AD<double> epsi1 = vars[o.epsi_start + i + 1];

      // The state at time t.
      AD<double> x0 = vars[o.x_start + i];
      AD<double> y0 = vars[o.y_start + i];
      AD<double> psi0 = vars[o.psi_start + i];
      AD<double> v0 = vars[o.v_start + i];
      AD<double> cte0 = vars[o.cte_start + i];
      AD<double> epsi0 = vars[o.epsi_start + i];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[o.delta_start + i];
      AD<double> a0 = vars[o.a_start + i];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0 + 3*coeffs[3]*x0*x0);

      // Setup the rest of the model constraints
      fg[2 + o.x_start + i]    = x1     - (x0 + v0 * CppAD::cos(psi0) * p.dt);  // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      fg[2 + o.y_start + i]    = y1     - (y0 + v0 * CppAD::sin(psi0) * p.dt);  // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      fg[2 + o.psi_start + i]  = psi1   - (psi0 - v0 * delta0 / Lf * p.dt);     // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      fg[2 + o.v_start + i]    = v1     - (v0 + a0 * p.dt);                     // v_[t+1] = v[t] + a[t] * dt
      fg[2 + o.cte_start + i]  = cte1   - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * p.dt));  // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      fg[2 + o.epsi_start + i] = epsi1  - ((psi0 - psides0) + v0 * delta0 / Lf * p.dt);   // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC(int N_, const Params& p_) : o(N_), p(p_) {
  mpc_x_vals.reserve(N_ / p.MPC_pts_offsets + 1);
  mpc_y_vals.reserve(N_ / p.MPC_pts_offsets + 1);
  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  // Uncomment this if you'd like more print information
  solver_options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  solver_options += "Sparse  true        forward\n";
  solver_options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  solver_options += "Numeric max_cpu_time          0.5\n";
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double delta, double a) {
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  size_t n_vars = o.N * state.size() + 2*(o.N - 1);   // number of model variables (includes both states and inputs).
  size_t n_constraints = o.N * state.size();          // number of constraints

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  vars[o.x_start] = x;
  vars[o.y_start] = y;
  vars[o.psi_start] = psi;
  vars[o.v_start] = v;
  vars[o.cte_start] = cte;
  vars[o.epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < o.delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  static const auto STEERING_LIMIT = deg2rad(25.0);
  for (int i = o.delta_start, j=0; i < o.a_start; i++, ++j) {
    if (j < p.projected) { // constrain until latency window
      vars_lowerbound[i] = vars_upperbound[i] = delta; 
      continue;
    }
    // steering limit further constrained based on current speed
    vars_lowerbound[i] = -STEERING_LIMIT; ///(v/25); //- 0.436332;
    vars_upperbound[i] =  STEERING_LIMIT; ///(v/25); //0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = o.a_start, j=0; i < n_vars; i++, ++j) {
    if (j < p.projected) { // constrain until latency window
      vars_lowerbound[i] = vars_upperbound[i] = a;
      continue;
    }
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[o.x_start] = x;
  constraints_lowerbound[o.y_start] = y;
  constraints_lowerbound[o.psi_start] = psi;
  constraints_lowerbound[o.v_start] = v;
  constraints_lowerbound[o.cte_start] = cte;
  constraints_lowerbound[o.epsi_start] = epsi;

  constraints_upperbound[o.x_start] = x;
  constraints_upperbound[o.y_start] = y;
  constraints_upperbound[o.psi_start] = psi;
  constraints_upperbound[o.v_start] = v;
  constraints_upperbound[o.cte_start] = cte;
  constraints_upperbound[o.epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, o, p);

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      solver_options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  mpc_x_vals.clear(); mpc_y_vals.clear();
  for (i = p.projected; i < o.N; i = i + p.MPC_pts_offsets)
  {
    mpc_x_vals.push_back(solution.x[o.x_start + i]);
    mpc_y_vals.push_back(solution.x[o.y_start + i]);
  }
  
  return {solution.x[o.delta_start + p.projected], solution.x[o.a_start + p.projected],
          solution.x[o.x_start + p.projected], solution.x[o.y_start + p.projected],
          solution.x[o.psi_start + p.projected], solution.x[o.v_start + p.projected],
          solution.x[o.cte_start + p.projected], solution.x[o.epsi_start + p.projected],
          };
}
