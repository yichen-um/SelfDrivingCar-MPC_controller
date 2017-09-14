#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.07;

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
const double Lf = 2.67;
double ref_v = 70;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

// Define FG object
class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  vector<double> previous_actuations;

  FG_eval(Eigen::VectorXd coeffs, vector<double> previous_actuations) { 
  	this->coeffs = coeffs;
  	this->previous_actuations = previous_actuations;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    // `fg` a vector of the cost constraints. Its first location is the cost, the following are the values states and actuations. The values of the constrains for the states are bounded by the state equations (vehicle dynamics, etc), and the actuation values are bounded by the performance of the actuators
      
    // `vars` is a vector of variable values (state & actuators). It contains all the current and model predicted states and control authorities in the future N steps.
      
      // Part 1, set up cost to be the first element of fg.
      fg[0] = 0;
      
      // Add the cost based on reference track
      for (int t = 0; t < N; t++) {
          fg[0] += CppAD::pow(vars[cte_start + t], 2);
          fg[0] += CppAD::pow(vars[epsi_start + t], 2);
          fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
      }
      // Add the cost based on control authorities
      for (int t = 0; t < N - 1; t++) {
          fg[0] += 1 * CppAD::pow(vars[delta_start + t], 2);
          fg[0] += 5 * CppAD::pow(vars[a_start + t], 2);
      }
      // Add the cost based on the change of control authorities
      for (int t = 0; t < N - 2; t++) {
          fg[0] += 2500 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
          fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
      }
      
      // Part 2, set up constrains, the first step in each segment is defined as the current states
      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];
      
      // the rest steps in each segment
      for (int t = 1; t < N; t++) {
          // The state at time t + 1
          AD<double> x1 = vars[x_start + t];
          AD<double> y1 = vars[y_start + t];
          AD<double> psi1 = vars[psi_start + t];
          AD<double> v1 = vars[v_start + t];
          AD<double> cte1 = vars[cte_start + t];
          AD<double> epsi1 = vars[epsi_start + t];
          
          // The state at time t
          AD<double> x0 = vars[x_start + t - 1];
          AD<double> y0 = vars[y_start + t - 1];
          AD<double> psi0 = vars[psi_start + t - 1];
          AD<double> v0 = vars[v_start + t - 1];
          AD<double> cte0 = vars[cte_start + t - 1];
          AD<double> epsi0 = vars[epsi_start + t - 1];
          
          // The actuation at time t
          AD<double> delta0 = vars[delta_start + t - 1];
          AD<double> a0 = vars[a_start + t - 1];
          
          // if (t > 1) {   // use previous actuations (to account for latency)
          //     a0 = vars[a_start + t - 2];
          //     delta0 = vars[delta_start + t - 2];
          // }
          
          // The reference trajectory
          AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
          AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));
          
          // The constrains in the rest steps are defined by the state esquations
          fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
          fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
          fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
          fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
          fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
          fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
      }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

Solution MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  // size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
    // Parse states
    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];
    
  // Set the number of steps prediction in each step
    size_t n_vars = N * 6 + (N - 1) * 2;
    
  // Set the number of constraints
    size_t n_constraints = N * 6;

  // Initial value of the independent variables. SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

  // Part 1, Set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
    // For non-actuators, the upper and lowerlimits are set to be the max negative and positive values. In such setting, the states are not bounded
    for (int i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }
    
    // For sterring angle, the upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    for (int i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

      // constraint delta to be the previous control with specified latency time
  	for (int i = delta_start; i < delta_start + latency_ind; i++) {
    	vars_lowerbound[i] = delta_prev;
    	vars_upperbound[i] = delta_prev;
  	}
    
    // For acceleration/decceleration upper and lower limits.
    for (int i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }
    
	// constaint acc to be the previous control with specified latency time
  	for (int i = a_start; i < a_start + latency_ind; i++) {
    	vars_lowerbound[i] = a_prev;
    	vars_upperbound[i] = a_prev;     
  	}
  // Part 2, Set lower and upper limits for the constraints. At the positions of the initial state it should assigned with value of initial states. At other positions, it should be 0 so that the state equations are enforced.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
    
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;
    
    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;
    
  // object that computes objective and constraints
  vector<double> previous_actuations = {delta_prev, a_prev};
  FG_eval fg_eval(coeffs, previous_actuations);

  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
    Solution sol;
	for (auto i = 0; i < N-1 ; i++) {
		cout << i << ": " << "solution.x[x_start+i]: " << solution.x[x_start+i] << "solution.x[y_start+i]: " << solution.x[y_start+i] << endl;
	    sol.X.push_back(solution.x[x_start+i]);
	    sol.Y.push_back(solution.x[y_start+i]);
	    sol.Delta.push_back(solution.x[delta_start+i]);
	    sol.A.push_back(solution.x[a_start+i]);
	}
  return sol;
}







