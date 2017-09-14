#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

const int latency_ind = 2; // Latency 

struct Solution {
	vector<double> X;
	vector<double> Y;
	vector<double> Psi;
	vector<double> V;
	vector<double> Cte;
	vector<double> Epsi;
	vector<double> Delta;
	vector<double> A;
	};

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  double delta_prev {0};
  double a_prev {0.1};
};

#endif /* MPC_H */
