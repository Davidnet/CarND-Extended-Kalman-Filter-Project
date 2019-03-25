#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // Sanity checks
  if (estimations.size() == 0) {
     return rmse;
  }

  VectorXd residuals(4);
  VectorXd residuals_sq(4);

  // accumulating squared residuals 
  for (int i = 0; i < estimations.size(); ++i) {
     residuals = estimations[i] - ground_truth[i];
     residuals_sq = residuals.array().square();
     rmse += residuals_sq;
  }

  // mean
  rmse = rmse / estimations.size();

  // squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
   MatrixXd Hj(3, 4);
   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);

   float den = px * px + py * py;

   if(fabs(den) < 0.0001) {

      Hj << 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;

      return Hj;
   }

	float inv_den = 1 / den;
	float denom_sqrt = sqrt(inv_den);

	//compute the Jacobian matrix
	Hj << px * denom_sqrt, py * denom_sqrt, 0, 0,
        -py * inv_den, px * inv_den, 0, 0,
		   py * (vx*py - vy*px) * pow(denom_sqrt, 3), px * (vy*px - vx*py) * pow(denom_sqrt, 3), px * denom_sqrt, py * denom_sqrt;

	return Hj;

}
