#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * fn that predict the state
   */

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * fn that updates the state by using Kalman Filter equations
   */

  VectorXd z_predict = H_ * x_;
  VectorXd y = z - z_predict;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = PHt * Sinv;

  // Estimation part
  x_ = x_ + (K * y);
  int i_size = x_.size();
  MatrixXd I = MatrixXd::Identity(i_size, i_size);

  // Update
  P_ -= K * H_ * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Numerical stability return in case
  float den = sqrt(px*px + py*py);
  if ( abs(den) < 0.0001) {
    return;
  }

  // compute H(x)
  VectorXd hx = VectorXd(3);
  hx << den,
        atan2(py, px),
        (px*vx + py*vy) / den;

  VectorXd y = z - hx;
  float phi = y(1);
  y(1) = atan2(sin(phi), cos(phi));

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = PHt * Sinv;

  // Estimation part
  x_ = x_ + (K * y);
  int i_size = x_.size();
  MatrixXd I = MatrixXd::Identity(i_size, i_size);

  // Update
  P_ -= K * H_ * P_;


}
