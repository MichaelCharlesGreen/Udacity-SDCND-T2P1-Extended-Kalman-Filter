#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

// Prediction is the same for both laser and radar. // L5M13
void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_; // L5M8
  MatrixXd Ft = F_.transpose(); // L5M9
  P_ = F_ * P_ * Ft + Q_;
}

// Update for laser. // L5M13
void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // L5M7
  // 29:18
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
    
  // The new estimate.
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

// Helper function...
// From: Project: Extended Kalman Filters M8 - Tips and Tricks:
// atan2() returns values between -pi and pi. When calculating phi in y = z - h(x) for
// radar measurements, the resulting angle phi in the y vector should be adjusted so that it is
// between -pi and pi. The Kalman filter is expecting small angle values between the range -pi and pi.
// HINT: when working in radians, you can add 2pi or subtract 2pi until the angle is within the desired
// range.
VectorXd NormalizeAngle(VectorXd &y) {
    std::cout << "y(1):" << y(1) << std::endl;
    if (y(1) > M_PI)
        y(1) = fmod(y(1), M_PI);
    if (y(1) < - M_PI)
        y(1) = -fmod(y(1), M_PI);
    
    return y;
};

// Update for radar. // L5M14 & L5M7

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    
    float px = x_(0); // position x
    float py = x_(1); // position y
    float vx = x_(2); // velocity x
    float vy = x_(3); // velocity y
    
    float ro = 0.0;
    float phi = 0.0;
    float ro_dot = 0.0;
    
    VectorXd z_pred(3);
    
    if (fabs(px) > 0.0001 && fabs(py) > 0.0001) {
        ro = sqrt(px * px + py * py);
        phi = atan2(py, px);
        ro_dot = (px * vx + py * vy) / ro;
    }
    
    z_pred << ro, phi, ro_dot;
    VectorXd y = z - z_pred;
    
    if (fabs(y(1)) > M_PI) {
        y = NormalizeAngle(y);
    }
    
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    
    // New estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
