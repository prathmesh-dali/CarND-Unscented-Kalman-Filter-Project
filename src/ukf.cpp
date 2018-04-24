#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  H_Laser = MatrixXd::Zero(2,5);
  H_Laser(0,0) = 1;
  H_Laser(1,1) = 1;
  R_Laser = MatrixXd::Zero(2,2);
  R_Laser(0,0) = std_laspx_*std_laspx_;
  R_Laser(1,1) = std_laspy_*std_laspy_;

  R_Radar = MatrixXd::Zero(3,3);
  R_Radar<<std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      P_ <<std_radr_*std_radr_, 0, 0, 0, 0,
			  0, std_radr_*std_radr_, 0, 0, 0,
			  0, 0, 1.0, 0, 0,
			  0, 0, 0, 1.0, 0,
        0, 0, 0, 0, 1.0;
      x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]), meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]), 3, meas_package.raw_measurements_[1], 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      P_ <<std_laspx_*std_laspx_, 0, 0, 0, 0,
			  0, std_laspx_*std_laspx_, 0, 0, 0,
			  0, 0, 1.0, 0, 0,
			  0, 0, 0, 1.0, 0,
        0, 0, 0, 0, 1.0;
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 3, 0.5, 0;
    }
    // done initializing, no need to predict or update
  
		time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    weights_(0) = lambda_/(lambda_+n_aug_);
    for(int i=1; i<2 * n_aug_ + 1; i++){
      weights_(i) = 1/(2*(lambda_+n_aug_));
    }
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

  Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  // cout << "Xsig_pred_ = " << Xsig_pred_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  //create augmented mean vector
  VectorXd x_aug_ = VectorXd::Zero(7);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd::Zero(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  x_aug_.head(x_.size()) = x_;
  P_aug_.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug_(P_aug_.rows()-2, P_aug_.cols()-2) = std_a_*std_a_;
  P_aug_(P_aug_.rows()-1, P_aug_.cols()-1) = std_yawdd_*std_yawdd_;

  
  
  MatrixXd A = P_aug_.llt().matrixL();
  
  Xsig_aug_.block(0,0,n_aug_,1) = x_aug_;
  Xsig_aug_.block(0,1,n_aug_,n_aug_) = (sqrt(lambda_+n_aug_)*A).colwise()+x_aug_;
  Xsig_aug_.block(0,n_aug_+1,n_aug_,n_aug_) = (-sqrt(lambda_+n_aug_)*A).colwise()+x_aug_;

  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x_ = Xsig_aug_(0,i);
    double p_y_ = Xsig_aug_(1,i);
    double v_ = Xsig_aug_(2,i);
    double yaw_ = Xsig_aug_(3,i);
    double yawd_ = Xsig_aug_(4,i);
    double nu_a_ = Xsig_aug_(5,i);
    double nu_yawdd_ = Xsig_aug_(6,i);

    //predicted state values
    double px_p_, py_p_;

    //avoid division by zero
    if (fabs(yawd_) > 0.001) {
        px_p_ = p_x_ + v_/yawd_ * ( sin (yaw_ + yawd_*delta_t) - sin(yaw_));
        py_p_ = p_y_ + v_/yawd_ * ( cos(yaw_) - cos(yaw_+yawd_*delta_t) );
    }
    else {
        px_p_ = p_x_ + v_*delta_t*cos(yaw_);
        py_p_ = p_y_ + v_*delta_t*sin(yaw_);
    }

    double v_p_ = v_;
    double yaw_p_ = yaw_ + yawd_*delta_t;
    double yawd_p_ = yawd_;

    //add noise
    px_p_ = px_p_ + 0.5*nu_a_*delta_t*delta_t * cos(yaw_);
    py_p_ = py_p_ + 0.5*nu_a_*delta_t*delta_t * sin(yaw_);
    v_p_ = v_p_ + nu_a_*delta_t;

    yaw_p_ = yaw_p_ + 0.5*nu_yawdd_*delta_t*delta_t;
    yawd_p_ = yawd_p_ + nu_yawdd_*delta_t;
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p_;
    Xsig_pred_(1,i) = py_p_;
    Xsig_pred_(2,i) = v_p_;
    Xsig_pred_(3,i) = yaw_p_;
    Xsig_pred_(4,i) = yawd_p_;
  }

  x_= (Xsig_pred_.array().rowwise() * weights_.transpose().array()).rowwise().sum();
  MatrixXd intermediate_ = Xsig_pred_.colwise()-x_;
  for(int i=0; i<2 * n_aug_+ 1; i++ ){
     intermediate_.col(i)[3]= atan2(sin(intermediate_.col(i)(3)),cos(intermediate_.col(i)(3))); 
  }
  MatrixXd intermediateT_ = intermediate_.transpose();
  
  MatrixXd rowMulti_ = (intermediate_.array().rowwise() * weights_.transpose().array());
  P_ = rowMulti_*intermediateT_;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  VectorXd z_pred = H_Laser * x_;  
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H_Laser.transpose();
	MatrixXd S = H_Laser * P_ * Ht + R_Laser;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - (K * H_Laser)) * P_;
  NIS_L = y.transpose()*Si*y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  double rho, phi, rho_dot;
  for(int i=0; i<2 * n_aug_ + 1; i++ ){
    VectorXd data = Xsig_pred_.col(i);
    rho = sqrt(data[0]*data[0] + data[1]*data[1]);
    phi = atan2(data[1], data[0]);
    rho_dot = rho>0.00001?(data[0]* cos(data[3])*data[2]+data[1]* sin(data[3])*data[2])/rho:0.0;
    Zsig.col(i)<<rho,phi,rho_dot;
  }
  z_pred.fill(0.0);
  z_pred = (Zsig.array().rowwise() * weights_.transpose().array()).rowwise().sum();
  MatrixXd z_diff = Zsig.colwise()-z_pred;
  for(int i=0; i<2 * n_aug_ + 1; i++ ){
    z_diff.col(i)[1]= atan2(sin(z_diff.col(i)(1)),cos(z_diff.col(i)(1))); 
  }
  
  MatrixXd z_diff_transpose = z_diff.transpose();
  MatrixXd weight_z_diff = (z_diff.array().rowwise() * weights_.transpose().array());
  S.fill(0.0);
  S = weight_z_diff*z_diff_transpose+R_Radar;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd x_diff = Xsig_pred_.colwise()-x_;
  for(int i=0; i<2 * n_aug_ + 1; i++ ){
     x_diff.col(i)[3]= atan2(sin(x_diff.col(i)(3)),cos(x_diff.col(i)(3))); 
  }
  MatrixXd weight_x_diff = (x_diff.array().rowwise() * weights_.transpose().array());
  Tc.fill(0.0);
  Tc = Tc + weight_x_diff*z_diff_transpose;
  MatrixXd K = Tc*S.inverse();
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();
  NIS_R = y.transpose()*S.inverse()*y;
}
