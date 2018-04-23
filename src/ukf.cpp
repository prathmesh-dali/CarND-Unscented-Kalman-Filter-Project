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
    P_ <<0.15, 0, 0, 0, 0,
			  0, 0.15, 0, 0, 0,
			  0, 0, 1.0, 0, 0,
			  0, 0, 0, 1.0, 0,
        0, 0, 0, 0, 1.0;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]), meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]), 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
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
  MatrixXd H_ = MatrixXd::Zero(2,5);
  H_(0,0) = 1;
  H_(1,1) = 1;
  MatrixXd R_ = MatrixXd::Zero(2,2);
  R_(0,0) = std_laspx_;
  R_(1,1) = std_laspy_;
  VectorXd z_pred = H_ * x_;  
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - (K * H_)) * P_;
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


  MatrixXd R = MatrixXd::Zero(3,3);
  R<<std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
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
  MatrixXd intermediate = Zsig.colwise()-z_pred;
  for(int i=0; i<2 * n_aug_ + 1; i++ ){
    VectorXd temp1 = intermediate.col(i);
    intermediate.col(i)[1]= atan2(sin(temp1(1)),cos(temp1(1))); 
  }
  
  MatrixXd intermediateT = intermediate.transpose();
  MatrixXd rowMulti = (intermediate.array().rowwise() * weights_.transpose().array());
  S.fill(0.0);
  S = rowMulti*intermediateT+R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd intermediate1 = Xsig_pred_.colwise()-x_;
  MatrixXd intermediate2 = (Zsig.colwise()-z_pred).transpose();
  for(int i=0; i<2 * n_aug_ + 1; i++ ){
     intermediate1.col(i)[3]= atan2(sin(intermediate1.col(i)(3)),cos(intermediate1.col(i)(3))); 
     intermediate2.row(i)[1]= atan2(sin(intermediate2.row(i)(1)),cos(intermediate2.row(i)(1))); 
  }
  MatrixXd intermediate3 = (intermediate1.array().rowwise() * weights_.transpose().array());
  Tc.fill(0.0);
  Tc = Tc + intermediate3*intermediate2;
  MatrixXd K = Tc*S.inverse();
  x_ = x_ + K*(meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K * S * K.transpose();
}
