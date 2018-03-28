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
  is_initialized_ = false;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  n_x_ = 5;
  
  n_aug_ = 7;
  
  lambda_ = 3 - n_aug_;
  
  no_of_sigma_points_ = 2 * n_aug_ + 1;
  
  x_aug_ = VectorXd(n_aug_);
  
  P_aug_ = MatrixXd(n_aug_,n_aug_);
  
  Xsig_pred_ = MatrixXd(n_x_, no_of_sigma_points_);
  
  weights_ = VectorXd(no_of_sigma_points_);
  
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_,0,0,
			  0,std_radphi_*std_radphi_,0,
			  0,0,std_radrd_*std_radrd_;
  
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_*std_laspx_,0,
			  0,std_laspy_*std_laspy_;
			  
  weights_(0)=lambda_/(lambda_+n_aug_);
  for (int i=1;i<no_of_sigma_points_;i++) {
      weights_(i) = 1 / (2*(lambda_+n_aug_));
  }
  
  H_laser_ = MatrixXd(2,n_x_);
  
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
	
    // first measurement
    cout << "UKF: " << endl;
    x_ << 0, 0, 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  x_(0) = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
	  x_(1) = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  x_(0) = meas_package.raw_measurements_[0];
	  x_(1) = meas_package.raw_measurements_[1];
    }
	
    P_ << 1, 0, 0, 0, 0,
		  0, 1, 0, 0, 0,
		  0, 0, 10, 0, 0,
		  0, 0, 0, 50, 0,
		  0, 0, 0, 0, 3;

    time_us_ = meas_package.timestamp_;
	
    // done initializing, no need to predict or update
    is_initialized_ = true;
	
	// print the output
    cout <<"-------\n";
    cout << "x_ = \n" << x_ << endl;
    cout << "P_ = \n" << P_ << endl;
	
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; 
  time_us_ = meas_package.timestamp_;
  
  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
   
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_)) {
    // Radar updates
	UpdateRadar(meas_package);
  } 
  else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_)) {
	// Laser updates
    UpdateLidar(meas_package);
  }
  
  // print the output
  cout <<"-------\n";
  cout << "x_ = \n" << x_ << endl;
  cout << "P_ = \n" << P_ << endl; 
  
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
  // Generating the sigma points
  x_aug_.head(5) = x_ ;
  x_aug_(5)=0;
  x_aug_(6)=0;
  
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  MatrixXd A = P_aug_.llt().matrixL();
  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, no_of_sigma_points_);
  
  Xsig_aug.col(0) = x_aug_;
  for (int i=0; i < n_aug_; i++) {
      Xsig_aug.col(i+1) = x_aug_ + sqrt(lambda_+n_aug_) * A.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * A.col(i);
  }
    
  //Predict the sigma points by passing them through the process model
  for(int i=0 ; i<no_of_sigma_points_; i++) {
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double psi = Xsig_aug(3,i);
      double psi_dot = Xsig_aug(4,i);
      double noise_a = Xsig_aug(5,i);
      double noise_psi_dot_dot = Xsig_aug(6,i);
      
      if (fabs(psi_dot) < 0.0001) {
        Xsig_pred_(0,i) = px + v*cos(psi)*delta_t + delta_t*delta_t*cos(psi)*noise_a/2;
        Xsig_pred_(1,i) = py + v*sin(psi)*delta_t + delta_t*delta_t*sin(psi)*noise_a/2;
      }
      else{
        Xsig_pred_(0,i) = px + (v*(sin(psi+psi_dot*delta_t)-sin(psi)))/psi_dot + delta_t*delta_t*cos(psi)*noise_a/2;
        Xsig_pred_(1,i) = py + (v*(-1*cos(psi+psi_dot*delta_t)+cos(psi)))/psi_dot + delta_t*delta_t*sin(psi)*noise_a/2;
      }
	  
	  Xsig_pred_(2,i) = v + delta_t*noise_a;
      Xsig_pred_(3,i) = psi + psi_dot*delta_t +delta_t*delta_t*noise_psi_dot_dot/2;
	  Xsig_pred_(4,i) = psi_dot +delta_t*noise_psi_dot_dot;
  }
  
  //Predict mean and covariance matrix using predicted sigma points  
  x_.fill(0.0);
  for(int i = 0; i<no_of_sigma_points_ ; i++) {
      x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for(int i = 0; i<no_of_sigma_points_ ; i++) {
      VectorXd diff = Xsig_pred_.col(i)-x_ ;
      diff(3)=atan2(sin(diff(3)),cos(diff(3)));
      P_ = P_ + weights_(i)*diff*diff.transpose();
  }
  
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
  H_laser_ << 1, 0, 0, 0, 0,
	          0, 1, 0, 0, 0;  
  
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;

  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;

  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
  
  //Calculating NIS
  NIS_laser_ = y.transpose() * S.inverse() * y;
  
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
  //measurement dimension
  int n_z = 3;
  
  //matrix to store sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, no_of_sigma_points_);
  
  //matrix to store predicted sigma points
  VectorXd z_pred = VectorXd(n_z);
  
  //covariance matrix with predicted covariances
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //converting predicted state space to measurement space
  for (int i=0; i<no_of_sigma_points_;i++) {
      float rho = sqrt(Xsig_pred_(0,i)*Xsig_pred_(0,i) + Xsig_pred_(1,i)*Xsig_pred_(1,i));
      float phi = atan2(Xsig_pred_(1,i),Xsig_pred_(0,i));
      float rho_dot;
      
      if(fabs(rho)<0.0001) {
          if(rho>0) {
			  rho += 0.0001;
		  }
		  else {
			  rho -= 0.0001;
		  }
      }
      
      rho_dot = (Xsig_pred_(0,i)*cos(Xsig_pred_(3,i))*Xsig_pred_(2,i) + Xsig_pred_(1,i)*sin(Xsig_pred_(3,i))*Xsig_pred_(2,i))/rho;

      Zsig(0,i) = rho;
	  Zsig(1,i) = phi;
	  Zsig(2,i) = rho_dot;
  }
  
  z_pred.fill(0.0);
  for (int i=0; i< no_of_sigma_points_;i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  
  S.fill(0.0);   
  for (int i=0; i < no_of_sigma_points_; i++) {
    VectorXd diff = Zsig.col(i) - z_pred;
    diff(1) = atan2(sin(diff(1)),cos(diff(1)));
          
    S = S + weights_(i)*diff*diff.transpose();
  }
  
  S = S + R_radar_ ;
  
  //Matrix to store cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //Update state and covariance matrices based on measurements  
  Tc.fill(0.0);
  for(int i=0; i<no_of_sigma_points_;i++) {
      VectorXd diff_state = Xsig_pred_.col(i) - x_ ;
      diff_state(3) = atan2(sin(diff_state(3)),cos(diff_state(3)));
      
      VectorXd diff_meas = Zsig.col(i) - z_pred ;
      diff_meas(1) = atan2(sin(diff_meas(1)),cos(diff_meas(1)));
      
      Tc = Tc + weights_(i) * diff_state * diff_meas.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();
  
  VectorXd diff_meas = meas_package.raw_measurements_ - z_pred;
  diff_meas(1) = atan2(sin(diff_meas(1)),cos(diff_meas(1)));  
  
  x_ = x_ + K * diff_meas;
  
  P_ = P_ - K * S * K.transpose();
  
  
  //Calculating NIS
  NIS_radar_ = diff_meas.transpose() * S.inverse() * diff_meas;
  
}
