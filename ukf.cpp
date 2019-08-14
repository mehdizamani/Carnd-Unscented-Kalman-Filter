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
  is_initialized_ = false;
  
  
  
  
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  

 
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  
  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_, n_x_);
  x_ << 1.0, 1.0, 1.0, 1.0, 0.1;
	//create example covariance matrix
  P_ <<     1.0,    0.0,    0.0,    0.0,   0.0,
	         0.0,    1.0,    0.0,    0.0,   0.0,
		     0.0,    0.0,    1.0,    0.0,   0.0,
		     0.0,    0.0,    0.0,    1.0,   0.0,
		     0.0,    0.0,    0.0,    0.0,   1.0;
	
  x_aug_     = VectorXd(n_aug_);
  P_aug_     = MatrixXd(n_aug_, n_aug_);
  Xsig_aug_  = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  L_ = MatrixXd(n_aug_,n_aug_);
  A_ = MatrixXd(n_x_,n_x_);
   // set weights  NOTICE: since the weights do not change during the iterations, 
  // it is better to initiate them only once
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
	weights_(i) = weight;
  }
  nis_radar_list_ = list<double>();
  nis_lidar_list_ = list<double>();
  
  
  lidar_data_elem_ = 2;
  radar_data_elem_ = 3;
  R_Lidar_ = MatrixXd::Zero(lidar_data_elem_, lidar_data_elem_);
  R_Lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  R_Radar_ = MatrixXd::Zero(radar_data_elem_, radar_data_elem_);
  R_Radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
	 
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
    if (!is_initialized_) 
	{
		time_us_ = measurement_pack.timestamp_;
		cout<<time_us_<<endl;
		is_initialized_ = true;
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
		{
		  //Convert radar from polar to state space x,y,v,rhodot,rho
     		double rho = measurement_pack.raw_measurements_(0);
			double phi = measurement_pack.raw_measurements_(1);
			double rhodot = measurement_pack.raw_measurements_(2);
			// convrt to cartesian space
			double dv = sqrt (pow(rhodot*cos(phi),2)+pow(rhodot*sin(phi),2));
			double dyaw = 0.0;
			if (rhodot*cos(phi)!=0)
				dyaw = rhodot*sin(phi)/rhodot*cos(phi);
			x_ << rho * cos(phi),rho * sin(phi),dv,dyaw,0.0;
		} 
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
		{
			//Initialize state. set the state with the initial location and zero velocity
			x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0,0.0; 
		}
	}

	//compute the time elapsed between the current and previous measurements
	double dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
	//dt - expressed in seconds
	time_us_ = measurement_pack.timestamp_;
	Prediction(dt);
	if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	  UpdateRadar(measurement_pack);
	else if (use_laser_&& +measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
	  UpdateLidar(measurement_pack);
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
  
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  L_ = P_aug_.llt().matrixL();
  A_ = P_.llt().matrixL();
  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }
  
  double dt=delta_t;
  Xsig_pred_ = SigmaPointMotionPrediction(Xsig_aug_,dt);
 
 //predicted state mean
 // x_.fill(0.0);
 // for (int i = 0; i < 2 * n_aug_ + 1; i++)   
  //   x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  x_ = Xsig_pred_*weights_;
  
 //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
	double diff=x_diff(3);
	x_diff(3) = AngelNormalization(diff);
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
 
    
}

MatrixXd UKF::SigmaPointMotionPrediction(const MatrixXd Xsig_out,const double delt)
{
	MatrixXd ret = MatrixXd(n_x_,2*n_aug_+1); 
	ret.fill(0.0);
	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		double p_x 		= Xsig_out(0,i);
		double p_y 		= Xsig_out(1,i);
		double v 		= Xsig_out(2,i);
		double yaw		= Xsig_out(3,i);
		double yawd		= Xsig_out(4,i); 
		double nu_a 	= Xsig_out(5,i);
		double nu_yawdd = Xsig_out(6,i);
		double dt2 = delt * delt;
		//to avoid any division by zero
		if (fabs(yawd) > 0.001) {
			ret(0,i) = p_x + v/yawd * (sin(yaw + yawd * delt) - sin(yaw)) + (0.5 * dt2 * cos(yawd) * nu_a);
			ret(1,i) = p_y + v/yawd * (-cos(yaw + yawd * delt) + cos(yaw)) + (0.5 * dt2 * sin(yawd) * nu_a);
		}
		else {
			ret(0,i) = p_x + v*delt*cos(yaw) + (0.5 * dt2 * cos(yawd) * nu_a);
			ret(1,i) = p_y + v*delt*sin(yaw) + (0.5 * dt2 * sin(yawd) * nu_a);
		}

		ret(2,i) = v + nu_a * delt;
		ret(3,i) = yaw + yawd * delt + 0.5 * nu_yawdd * dt2;
		ret(4,i) = yawd + nu_yawdd * delt;
    }
	return ret;
}


double UKF::AngelNormalization(const double diff)
{
	double inut = diff;
	while (inut> M_PI) 
		inut -= 2.0*M_PI;
	while (inut<-M_PI)
			inut += 2.0*M_PI;
	return inut;
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
    VectorXd z = VectorXd(lidar_data_elem_);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
	 // initializing the z_sig
    MatrixXd z_sig = MatrixXd(lidar_data_elem_,2*n_aug_+1);
    // converting the sigma points from state space to the measurement space
    // and fill up the z_sig	
    for (int i=0; i < 2*n_aug_+1; i++)
	{
		z_sig(0,i)  		= Xsig_pred_(0,i);
		z_sig(1,i)  		= Xsig_pred_(1,i);
    }	
    // calculate the weighted of the z sigma points
    VectorXd z_sig_weighted= VectorXd(lidar_data_elem_);
	z_sig_weighted = z_sig*weights_;
	
	// calculate the diference between the weighted sigma ponits z
	// and the sigma point z
	VectorXd z_diff = VectorXd(lidar_data_elem_);
	z_diff.fill(0.0);
	
	// calculate the covarince of sigma point in z space
	MatrixXd S = MatrixXd(lidar_data_elem_,lidar_data_elem_);
	S.fill(0.0);
	
	VectorXd x_diff = VectorXd(n_x_);
	x_diff.fill(0.0);
	
	//create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, lidar_data_elem_); 
    Tc.fill(0.0);
	
	for (int i=0; i< 2*n_aug_+1; i++)
	{
		z_diff = z_sig_weighted - z_sig.col(i);
		x_diff =  x_- Xsig_pred_.col(i);
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	S = S + R_Lidar_;
	
	//Kalman gain K;
    MatrixXd K = Tc * S.inverse();

	
	//update state mean and covariance matrix
    VectorXd z_difff = z - z_sig_weighted ; 
	
	x_ = x_ + K * z_difff;
    P_ = P_ - K * S * K.transpose();
    double nis_raser = z_difff.transpose() * S.inverse() * z_difff;
	//cout << "NIS_laser" << nis_raser << endl;
	nis_lidar_list_.push_back(nis_raser);
  
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
    // obtaining the z from the measuement
	VectorXd z = VectorXd(radar_data_elem_);
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
	z(2) = meas_package.raw_measurements_(2);
    // initializing the z_sig
    MatrixXd z_sig = MatrixXd(radar_data_elem_,2*n_aug_+1);
    // converting the sigma points from state space to the measurement space
    // and fill up the z_sig	
    for (int i=0; i<2*n_aug_+1;i++)
	{
		double p_x 		= Xsig_pred_(0,i);
		double p_y 		= Xsig_pred_(1,i);
		double v 		= Xsig_pred_(2,i);
		double yaw		= Xsig_pred_(3,i);
		double yawd		= Xsig_pred_(4,i); 
	    z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);     
        z_sig(1,i) = atan2(p_y, p_x);              
        z_sig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/sqrt(p_x*p_x + p_y*p_y);	
	}	
   
	// calculate the weighted of the z sigma points
    VectorXd z_sig_weighted= VectorXd(radar_data_elem_);
	z_sig_weighted = z_sig*weights_;
	
	
	// calculate the diference between the weighted sigma ponits z
	// and the sigma point z
	VectorXd z_diff = VectorXd(radar_data_elem_);
	z_diff.fill(0.0);
	
	// calculate the covarince of sigma point in z space
	MatrixXd S = MatrixXd(radar_data_elem_,radar_data_elem_);
	S.fill(0.0);
	
	VectorXd x_diff = VectorXd(n_x_);
	x_diff.fill(0.0);
	
	//create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, radar_data_elem_); 
    Tc.fill(0.0);
	
	
	for (int i=0; i< 2*n_aug_+1; i++)
	{
		z_diff = z_sig_weighted - z_sig.col(i);
		z_diff(1) = AngelNormalization(z_diff(1));	
		x_diff =  x_- Xsig_pred_.col(i);
		x_diff(3) = AngelNormalization(x_diff(3));	
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	S = S + R_Radar_;

	
	//Kalman gain K;
    MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z_difff = z - z_sig_weighted ; 
	
    x_ = x_ + K * z_difff;
    P_ = P_ - K * S * K.transpose();
    double nis_radar = z_difff.transpose() * S.inverse() * z_difff;
	//cout << "NIS_radar" << nis_radar << endl;
	nis_radar_list_.push_back(nis_radar);
}

/*
cd /opt/carma
rosbag play cmd_speed1_no_can_speed.bag /saxton_cav/drivers/srx_controller/control/cmd_speed:=/srx_controller/control/cmd_speed
rostopic pub -r 10 /srx_controller/control/cmd_speed cav_msgs/speedAccel "speed: 0.0 maxaccel: 0.0"
*/