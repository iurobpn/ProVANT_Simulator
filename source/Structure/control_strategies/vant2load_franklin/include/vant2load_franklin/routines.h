/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file routines.h
 * @brief This file contains the declaration of the routines class.
 *
 * @author Brenner Santana Rego
 */

#pragma once

#include <vector>

#include <eigen3/Eigen/Eigen>

#include <simulator_msgs/Sensor.h>

class routines
{
public:
  routines() = default;
  ~routines() = default;

  static void generate_Trajectory(int k, double Ts, int which_trajectory, Eigen::MatrixXd& csiref,
                                  Eigen::MatrixXd& csirefdot, Eigen::MatrixXd& csirefddot);
  static Eigen::MatrixXd generate_Disturbances(int k, double Ts, int which_disturbance);
  static void generate_Measurement_clean(int k, double Ts, simulator_msgs::Sensor mainbody,
                                         simulator_msgs::Sensor rightpropeller, simulator_msgs::Sensor leftpropeller,
                                         simulator_msgs::Sensor load, simulator_msgs::Sensor rightservo,
                                         simulator_msgs::Sensor leftservo, simulator_msgs::Sensor rodx,
                                         simulator_msgs::Sensor rody, std::vector<double>& y, std::vector<int>& I);
  static double first_order_filter(double in_prev, double in);
  static void generate_Measurement(int k, double Ts, simulator_msgs::Sensor mainbody,
                                   simulator_msgs::Sensor rightpropeller, simulator_msgs::Sensor leftpropeller,
                                   simulator_msgs::Sensor load, simulator_msgs::Sensor rightservo,
                                   simulator_msgs::Sensor leftservo, simulator_msgs::Sensor rodx,
                                   simulator_msgs::Sensor rody, std::vector<double>& y, std::vector<int>& I);
  static std::vector<double> wIIB2EtaDot(double in_a, double in_b, double in_c, double phi, double theta, double psii);
  static std::vector<double> pqr2EtaDot(double in_a, double in_b, double in_c, double phi, double theta, double psii);
  static std::vector<double> wIIL2pqr(double in_a, double in_b, double in_c, double phi, double theta, double psii);
  static Eigen::MatrixXd rotx(double angle);
  static Eigen::MatrixXd roty(double angle);
  static Eigen::MatrixXd rotz(double angle);
  static Eigen::MatrixXd get_LKFconfidencelimits(Eigen::MatrixXd Pxx);
  static Eigen::MatrixXd get_ZSEconfidencelimits(Eigen::MatrixXd XG);
  static double get_ZSEfrobnorm(Eigen::MatrixXd XG);
};
