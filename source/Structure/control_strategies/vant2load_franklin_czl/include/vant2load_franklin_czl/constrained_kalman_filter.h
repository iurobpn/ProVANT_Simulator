/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file constrained_kalman_filter.h
 * @brief This file contains the declaration of the constrained_kalman_filter class.
 *
 * @author Brenner S. Rego
 */

#pragma once

#include <iostream>

#include <eigen3/Eigen/Eigen>

#include "math.h"
#include "czonotope.h"

class constrained_kalman_filter
{
public:
  constrained_kalman_filter() = default;
  ~constrained_kalman_filter() = default;
  static void execute(int mode, Eigen::MatrixXd xhat_plus_prev, Eigen::MatrixXd P_plus_prev, Eigen::MatrixXd u_prev,
                      std::vector<double> y, std::vector<int> I, Eigen::MatrixXd Q, Eigen::MatrixXd R, cz Xhat,
                      Eigen::MatrixXd& xhat_plus, Eigen::MatrixXd& P_plus);

private:
  static void prediction(Eigen::MatrixXd xhat_plus_prev, Eigen::MatrixXd P_plus_prev, Eigen::MatrixXd u_prev,
                         Eigen::MatrixXd Q, Eigen::MatrixXd& xhat_minus, Eigen::MatrixXd& P_minus);
  static void correction(Eigen::MatrixXd xhat_minus, Eigen::MatrixXd P_minus, std::vector<double> y, std::vector<int> I,
                         Eigen::MatrixXd R, cz Xhat, Eigen::MatrixXd& xhat_plus, Eigen::MatrixXd& P_plus);
};
