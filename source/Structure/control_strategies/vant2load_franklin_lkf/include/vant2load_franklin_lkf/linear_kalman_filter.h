/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file linear_kalman_filter.h
 * @brief This file contains the declaration of the linear_kalman_filter class.
 *
 * @author Brenner S. Rego
 */

#pragma once

#include <eigen3/Eigen/Eigen>

class linear_kalman_filter
{
public:
  linear_kalman_filter() = default;
  ~linear_kalman_filter() = default;

  static void execute(Eigen::MatrixXd xhat_plus_prev, Eigen::MatrixXd u_prev, std::vector<double> y, std::vector<int> I,
                      Eigen::MatrixXd Pxx_plus_prev, Eigen::MatrixXd Q, Eigen::MatrixXd R, Eigen::MatrixXd& xhat_plus,
                      Eigen::MatrixXd& Pxx_plus);
};
