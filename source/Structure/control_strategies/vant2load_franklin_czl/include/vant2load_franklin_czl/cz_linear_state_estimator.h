/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file cz_linear_state_estimator.h
 * @brief This file contains the declaration of the cz_linear_state_estimator class.
 *
 * @author Brenner S. Rego
 */

#pragma once

#include <eigen3/Eigen/Eigen>

#include <vector>

#include "czonotope.h"

class cz_linear_state_estimator
{
public:
  cz_linear_state_estimator() = default;
  ~cz_linear_state_estimator() = default;

  static cz execute(int mode, cz Xhat_prev, Eigen::MatrixXd u, std::vector<double> y, std::vector<int> I, cz W, cz V);

private:
  static cz prediction(cz X, Eigen::MatrixXd u, cz W);
  static cz correction(cz Xbar, std::vector<double> y, std::vector<int> I, cz V);
};
