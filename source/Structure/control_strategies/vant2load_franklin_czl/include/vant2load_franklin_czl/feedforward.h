/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file feedforward.h
 * @brief This file contains the declaration of the feedforward class.
 *
 * @author Brenner S. Rego
 */

#pragma once

#include <eigen3/Eigen/Eigen>

class feedforward
{
public:
  feedforward() = default;
  ~feedforward() = default;
  static Eigen::MatrixXd compute(Eigen::MatrixXd qref, Eigen::MatrixXd qrefdot, Eigen::MatrixXd qrefddot);
};
