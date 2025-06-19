/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @brief This file contains the declaration of the function that returns the
 * a skew symmetric matrix for a given vector.
 * @author Daniel Cardoso
 */

#include <eigen3/Eigen/Eigen>
#include <cmath>

Eigen::MatrixXd SkewSymmetricMatrix(Eigen::VectorXd Vector)
{
  // Place Vet in the Skew Symmetric matrix S
  Eigen::MatrixXd SkewMatrix(3, 3);
  SkewMatrix << 0, -Vector(2), Vector(1), Vector(2), 0, -Vector(0), -Vector(1), Vector(0), 0;

  return SkewMatrix;
}
