/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file zonotopic_state_estimator.h
 * @brief This file contains the declaration of the zonotopic_state_estimator class.
 *
 * @author Brenner S. Rego
 */

#include <eigen3/Eigen/Eigen>

#include <vector>

class zonotopic_state_estimator
{
public:
  zonotopic_state_estimator() = default;
  ~zonotopic_state_estimator() = default;

  static void execute(int mode, Eigen::MatrixXd Xhatprev_c, Eigen::MatrixXd Xhatprev_G, Eigen::MatrixXd W_c,
                      Eigen::MatrixXd W_G, Eigen::MatrixXd V_c, Eigen::MatrixXd V_G, Eigen::MatrixXd deltau,
                      std::vector<double> y, std::vector<int> I, int max_order, Eigen::MatrixXd& Xhat_c,
                      Eigen::MatrixXd& Xhat_G);

private:
  static void prediction(Eigen::MatrixXd X_c, Eigen::MatrixXd X_G, Eigen::MatrixXd W_c, Eigen::MatrixXd W_G,
                         Eigen::MatrixXd deltau, Eigen::MatrixXd Anu, Eigen::MatrixXd Bnu, Eigen::MatrixXd& Xbar_c,
                         Eigen::MatrixXd& Xbar_G);
  static void strip_computation(Eigen::MatrixXd Hnu, Eigen::MatrixXd V_c, Eigen::MatrixXd V_G, int i,
                                Eigen::MatrixXd& rho, double& s, double& sigma);
  static void intersection(Eigen::MatrixXd Xtilde_c, Eigen::MatrixXd Xtilde_G, Eigen::MatrixXd rho, double s,
                           double sigma, double y_i, int i, Eigen::MatrixXd& Xtildenew_c, Eigen::MatrixXd& Xtildenew_G);
  static void order_reduction(Eigen::MatrixXd X_c, Eigen::MatrixXd X_G, int r_max, Eigen::MatrixXd& Xred_c,
                              Eigen::MatrixXd& Xred_G);
};
