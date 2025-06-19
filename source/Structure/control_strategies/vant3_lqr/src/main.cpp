/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file This file contains the implementation of the vant3_lqr class for
 * UAV 3.0, and is used to test the UAV.
 *
 * @author Daniel Neri Cardoso
 */

#include <control_strategies_base/icontroller.hpp>

#include <eigen3/Eigen/Eigen>
#include "simulator_msgs/Sensor.h"

class vant3_lqr : public Icontroller
{
private:
  Eigen::VectorXd Xref;
  Eigen::VectorXd Erro;
  Eigen::VectorXd Input;
  Eigen::MatrixXd K;
  Eigen::VectorXd X;
  Eigen::VectorXd Ur;
  double T;

public:
  vant3_lqr() : Xref(16), K(6, 16), X(16), Erro(16), Input(6), Ur(6)
  {
    double T = 0.012;
  }
  ~vant3_lqr()
  {
  }

  void config()
  {
    K << -0.0027365, 7.3015, 7.3029, -26.668, -0.019972, 0.11216, -0.62175, 0.57331, -0.0038204, 7.9694, 6.0275,
        -4.4679, -0.0056377, 0.066199, 0.023509, -0.027802, -0.0030025, -7.3015, 7.303, 26.668, -0.018014, -0.11216,
        0.57888, -0.61622, -0.003877, -7.9693, 6.0276, 4.4679, -0.0044893, -0.066194, -0.027294, 0.02401, 0.035355,
        -0.00071141, 1.3878e-05, 0.0056206, 0.24147, 0.027004, 0.40342, 0.06954, 0.047936, -0.0010205, 1.1657e-05,
        0.0018148, 0.066616, 0.02518, 0.042581, 0.0014621, 0.035355, 0.00071012, 1.3906e-05, -0.0055989, 0.24147,
        -0.027004, 0.069545, 0.40342, 0.047936, 0.0010176, 1.168e-05, -0.0018076, 0.066615, -0.02518, 0.0014623,
        0.042581, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Xref << 0, 0, 1, 0.000144, -0.196, 0, 0.196, 0.196, 0, 0, 0, 0, 0, 0, 0, 0;

    Ur << 8.27, 8.24, 7.71e-05, 6.67e-05, 0, 0;
  }

  std::vector<double> execute(simulator_msgs::SensorArray arraymsg)
  {
    simulator_msgs::Sensor msg;
    msg = arraymsg.values.at(0);

    X << msg.values.at(0), msg.values.at(1), msg.values.at(2), msg.values.at(3), msg.values.at(4), msg.values.at(5),
        msg.values.at(6), msg.values.at(7), msg.values.at(8), msg.values.at(9), msg.values.at(10), msg.values.at(11),
        msg.values.at(12), msg.values.at(13), msg.values.at(14), msg.values.at(15);

    Erro = X - Xref;
    Input = -K * Erro;
    Input = Input + Ur;

    std::vector<double> out(Input.data(), Input.data() + Input.rows() * Input.cols());
    return out;
  }

  std::vector<double> Reference()
  {
    std::vector<double> out(Xref.data(), Xref.data() + Xref.rows() * Xref.cols());
    return out;
  }

  std::vector<double> Error()
  {
    std::vector<double> out(Erro.data(), Erro.data() + Erro.rows() * Erro.cols());
    return out;
  }

  std::vector<double> State()
  {
    std::vector<double> out(X.data(), X.data() + X.rows() * X.cols());
    return out;
  }
};

extern "C" {
Icontroller* create(void)
{
  return new vant3_lqr;
}
void destroy(Icontroller* p)
{
  delete p;
}
}
