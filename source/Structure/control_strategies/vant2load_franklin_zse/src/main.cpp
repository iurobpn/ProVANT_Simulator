/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file main.cpp
 * @brief This file contains the implementation of the control law.
 *
 * The mathematical background used to design this controller is presented in
 * the paper "Suspended Load Path Tracking Control Using a Tilt-rotor UAV Based
 * on Zonotopic State Estimation" publised in Journal of the Franklin Institute
 * from authors Brenner S. Rego and Guilherme V. Raffo, and is available at
 *  https://arxiv.org/pdf/1809.07870.pdf.
 *
 * @author Brenner S. Rego
 */

#include <control_strategies_base/icontroller.hpp>

#include <simulator_msgs/Sensor.h>

#include <eigen3/Eigen/Eigen>

#include <cmath>

#include "vant2load_franklin_zse/feedforward.h"
#include "vant2load_franklin_zse/routines.h"
#include "vant2load_franklin_zse/zonotopic_state_estimator.h"

class vant2load_Franklin_ZSE : public Icontroller
{
private:
  Eigen::VectorXd Xref;
  Eigen::VectorXd Erro;
  Eigen::VectorXd Input;  // This vector includes the applied disturbances
  Eigen::MatrixXd K;
  Eigen::VectorXd X;
  Eigen::VectorXd Xstates;
  Eigen::MatrixXd ZSE_X0c;
  Eigen::MatrixXd ZSE_X0G;
  Eigen::MatrixXd ZSE_Wc;
  Eigen::MatrixXd ZSE_WG;
  Eigen::MatrixXd ZSE_Vc;
  Eigen::MatrixXd ZSE_VG;
  int ZSE_maxorder;
  double T;
  double equil_phi;
  double equil_theta;
  double equil_psii;
  double equil_g1;
  double equil_g2;
  double equil_aR;
  double equil_aL;
  double equil_fR;
  double equil_fL;
  double equil_tauaR;
  double equil_tauaL;
  Eigen::MatrixXd equil_output;
  double PI;

public:
  vant2load_Franklin_ZSE()
    : Xref(24)
    , Erro(24)
    , Input(7)
    , K(4, 24)
    , X(24 + 23 + 46 + 1 + 1 + 1 + 3 + 3)
    , Xstates(24)
    , ZSE_X0c(23, 1)
    , ZSE_X0G(23, 23)
    , ZSE_Wc(23, 1)
    , ZSE_WG(23, 23)
    , ZSE_Vc(16, 1)
    , ZSE_VG(16, 16)
    , equil_output(16, 1)
  {
    // Controller sampling time
    T = 0.012;

    // Equilibrium values
    equil_phi = 0;
    equil_theta = 0;
    equil_psii = 0;
    equil_g1 = 1.316967145565129e-04;
    equil_g2 = 0.013960153151247;
    equil_aR = 0.014005280124626;
    equil_aL = 0.013809089793911;
    equil_fR = 11.732256727680697;
    equil_fL = 11.767602459917352;
    equil_tauaR = 4.138868162580376e-07;
    equil_tauaL = 1.012099978192413e-05;

    equil_output << -0.001661204266230097, 1.567038189486787e-05, 0.6189884034468275, -0.000131709548513303,
        -0.0139601530301764, 1.838625741689666e-06, 0, 0, 0, -0.006979849797580309, 6.584835708791055e-05,
        -0.4999512749866701, 0.01400528012462612, 0.01380908979391114, 0, 0;

    // Mathematical constant PI
    PI = 3.14159265358979323846;
  }

  ~vant2load_Franklin_ZSE()
  {
  }

  void config()
  {
    // [NEW] H2Hinf Trajectory 2 (reduced domain of attraction): Gain Matrix 8 [TRY THIS ONE]
    K << -0.9132946313073804, 5.608219057174488, 6.049704650196739, -23.20025951654861, -1.630546515546394,
        0.02148077682426153, 14.65912470456652, 0.6728633659104601, -0.01235464436009421, -1.182733819870061,
        -0.6795341654247624, 5.631597220118333, 4.260767271616122, -5.301088138961158, -0.3011658470382209,
        0.03028996800792656, 1.651747291855211, 0.05700756111732627, -0.008342173470425538, -0.02199791243939212,
        -0.474869541799516, 2.398147698508092, 3.206342458002378, 0.01756267237423986, -0.8371352936136646,
        -6.451175152682915, 5.714997803891847, 23.74476540619154, -1.359244288115006, -0.03501065071649617,
        -14.7614012535182, 0.5296927819666141, -1.177068417141205, 0.06488221914310917, -0.5996191351738775,
        -6.052186540020592, 4.333474628422141, 5.368238675233806, -0.252722863927784, -0.03743399824551547,
        -1.660067835952305, 0.0424408621650129, -0.02188313316014417, -0.00462821048595196, -0.4158963883285703,
        -2.944596222015817, 2.795327009132437, -0.01736331816698299, 0.01331388589603127, 0.01654085920647495,
        0.0003927350302298688, -0.0617037094224527, 0.05486112406384938, 0.004539640638941521, 0.03704954900985658,
        -0.0780765301973939, -0.2945449582188707, 0.06945469059894138, 0.01269067517417017, 0.01581301280421729,
        -0.0005887749746805169, -0.01384427260989705, 0.01345076981568556, 0.002619426073617202, 0.004222211097536802,
        -0.005196568005687351, -0.0007507059886165572, -0.0001756837151079768, 0.005919234919065153,
        0.007306907668689294, 0.0007442649416933643, 0.002440815396601709, 0.01314508293317239, -0.01624201547570362,
        -0.0004436165517214821, 0.06118101463472548, 0.05573754272493357, -0.004697199104652598, -0.03668599270697293,
        -0.07932540983009606, 0.06946780394207275, -0.2940573375792874, 0.01276919356825688, -0.0156726380548064,
        -0.0003891507859917382, 0.01371836392024517, 0.01374154472879415, -0.002671458920751988, -0.004166211578938072,
        -0.005285545502828642, -0.0001738078621249789, -0.0007143633361166265, 0.005788438251821402,
        -0.007158529385822105, -0.0003383310614920475, -0.002509037396835588;

    // ZSE Adjust Parameters
    ZSE_X0c << 0, 0, 1, equil_phi, equil_theta, equil_psii, equil_g1, equil_g2, equil_aR, equil_aL, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0;  // x hat inicial
    ZSE_X0G << Eigen::MatrixXd::Zero(23, 23);
    ZSE_X0G(0, 0) = 0.5;  // NEW adjust for final version (30 deg)
    ZSE_X0G(1, 1) = 0.5;
    ZSE_X0G(2, 2) = 0.5;
    ZSE_X0G(3, 3) = PI / 4.0;
    ZSE_X0G(4, 4) = PI / 4.0;
    ZSE_X0G(5, 5) = PI / 4.0;
    ZSE_X0G(6, 6) = PI / 4.0;
    ZSE_X0G(7, 7) = PI / 4.0;
    ZSE_X0G(8, 8) = PI / 4.0;
    ZSE_X0G(9, 9) = PI / 4.0;
    ZSE_X0G(10, 10) = 1.0;
    ZSE_X0G(11, 11) = 1.0;
    ZSE_X0G(12, 12) = 1.0;
    ZSE_X0G(13, 13) = 1.0;
    ZSE_X0G(14, 14) = 1.0;
    ZSE_X0G(15, 15) = 1.0;
    ZSE_X0G(16, 16) = 1.0;
    ZSE_X0G(17, 17) = 1.0;
    ZSE_X0G(18, 18) = 1.0;
    ZSE_X0G(19, 19) = 1.0;
    ZSE_X0G(20, 20) = 1.0;
    ZSE_X0G(21, 21) = 1.0;
    ZSE_X0G(22, 22) = 1.0;

    ZSE_Wc << Eigen::MatrixXd::Zero(23, 1);
    ZSE_WG << Eigen::MatrixXd::Zero(23, 23);
    ZSE_WG(0, 0) = 1.0E-04;
    ZSE_WG(1, 1) = 1.0E-04;
    ZSE_WG(2, 2) = 1.0E-04;
    ZSE_WG(3, 3) = 1.0E-04;
    ZSE_WG(4, 4) = 1.0E-04;
    ZSE_WG(5, 5) = 1.0E-04;
    ZSE_WG(6, 6) = 1.0E-04;
    ZSE_WG(7, 7) = 1.0E-04;
    ZSE_WG(8, 8) = 1.5E-04;
    ZSE_WG(9, 9) = 1.5E-04;
    ZSE_WG(10, 10) = 1.0E-04;
    ZSE_WG(11, 11) = 1.0E-04;
    ZSE_WG(12, 12) = 1.0E-04;
    ZSE_WG(13, 13) = 0.01;
    ZSE_WG(14, 14) = 0.01;
    ZSE_WG(15, 15) = 0.01;
    ZSE_WG(16, 16) = 0.05;
    ZSE_WG(17, 17) = 0.05;
    ZSE_WG(18, 18) = 1.0E-04;
    ZSE_WG(19, 19) = 1.0E-04;
    ZSE_WG(20, 20) = 0.01;
    ZSE_WG(21, 21) = 0.01;
    ZSE_WG(22, 22) = 0.01;

    ZSE_Vc << equil_output;
    ZSE_VG << Eigen::MatrixXd::Zero(16, 16);
    ZSE_VG(0, 0) = 0.18;
    ZSE_VG(1, 1) = 0.18;
    ZSE_VG(2, 2) = 0.612;
    ZSE_VG(3, 3) = 3.1416E-03;
    ZSE_VG(4, 4) = 3.1416E-03;
    ZSE_VG(5, 5) = 0.03;
    ZSE_VG(6, 6) = 19.872E-03;
    ZSE_VG(7, 7) = 19.872E-03;
    ZSE_VG(8, 8) = 0.24;
    ZSE_VG(9, 9) = 0.006;
    ZSE_VG(10, 10) = 0.006;
    ZSE_VG(11, 11) = 0.06;
    ZSE_VG(12, 12) = 6.8067E-03;
    ZSE_VG(13, 13) = 6.8067E-03;
    ZSE_VG(14, 14) = 0.6093;
    ZSE_VG(15, 15) = 0.6093;

    ZSE_maxorder = 23 * 75;
  }

  std::vector<double> execute(simulator_msgs::SensorArray arraymsg)
  {
    // Instant time k counter
    static double count = 0;

    // Integrator variables
    static double xint, x_ant = 0;
    static double yint, y_ant = 0;
    static double zint, z_ant = 0;
    static double yawint, yaw_ant = 0;

    // Equilibrium states and inputs (for state estimation)
    Eigen::MatrixXd Xeq(23, 1);
    Eigen::MatrixXd Ueq(4, 1);
    Xeq << 0, 0, 0, equil_phi, equil_theta, equil_psii, equil_g1, equil_g2, equil_aR, equil_aL, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0;
    Ueq << equil_fR, equil_fL, equil_tauaR, equil_tauaL;

    static Eigen::MatrixXd Inputprev(4, 1);
    if (count == 0)
    {
      Inputprev << Ueq;  // u anterior
    }

    // ZSE variables
    static Eigen::MatrixXd ZSE_Xc(23, 1);
    static Eigen::MatrixXd ZSE_Xcprev(23, 1);
    static Eigen::MatrixXd ZSE_XG(23, 23);
    static Eigen::MatrixXd ZSE_XGprev(23, 23);

    // Getting sensor data
    simulator_msgs::Sensor msgstates;
    simulator_msgs::Sensor mainbody_sensordata;
    simulator_msgs::Sensor rightpropeller_sensordata;
    simulator_msgs::Sensor leftpropeller_sensordata;
    simulator_msgs::Sensor load_sensordata;
    simulator_msgs::Sensor rightservo_sensordata;
    simulator_msgs::Sensor leftservo_sensordata;
    simulator_msgs::Sensor rodx_sensordata;
    simulator_msgs::Sensor rody_sensordata;

    msgstates = arraymsg.values.at(0);
    mainbody_sensordata = arraymsg.values.at(1);
    rightpropeller_sensordata = arraymsg.values.at(2);
    leftpropeller_sensordata = arraymsg.values.at(3);
    load_sensordata = arraymsg.values.at(4);
    rightservo_sensordata = arraymsg.values.at(5);
    leftservo_sensordata = arraymsg.values.at(6);
    rodx_sensordata = arraymsg.values.at(7);
    rody_sensordata = arraymsg.values.at(8);

    // Reference

    // Generate trajectory
    Eigen::MatrixXd csiref(3, 1);
    Eigen::MatrixXd csirefdot(3, 1);
    Eigen::MatrixXd csirefddot(3, 1);
    int which_trajectory = 2;  // Second trajectory from Rego's Master thesis
    routines::generate_Trajectory(count, T, which_trajectory, csiref, csirefdot, csirefddot);

    Xref << csiref, equil_phi, equil_theta, equil_psii, equil_g1, equil_g2, equil_aR, equil_aL, csirefdot, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0;

    // Disturbance forces
    Eigen::MatrixXd dtrb(3, 1);
    static Eigen::MatrixXd dtrb_prev(3, 1);
    int which_disturbance = 0;  // No disturbance forces
    dtrb = routines::generate_Disturbances(count, T, which_disturbance);

    // Wind disturbances
    Eigen::MatrixXd wind(3, 1);
    static Eigen::MatrixXd wind_prev(3, 1);
    // int which_wind = 0; // No wind disturbances
    int which_wind = 1;  // Wind disturbances based on Rego's Master thesis, second trajectory
    wind = routines::generate_Wind(count, T, which_wind);

    // Disturbance forces filtering
    if (count == 0)
    {
      dtrb_prev << 0, 0, 0;
    }
    dtrb(0, 0) = routines::first_order_filter(dtrb_prev(0, 0), dtrb(0, 0));
    dtrb(1, 0) = routines::first_order_filter(dtrb_prev(1, 0), dtrb(1, 0));
    dtrb(2, 0) = routines::first_order_filter(dtrb_prev(2, 0), dtrb(2, 0));
    dtrb_prev << dtrb;

    // Wind filtering
    if (count == 0)
    {
      wind_prev << 0, 0, 0;
    }
    wind(0, 0) = routines::first_order_filter(wind_prev(0, 0), wind(0, 0));
    wind(1, 0) = routines::first_order_filter(wind_prev(1, 0), wind(1, 0));
    wind(2, 0) = routines::first_order_filter(wind_prev(2, 0), wind(2, 0));
    wind_prev << wind;

    // Aerodynamic disturbances
    Eigen::MatrixXd aerodtrb(3, 1);
    aerodtrb = routines::generate_AerodynDisturbances(
        msgstates.values.at(10), msgstates.values.at(11), msgstates.values.at(12),  // xdot, ydot, zdot
        msgstates.values.at(3), msgstates.values.at(4), msgstates.values.at(5),     // phi, theta, psi
        wind);                                                                      // Wind disturbances

    // Generate current measurement
    std::vector<double> y;
    std::vector<int> I;
    routines::generate_Measurement_truncated(count, T, mainbody_sensordata, rightpropeller_sensordata,
                                             leftpropeller_sensordata, load_sensordata, rightservo_sensordata,
                                             leftservo_sensordata, rodx_sensordata, rody_sensordata, y, I);

    // Zonotopic state estimator estimation
    if (count == 0)
    {
      ZSE_Xcprev << ZSE_X0c;
      ZSE_XGprev << ZSE_X0G;
      ZSE_XG.resize(23, ZSE_XGprev.cols() + 16);
      zonotopic_state_estimator::execute(1, ZSE_Xcprev - Xeq, ZSE_XGprev, ZSE_Wc, ZSE_WG, ZSE_Vc, ZSE_VG,
                                         Inputprev - Ueq, y, I, ZSE_maxorder, ZSE_Xc, ZSE_XG);  // Initial correction
    }
    else
    {
      int next_order = ZSE_XGprev.cols() + 23 + I.size();
      if (next_order > ZSE_maxorder)
      {
        ZSE_XG.resize(23, ZSE_maxorder);
      }
      else
      {
        ZSE_XG.resize(23, next_order);
      }

      zonotopic_state_estimator::execute(0, ZSE_Xcprev - Xeq, ZSE_XGprev, ZSE_Wc, ZSE_WG, ZSE_Vc, ZSE_VG,
                                         Inputprev - Ueq, y, I, ZSE_maxorder, ZSE_Xc, ZSE_XG);  // Estimation
    }
    ZSE_Xc = ZSE_Xc + Xeq;

    // Zonotopic state estimator confidence limits
    Eigen::MatrixXd ZSE_conflimits(46, 1);
    ZSE_conflimits = routines::get_ZSEconfidencelimits(ZSE_XG);

    //		// Zonotopic state estimator Frobenius norm
    double ZSE_frobnorm = routines::get_ZSEfrobnorm(ZSE_XG);

    // Convertendo velocidade angular
    std::vector<double> etadot = routines::wIIB2EtaDot(msgstates.values.at(13),  // wIIL_x
                                                       msgstates.values.at(14),  // wIIL_y
                                                       msgstates.values.at(15),  // wIIL_z
                                                       msgstates.values.at(3),   // phi
                                                       msgstates.values.at(4),   // theta
                                                       msgstates.values.at(5));  // psii

    // Integrador Trapezoidal (usando estados estimados pelo ZSE)
    double x_atual = ZSE_Xc(0) - Xref(0);
    xint = xint + (T / 2) * (x_atual + x_ant);
    x_ant = x_atual;
    double y_atual = ZSE_Xc(1) - Xref(1);
    yint = yint + (T / 2) * (y_atual + y_ant);
    y_ant = y_atual;
    double z_atual = ZSE_Xc(2) - Xref(2);
    zint = zint + (T / 2) * (z_atual + z_ant);
    z_ant = z_atual;
    double yaw_atual = ZSE_Xc(5) - Xref(5);
    yawint = yawint + (T / 2) * (yaw_atual + yaw_ant);
    yaw_ant = yaw_atual;

    // State vector (augmented)
    Xstates << msgstates.values.at(0),  // x
        msgstates.values.at(1),         // y
        msgstates.values.at(2),         // z
        msgstates.values.at(3),         // roll
        msgstates.values.at(4),         // pitch
        msgstates.values.at(5),         // yaw
        msgstates.values.at(8),         // g1 x
        msgstates.values.at(9),         // g2 y
        msgstates.values.at(6),         // aR
        msgstates.values.at(7),         // aL
        msgstates.values.at(10),        // vx
        msgstates.values.at(11),        // vy
        msgstates.values.at(12),        // vz
        etadot.at(0),                   // droll
        etadot.at(1),                   // pitch
        etadot.at(2),                   // yaw
        msgstates.values.at(18),        // g1dot
        msgstates.values.at(19),        // g2dot
        msgstates.values.at(16),        // aRdot
        msgstates.values.at(17),        // aLdot
        xint, yint, zint, yawint;

    // State vector ZSE (augmented)
    Eigen::MatrixXd Xstates_ZSE(24, 1);
    Xstates_ZSE << ZSE_Xc.block(0, 0, 20, 1),  // xhat
        xint, yint, zint, yawint;

    // Tracking error
    Erro = Xstates - Xref;

    // Tracking error (using estimation from ZSE)
    Eigen::MatrixXd Erro_ZSE(24, 1);
    Erro_ZSE = Xstates_ZSE - Xref;

    Input << -K * Erro_ZSE, dtrb + aerodtrb;

    // Reference feedfoward
    Eigen::MatrixXd qref(10, 1);
    qref << csiref, equil_phi, equil_theta, equil_psii, equil_g1, equil_g2, equil_aR, equil_aL;
    Eigen::MatrixXd qrefdot(10, 1);
    qrefdot << csirefdot, 0, 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd qrefddot(10, 1);
    qrefddot << csirefddot, 0, 0, 0, 0, 0, 0, 0;
    Eigen::MatrixXd uref(7, 1);
    uref << feedforward::compute(qref, qrefdot, qrefddot), 0, 0, 0;
    Input = Input + uref;

    count++;

    // Storing values for the next iteration
    ZSE_Xcprev << ZSE_Xc;
    ZSE_XGprev.resize(23, ZSE_XG.cols());
    ZSE_XGprev << ZSE_XG;
    Inputprev << Input.block(0, 0, 4, 1);

    // States, ZSE data, size of measured vector and disturbances
    X << Xstates, ZSE_Xc, ZSE_conflimits, ZSE_frobnorm, y.size(), I.size(), dtrb, aerodtrb;

    // Stop criterion
    if (count >= 5001)
    {
      X << 1;  // Force error in Eigen, kills the controller process
    }

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
  return new vant2load_Franklin_ZSE;
}
void destroy(Icontroller* p)
{
  delete p;
}
}
