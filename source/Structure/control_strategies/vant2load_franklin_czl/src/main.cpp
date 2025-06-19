/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file main.cpp
 * @brief This file contains the main implementation of the control law.
 *
 * The mathematical background used to design this controller is presented in
 * the paper "Suspended Load Path Tracking Control Using a Tilt-rotor UAV
 * Based on Zonotopic State Estimation" publised in Journal of the Franklin
 * Institute from authors Brenner S. Rego and Guilherme V. Raffo, and is
 * available at https://arxiv.org/pdf/1809.07870.pdf
 *
 * @author Brenner S. Rego
 */

#include <control_strategies_base/icontroller.hpp>

#include <eigen3/Eigen/Eigen>

#include <simulator_msgs/Sensor.h>

#include "vant2load_franklin_czl/czonotope.h"
#include "vant2load_franklin_czl/constrained_kalman_filter.h"
#include "vant2load_franklin_czl/cz_linear_state_estimator.h"
#include "vant2load_franklin_czl/feedforward.h"
#include "vant2load_franklin_czl/routines.h"

#include <cmath>

class vant2load_Franklin_CZL : public Icontroller
{
private:
  Eigen::VectorXd Xref;
  Eigen::VectorXd Erro;
  Eigen::VectorXd Input;  // This vector includes the applied disturbances
  Eigen::MatrixXd K;
  Eigen::VectorXd X;
  Eigen::VectorXd Xstates;
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
  cz Xhat;
  cz Xhatcand;
  bool Xhatisempty;
  Eigen::MatrixXd Xhathull;
  Eigen::MatrixXd MidHull;
  Eigen::MatrixXd xA0;
  Eigen::MatrixXd xb0;
  Eigen::MatrixXd wA;
  Eigen::MatrixXd wb;
  Eigen::MatrixXd vA;
  Eigen::MatrixXd vb;
  cz X0;
  cz W;
  cz V;
  int ng_max;
  int nc_max;
  bool complexity_reduction;
  Eigen::MatrixXd ZSE_X0c;
  Eigen::MatrixXd ZSE_X0G;
  Eigen::MatrixXd ZSE_Wc;
  Eigen::MatrixXd ZSE_WG;
  Eigen::MatrixXd ZSE_Vc;
  Eigen::MatrixXd ZSE_VG;
  Eigen::MatrixXd CKF_xhat_plus;
  Eigen::MatrixXd CKF_P_plus;
  Eigen::MatrixXd CKF_xhat_plus_prev;
  Eigen::MatrixXd CKF_P_plus_prev;
  Eigen::MatrixXd CKF_xhat0;
  Eigen::MatrixXd CKF_Pxx0;
  Eigen::MatrixXd CKF_Q;
  Eigen::MatrixXd CKF_R;

public:
  vant2load_Franklin_CZL()
    : Xref(24)
    , Erro(24)
    , Input(7)
    , K(4, 24)
    , X(24 + 23 + 46 + 23 + 1 + 1 + 1 + 1 + 1 + 1 + 23)
    , Xstates(24)
    ,

    ZSE_X0c(23, 1)
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

  ~vant2load_Franklin_CZL()
  {
  }

  void config()
  {
    // ZSE Adjust Parameters
    // inicial
    ZSE_X0c << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  // deltax hat inicial
    ZSE_X0G << Eigen::MatrixXd::Zero(23, 23);
    ZSE_X0G(0, 0) = 0.5;
    ZSE_X0G(1, 1) = 0.5;
    ZSE_X0G(2, 2) = 0.5;
    ZSE_X0G(3, 3) = 0.2 * std::abs(equil_phi) + 0.02;
    ZSE_X0G(4, 4) = 0.2 * std::abs(equil_theta) + 0.02;
    ZSE_X0G(5, 5) = PI / 180;
    ZSE_X0G(6, 6) = 1.5 * std::abs(equil_g1);
    ZSE_X0G(7, 7) = 1.5 * std::abs(equil_g2);
    ZSE_X0G(8, 8) = 1.5 * std::abs(equil_aR);
    ZSE_X0G(9, 9) = 1.5 * std::abs(equil_aL);
    ZSE_X0G(10, 10) = 0.02;
    ZSE_X0G(11, 11) = 0.02;
    ZSE_X0G(12, 12) = 0.02;
    ZSE_X0G(13, 13) = 0.02;
    ZSE_X0G(14, 14) = 0.02;
    ZSE_X0G(15, 15) = 0.02;
    ZSE_X0G(16, 16) = 0.02;
    ZSE_X0G(17, 17) = 0.02;
    ZSE_X0G(18, 18) = 0.02;
    ZSE_X0G(19, 19) = 0.02;
    ZSE_X0G(20, 20) = 0.02;
    ZSE_X0G(21, 21) = 0.02;
    ZSE_X0G(22, 22) = 0.02;

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
    ZSE_WG(8, 8) = 1.0E-04;
    ZSE_WG(9, 9) = 1.0E-04;
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

    X0 = czonotope::create(ZSE_X0c, ZSE_X0G, xA0, xb0);
    W = czonotope::create(ZSE_Wc, ZSE_WG, wA, wb);
    V = czonotope::create(ZSE_Vc, ZSE_VG, vA, vb);

    ng_max = 50;
    nc_max = 10;
    complexity_reduction = true;

    // Constrained Kalman filter parameters
    CKF_xhat0 = Eigen::MatrixXd::Zero(23, 1);
    CKF_xhat0 << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    CKF_Pxx0 = 1.0 * Eigen::MatrixXd::Identity(23, 23);

    CKF_Q = Eigen::MatrixXd::Zero(23, 23);
    CKF_Q(0, 0) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(1, 1) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(2, 2) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(3, 3) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(4, 4) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(5, 5) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(6, 6) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(7, 7) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(8, 8) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(9, 9) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(10, 10) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(11, 11) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(12, 12) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(13, 13) = pow(0.01 / 3.0, 2);
    CKF_Q(14, 14) = pow(0.01 / 3.0, 2);
    CKF_Q(15, 15) = pow(0.01 / 3.0, 2);
    CKF_Q(16, 16) = pow(0.05 / 3.0, 2);
    CKF_Q(17, 17) = pow(0.05 / 3.0, 2);
    CKF_Q(18, 18) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(19, 19) = pow(1.0E-04 / 3.0, 2);
    CKF_Q(20, 20) = pow(0.01 / 3.0, 2);
    CKF_Q(21, 21) = pow(0.01 / 3.0, 2);
    CKF_Q(22, 22) = pow(0.01 / 3.0, 2);

    CKF_R = Eigen::MatrixXd::Zero(16, 16);
    CKF_R(0, 0) = pow(0.18 / 3.0, 2);
    CKF_R(1, 1) = pow(0.18 / 3.0, 2);
    CKF_R(2, 2) = pow(0.612 / 3.0, 2);
    CKF_R(3, 3) = pow(3.1416E-03 / 3.0, 2);
    CKF_R(4, 4) = pow(3.1416E-03 / 3.0, 2);
    CKF_R(5, 5) = pow(0.03 / 3.0, 2);
    CKF_R(6, 6) = pow(19.872E-03 / 3.0, 2);
    CKF_R(7, 7) = pow(19.872E-03 / 3.0, 2);
    CKF_R(8, 8) = pow(0.24 / 3.0, 2);
    CKF_R(9, 9) = pow(0.006 / 3.0, 2);
    CKF_R(10, 10) = pow(0.006 / 3.0, 2);
    CKF_R(11, 11) = pow(0.06 / 3.0, 2);
    CKF_R(12, 12) = pow(6.8067E-03 / 3.0, 2);
    CKF_R(13, 13) = pow(6.8067E-03 / 3.0, 2);
    CKF_R(14, 14) = pow(0.6093 / 3.0, 2);
    CKF_R(15, 15) = pow(0.6093 / 3.0, 2);
  }

  std::vector<double> execute(simulator_msgs::SensorArray arraymsg)
  {
    // printf("\nCheguei aqui\n");
    // std::cout << "Cheguei aqui" << std::endl;

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

    // Disturbances
    Eigen::MatrixXd dtrb(3, 1);
    static Eigen::MatrixXd dtrb_prev(3, 1);
    // int which_disturbance = 0; // No disturbances
    int which_disturbance = 1;  // Disturbances from Rego's Master thesis, second trajectory
    dtrb = routines::generate_Disturbances(count, T, which_disturbance);

    // Disturbance filtering
    if (count == 0)
    {
      dtrb_prev << 0, 0, 0;
    }
    dtrb(0, 0) = routines::first_order_filter(dtrb_prev(0, 0), dtrb(0, 0));
    dtrb(1, 0) = routines::first_order_filter(dtrb_prev(1, 0), dtrb(1, 0));
    dtrb(2, 0) = routines::first_order_filter(dtrb_prev(2, 0), dtrb(2, 0));
    dtrb_prev << dtrb;

    // Generate current measurement
    std::vector<double> y;
    std::vector<int> I;
    routines::generate_Measurement(count, T, mainbody_sensordata, rightpropeller_sensordata, leftpropeller_sensordata,
                                   load_sensordata, rightservo_sensordata, leftservo_sensordata, rodx_sensordata,
                                   rody_sensordata, y, I);

    // CZ linear state estimator
    if (count == 0)  // Initial correction
    {
      Xhatcand = cz_linear_state_estimator::execute(1, X0, Inputprev - Ueq, y, I, W, V);

      Xhatisempty = czonotope::isempty(Xhatcand);
      if (Xhatisempty)
      {
        Xhat = X0;
      }
      else
      {
        Xhat = Xhatcand;
      }
    }
    else  // State estimation
    {
      Xhatcand = cz_linear_state_estimator::execute(2, Xhat, Inputprev - Ueq, y, I, W, V);
      Xhatisempty = czonotope::isempty(Xhatcand);
      if (Xhatisempty)
      {
        Xhat = cz_linear_state_estimator::execute(0, Xhat, Inputprev - Ueq, y, I, W, V);
      }
      else
      {
        Xhat = Xhatcand;
      }
    }

    // Complexity reduction
    if (complexity_reduction)
    {
      Xhat = czonotope::creduction(Xhat, nc_max);
      Xhat = czonotope::greduction(Xhat, ng_max);
    }

    // Interval hull
    if (Xhat.A.rows() == 0)  // If Xhat is a zonotope
    {
      Eigen::MatrixXd ZonHullRad = Xhat.G.cwiseAbs().rowwise().sum();

      Eigen::MatrixXd ZonHull(23, 2);
      ZonHull << -ZonHullRad + Xhat.c, ZonHullRad + Xhat.c;

      Xhathull = ZonHull;
    }
    else
    {
      Xhathull = czonotope::intervalhull(Xhat);
    }
    MidHull = 0.5 * Xhathull.rowwise().sum();

    // Reshape
    Eigen::Map<Eigen::RowVectorXd> Xhathull_(Xhathull.data(), Xhathull.size());

    // Constrained Kalman filter
    if (count == 0)  // Initial estimate
    {
      CKF_xhat_plus = MidHull;
      CKF_P_plus = CKF_Pxx0;
      CKF_xhat_plus_prev = CKF_xhat_plus;
      CKF_P_plus_prev = CKF_P_plus;
    }
    else  // Constrained state estimation
    {
      constrained_kalman_filter::execute(2, CKF_xhat_plus_prev, CKF_P_plus_prev, Inputprev - Ueq, y, I, CKF_Q, CKF_R,
                                         Xhat, CKF_xhat_plus, CKF_P_plus);
      CKF_xhat_plus_prev = CKF_xhat_plus;
      CKF_P_plus_prev = CKF_P_plus;
    }

    // Convertendo velocidade angular
    std::vector<double> etadot = routines::wIIB2EtaDot(msgstates.values.at(13),  // wIIL_x
                                                       msgstates.values.at(14),  // wIIL_y
                                                       msgstates.values.at(15),  // wIIL_z
                                                       msgstates.values.at(3),   // phi
                                                       msgstates.values.at(4),   // theta
                                                       msgstates.values.at(5));  // psii

    // Integrador Trapezoidal (usando estados estimados pelo CKF)
    double x_atual = (CKF_xhat_plus + Xeq)(0) - Xref(0);
    xint = xint + (T / 2) * (x_atual + x_ant);
    x_ant = x_atual;
    double y_atual = (CKF_xhat_plus + Xeq)(1) - Xref(1);
    yint = yint + (T / 2) * (y_atual + y_ant);
    y_ant = y_atual;
    double z_atual = (CKF_xhat_plus + Xeq)(2) - Xref(2);
    zint = zint + (T / 2) * (z_atual + z_ant);
    z_ant = z_atual;
    double yaw_atual = (CKF_xhat_plus + Xeq)(5) - Xref(5);
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

    // State vector CKF (augmented)
    Eigen::MatrixXd Xstates_CKF(24, 1);
    Xstates_CKF << CKF_xhat_plus.block(0, 0, 20, 1),  // xhat
        xint, yint, zint, yawint;

    // Tracking error
    Erro = Xstates - Xref;

    // Tracking error (using estimation from CKF)
    Eigen::MatrixXd Erro_CKF(24, 1);
    Erro_CKF = (Xstates_CKF + (Eigen::MatrixXd(24, 1) << Xeq.block(0, 0, 20, 1), 0, 0, 0, 0).finished()) - Xref;

    Input << -K * Erro_CKF, dtrb;

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

    // Feedforward
    count++;

    // Storing values for the next iteration
    Inputprev << Input.block(0, 0, 4, 1);

    // States
    X << Xstates, Xhat.c, Xhathull_.transpose(), MidHull, Xhatisempty, y.size(), I.size(), Xhat.G.rows(), Xhat.G.cols(),
        Xhat.A.rows(), CKF_xhat_plus;

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
  return new vant2load_Franklin_CZL;
}
void destroy(Icontroller* p)
{
  delete p;
}
}
