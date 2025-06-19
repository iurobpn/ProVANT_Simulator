/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @brief This file contains the implementation of the Winfinity controller
 * with the vector of generalized coordinates split into q_c
 * (Controlled degrees of freedom), and q_s (Stabilized degrees of freedom)
 * for trajectory tracking in the helicopter flight mode of the UAV 4.0.
 *
 * @details The mathematical background used to design this controller is
 * presented in the Ph.D. Thesis "ROBUST CONTROL FRAMEWORK IN THE WEIGHTED
 *  SOBOLEV SPACE" available at https://www.ppgee.ufmg.br/bancodefesas.php.
 *
 * @author Daniel Cardoso
 */

#include <control_strategies_base/icontroller.hpp>
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include "simulator_msgs/Sensor.h"
#include <cmath>
#include <ros/ros.h>
#include "vant4_2_winfnonlinear/TiltrotorParameters.h"
#include "vant4_2_winfnonlinear/TiltrotorMatrices.h"
#include "vant4_2_winfnonlinear/TiltrotorInputCouplingMatrices.h"
#include "vant4_2_winfnonlinear/cplexutils.h"

class Vant42Winfnonlinear : public Icontroller
{
private:
  Eigen::MatrixXd K, M, C, V, invE, M1, B, Bo, BAero, Ba;

private:
  Eigen::VectorXd Input, Input2, Trajectory, Erro, X, Uref, Xref, G, q, qp, qpnew, Uopt, qdpp, M2, FUa;

private:
  double T, NSets, VelR, phiR, T_desR, DeltaRPS;

private:
  Eigen::VectorXd RPSA, qref, qrefp, qrefpp, ForcesFriction, SystemStates, EnvironmentWindI;

public:
  Vant42Winfnonlinear()
    : V(17, 17)
    , M1(17, 8)
    , qref(8)
    , qrefp(8)
    , invE(8, 8)
    , qrefpp(8)
    , ForcesFriction(8)
    , RPSA(2)
    , M(8, 8)
    , C(8, 8)
    , G(8)
    , q(8)
    , qp(8)
    , qpnew(8)
    , X(17)
    , Input(8)
    , Input2(4)
    , Trajectory(24)
    , SystemStates(16)
    , K(8, 20)
    , Uopt(8)
    , Bo(8, 4)
    , B(8, 8)
    , FUa(3)
    , BAero(8, 5)
    , Ba(8, 4)
    , M2(8)
    , qdpp(8)
    , EnvironmentWindI(3)
  {
    InitializeParameters();  // Initialize Tilt-rotor UAV parameters

    RPSA << 53.9, 53.9;  // Give initial RPS to improve computation of GetRPS function
    DeltaRPS = 10;       // This variable is a parameter in GetRPS function

    T = 0.012;  // Period of execution controller

    // Controller variables obtained by solving Riccati equations
    V << 18.974, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.44721, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1.3352, 0, 0, 0, 0, 0, 1.4142, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.3352, 5.3186e-16, 0, 0, 0, 0, 1.4142,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.3186e-16, 77.535, 0, 0, 0, 0, 0, 0.3873, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        525.21, 5.3931e-45, -1.1136e-60, 0, 0, 0, 75.842, 5.289e-45, -1.5979e-60, 5, 0, 0, 0, 0, 0, 0, 0, 5.3931e-45,
        525.21, 4.4824e-15, 0, 0, 0, 5.289e-45, 75.842, -6.2053e-15, 0, 5, 0, 0, 0, 0, 0, 0, -1.1136e-60, 4.4824e-15,
        116.1, 0, 0, 0, -1.5979e-60, -6.2053e-15, 34.8, 0, 0, 5, 0, 0, 1.4142, 0, 0, 0, 0, 0, 18.883, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1.4142, 0, 0, 0, 0, 0, 18.883, 7.5216e-15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3873, 0, 0, 0, 0,
        1.3733e-17, 2.0019, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 75.842, 5.289e-45, -1.5979e-60, 0, 0, 0, 74.666,
        6.3736e-45, -3.832e-60, 5.2521, 5.3931e-47, -1.1136e-62, 0, 0, 0, 0, 0, 5.289e-45, 75.842, -6.2053e-15, 0, 0, 0,
        6.3736e-45, 74.666, -1.3729e-14, 5.3931e-47, 5.2521, 4.4824e-17, 0, 0, 0, 0, 0, -1.5979e-60, -6.2053e-15, 34.8,
        0, 0, 0, -2.5204e-60, -3.3984e-15, 75.808, -1.1136e-61, 4.4824e-16, 11.61, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0,
        5.2521, 5.3931e-47, -1.1136e-62, 0.75842, 5.289e-47, -1.5979e-62, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 5.3931e-47,
        5.2521, 4.4824e-17, 5.289e-47, 0.75842, -6.2053e-17, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, -1.1136e-61, 4.4824e-16,
        11.61, -1.5979e-61, -6.2053e-16, 3.48;

    invE << 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0,
        0, 0, 0.066667, 0, 0, 0, 0, 0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, 0, 0, 0.02;
    M1 << invE, Eigen::MatrixXd::Zero(9, 8);
  }
  virtual ~Vant42Winfnonlinear()
  {
  }

  Eigen::VectorXd ForwardVelocitytoAlpha(double Ub, double Ubp)
  {
    Eigen::VectorXd out(2);

    if (Ub > 33.99)
    {
      Ub = 33.99;
    }
    else if (Ub < 0.0)
    {
      Ub = 0.0;
    }
    double Alpha = 3.0081e-08 * pow(Ub, 5) + 4.5551e-06 * pow(Ub, 4) - 0.00024816 * pow(Ub, 3) + 0.003874 * pow(Ub, 2) -
                   0.0061735 * Ub - 0.13061;
    double Alphap = (5 * 3.0081e-08 * pow(Ub, 4) + 4 * 4.5551e-06 * pow(Ub, 3) + 3 * (-0.00024816) * pow(Ub, 2) +
                     2 * 0.003874 * Ub + (-0.0061735)) *
                    Ubp;

    out << Alpha, Alphap;
    return out;
  }

  Eigen::MatrixXd SkewSymmetricMatrix(Eigen::VectorXd Vector)
  {
    // Place Vet in the Skew Symmetric matrix S
    Eigen::MatrixXd SkewMatrix(3, 3);
    SkewMatrix << 0, -Vector(2), Vector(1), Vector(2), 0, -Vector(0), -Vector(1), Vector(0), 0;

    return SkewMatrix;
  }

private:
  double sigmf(double Val, double Rate, double center)
  {
    return (1 / (pow(2.718281828459046, (-Rate * (Val - center))) + 1));
  }

private:
  double BumpFunction(double Theta)
  {
    double Val;
    if (abs(Theta) < 1.0)
    {
      Val = exp(-1.0 / (1.0 - pow(Theta, 2)));
    }
    else
    {
      Val = 0;
    }
    return Val;
  }

  Eigen::VectorXd MakeTrajectoryForwardAccelerationCircularDeceleration(double Tempo)
  {
    // Trajectory with Forward acceleration, circular path and deceleration
    Eigen::VectorXd Traj(12);

    // Trajectory parameters
    const double MaxVel = 33.0;         // Max Forward flight vel
    const double Acceleration = 0.5;    // Acceleration
    const double ZPos = 50;             // Max altitude
    const double ZPosF = 49;            // Max Landing position 49 = 1 m up to the ground
    const double TForwardFlight = 20;   // Time first forward flight
    const double TForwardFlight2 = 20;  // Time second forward flight
    const double TCircular = 80;        // Time circular path

    const double TMaxVel = MaxVel / Acceleration;
    const double R = MaxVel * TCircular / (2 * pi);
    double x, y, z, psi, xp, yp, zp, psip, xpp, ypp, zpp, psipp;

    if (Tempo <= TMaxVel)
    {
      x = (1.0 / 2.0) * Acceleration * pow(Tempo, 2);
      y = 0.0;
      z = (1.0 / 2.0) * (2.0 * ZPos / pow(TMaxVel, 2)) * pow(Tempo, 2) + 1.0;
      psi = 0.0;

      xp = Acceleration * Tempo;
      yp = 0.0;
      zp = (2.0 * ZPos / pow(TMaxVel, 2)) * Tempo;
      psip = 0.0;

      xpp = Acceleration;
      ypp = 0;
      zpp = (2.0 * ZPos / pow(TMaxVel, 2));
      psipp = 0;
    }
    else if (Tempo > TMaxVel && Tempo <= (TMaxVel + TForwardFlight))
    {
      Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(TMaxVel);
      Tempo = Tempo - TMaxVel;
      const double x0 = Traj(0);
      const double v0 = Traj(4);

      x = x0 + v0 * Tempo;
      y = 0;
      z = Traj(2);
      psi = 0;

      xp = v0;
      yp = 0;
      zp = 0;
      psip = 0;

      xpp = 0;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }
    else if (Tempo > (TMaxVel + TForwardFlight) && Tempo <= (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular))
    {
      Tempo = Tempo - (TMaxVel + TForwardFlight);
      Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(TMaxVel + TForwardFlight);

      const double xi = Traj(0);
      const double yi = Traj(1);
      const double zi = Traj(2);

      x = R * cos(2.0 * pi * Tempo / TCircular - pi / 2.0) + xi;
      y = R * sin(2.0 * pi * Tempo / TCircular - pi / 2.0) + R + yi;
      z = zi;
      psi = 2.0 * pi * Tempo / TCircular;

      xp = (2.0 * R * pi * sin(pi / 2.0 - (2.0 * pi * Tempo) / TCircular)) / TCircular;
      yp = (2.0 * R * pi * cos(pi / 2.0 - (2.0 * pi * Tempo) / TCircular)) / TCircular;
      zp = 0;
      psip = 2 * pi / TCircular;

      xpp = -(4.0 * pow(pi, 2) * R * cos((2.0 * pi * Tempo) / TCircular - pi / 2.0)) / pow(TCircular, 2);
      ypp = -(4.0 * pow(pi, 2) * R * sin((2.0 * pi * Tempo) / TCircular - pi / 2.0)) / pow(TCircular, 2);
      zpp = 0;
      psipp = 0;
    }
    else if (Tempo > (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular) &&
             Tempo <= (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular + TForwardFlight2))
    {
      Tempo = Tempo - (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular);
      Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular);

      x = Traj(0);
      y = Traj(1) + Traj(5) * Tempo;
      z = Traj(2);
      psi = Traj(3);

      xp = 0;
      yp = Traj(5);
      zp = 0;
      psip = 0;

      xpp = 0;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }
    else if (Tempo > (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular + TForwardFlight2) &&
             Tempo <= (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular + TForwardFlight2 + TMaxVel))
    {
      Tempo = Tempo - (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular + TForwardFlight2);
      Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular +
                                                                   TForwardFlight2);

      x = Traj(0);
      y = Traj(1) + Traj(5) * Tempo + (1.0 / 2.0) * Acceleration * pow(Tempo, 2);
      z = Traj(2) - (1.0 / 2.0) * (2 * ZPosF / pow(TMaxVel, 2)) * (pow(Tempo, 2));
      psi = Traj(3);

      xp = 0;
      yp = Traj(5) + Acceleration * Tempo;
      zp = -(2.0 * ZPosF / pow(TMaxVel, 2)) * Tempo;
      psip = 0;

      xpp = 0;
      ypp = Acceleration;
      zpp = -(2.0 * ZPosF / pow(TMaxVel, 2));
      psipp = 0;
    }
    else
    {
      Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(
          (TMaxVel + TForwardFlight + (3.0 / 4.0) * TCircular + TForwardFlight2 + TMaxVel));

      x = Traj(0);
      y = Traj(1);
      z = Traj(2);
      psi = Traj(3);

      xp = 0;
      yp = 0;
      zp = 0;
      psip = Traj(7);

      xpp = 0;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }

    Traj << x, y, z, psi, xp, yp, zp, psip, xpp, ypp, zpp, psipp;

    return Traj;
  }

  Eigen::VectorXd TrajectoryForwardAccelerationCircularDeceleration(double Tempo)
  {
    // This function normalize the function MakeTrajectoryForwardAccelerationCircularDeceleration to use in the
    // controller
    Eigen::VectorXd Traj(12), Trajectory(24), XpYpZp(3), XppYppZpp(3), UVW(3), UVWp(3), PhipThetapPsip(3);
    Eigen::MatrixXd RIB(3, 3), Wn(3, 3);

    Traj = MakeTrajectoryForwardAccelerationCircularDeceleration(Tempo);

    double phi = 0;
    double theta = 0;
    double psi = Traj(3);

    PhipThetapPsip << 0, 0, Traj(7);

    XpYpZp << Traj(4), Traj(5), Traj(6);

    Wn << 1.0, 0.0, -sin(theta), 0.0, cos(phi), cos(theta) * sin(phi), 0.0, -sin(phi), cos(phi) * cos(theta);

    RIB << (cos(psi) * cos(theta)), (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)),
        (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)), (cos(theta) * sin(psi)),
        (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)),
        (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)), (-sin(theta)), (cos(theta) * sin(phi)),
        (cos(phi) * cos(theta));

    UVW = RIB.transpose() * XpYpZp;
    UVWp = RIB.transpose() * XppYppZpp + SkewSymmetricMatrix(UVW) * Wn * PhipThetapPsip;

    Eigen::VectorXd out(2);
    out = ForwardVelocitytoAlpha(UVW(0), UVWp(0));

    double Alpha = out(0);
    double Alphap = out(1);

    Eigen::VectorXd qr(8);
    Eigen::VectorXd qpr(8);
    Eigen::VectorXd qppr(8);
    qr << Alpha, Alpha, 0, 0, Traj(3), Traj(0), Traj(1), Traj(2);
    qpr << Alphap, Alphap, 0, 0, Traj(7), Traj(4), Traj(5), Traj(6);
    qppr << 0, 0, 0, 0, Traj(11), Traj(8), Traj(9), Traj(10);

    Trajectory << qr, qpr, qppr;

    return Trajectory;
  }

  Eigen::VectorXd TrajetoriaReferenciaTransition(double Tempo)
  {
    // Trajectory with acceleration, forward fight, and decceleration
    Eigen::VectorXd Traj(24);

    // Trajectory parameters
    const double MaxVel = 33.0;       // Max forward velocity
    const double Acceleration = 2.5;  // Acceleration and decceleration
    double TForwardFlight = 10.0;     // Forward flight Time

    double TMaxVel = MaxVel / Acceleration;
    double x, y, z, psi, xp, yp, zp, psip, xpp;

    if (Tempo < TMaxVel)
    {
      x = (1.0 / 2.0) * Acceleration * pow(Tempo, 2.0);
      y = 0.0;
      z = 2.0;
      psi = 0.0;

      xp = Acceleration * Tempo;
      yp = 0.0;
      zp = 0.0;
      psip = 0.0;

      xpp = Acceleration;
    }
    else if (Tempo >= TMaxVel && Tempo <= (TMaxVel + TForwardFlight))
    {
      Tempo = Tempo - TMaxVel;

      x = (1.0 / 2.0) * Acceleration * pow(TMaxVel, 2.0) + Acceleration * TMaxVel * Tempo;
      y = 0.0;
      z = 2.0;
      psi = 0.0;

      xp = Acceleration * TMaxVel;
      yp = 0.0;
      zp = 0.0;
      psip = 0.0;

      xpp = 0.0;
    }
    else if (Tempo > (TMaxVel + TForwardFlight) && Tempo <= (2.0 * TMaxVel + TForwardFlight))
    {
      Tempo = Tempo - (TMaxVel + TForwardFlight);

      x = (1.0 / 2.0) * Acceleration * pow(TMaxVel, 2) + Acceleration * TMaxVel * (TForwardFlight) +
          Acceleration * TMaxVel * Tempo - (1.0 / 2.0) * Acceleration * pow(Tempo, 2);
      y = 0.0;
      z = 2.0;
      psi = 0.0;

      xp = Acceleration * TMaxVel - Acceleration * Tempo;
      yp = 0.0;
      zp = 0.0;
      psip = 0.0;

      xpp = -Acceleration;
    }
    else
    {
      double TempoFinal = TMaxVel;

      x = (1.0 / 2.0) * Acceleration * pow(TMaxVel, 2) + Acceleration * TMaxVel * (TForwardFlight) +
          Acceleration * TMaxVel * TempoFinal - (1.0 / 2.0) * Acceleration * pow(TempoFinal, 2);
      y = 0.0;
      z = 2.0;
      psi = 0.0;

      xp = 0.0;
      yp = 0.0;
      zp = 0.0;
      psip = 0.0;

      xpp = 0.0;
    }

    Eigen::VectorXd out(2);
    out = ForwardVelocitytoAlpha(xp, xpp);

    double Alpha = out(0);
    double Alphap = out(1);
    double Theta = 0;
    double Thetap = 0;

    Eigen::VectorXd qr(8);
    Eigen::VectorXd qpr(8);
    Eigen::VectorXd qppr(8);
    qr << Alpha, Alpha, 0, Theta, psi, x, y, z;
    qpr << Alphap, Alphap, 0, Thetap, psip, xp, yp, zp;
    qppr << 0, 0, 0, 0, 0, xpp, 0, 0;

    Traj << qr, qpr, qppr;

    return Traj;
  }

  Eigen::VectorXd TrajetoriaReferenciaHovering(double Tempo)
  {
    // Trajectory hovering
    Eigen::VectorXd Traj(24);
    Eigen::VectorXd qr(8);
    Eigen::VectorXd qpr(8);
    Eigen::VectorXd qppr(8);

    // Trajectory Hovering
    Eigen::VectorXd out(2);
    out = ForwardVelocitytoAlpha(0, 0);
    double Alpha = out(0);
    double Alphap = out(1);

    qr << Alpha, Alpha, 0, 0, 0, 0, 0, 2;
    qpr << Alphap, Alphap, 0, 0, 0, 0, 0, 0;
    qppr << 0, 0, 0, 0, 0, 0, 0, 0;

    Traj << qr, qpr, qppr;
    return Traj;
  }

  void config()
  {
  }

private:
  Eigen::VectorXd TrajetoriaReferenciaHelicoidal(double Tempo)
  {
    // Trajectory helicoidal
    Eigen::VectorXd Traj(12);  //[x y z psi xp yp zp psip xpp ypp zpp psipp]

    // The parameters of this trajectory can not be changed
    double Turns = 2;
    double T = 25;
    double f = 1 / T;
    double w = 2 * pi * f;
    double a = 20;
    double TakeoffVelocity = 0.5;
    double DeccelerationTime = 10.0;
    double LandingVelocity = -2;
    double LandingTime = -(1 + TakeoffVelocity * Turns * T) / LandingVelocity;
    double x, y, z, psi, xp, yp, zp, psip, xpp, ypp, zpp, psipp;

    if (Tempo <= Turns * T)
    {
      x = a * sin(w * Tempo);
      y = a * cos(w * Tempo);
      z = TakeoffVelocity * Tempo + 1;
      psi = 0;

      xp = a * w * cos(Tempo * w);
      yp = -a * w * sin(Tempo * w);
      zp = TakeoffVelocity;
      psip = 0;

      xpp = -a * pow(w, 2) * sin(Tempo * w);
      ypp = -a * pow(w, 2) * cos(Tempo * w);
      zpp = 0;
      psipp = 0;
    }
    else if (Tempo > Turns * T && Tempo <= (Turns * T + DeccelerationTime))
    {
      Tempo = Tempo - Turns * T;
      Traj = TrajetoriaReferenciaHelicoidal(Turns * T);
      double x0 = Traj(0);
      double y0 = Traj(1);
      double z0 = Traj(2);
      double psi0 = Traj(3);
      double xp0 = Traj(4);

      double Decceleration = -xp0 / DeccelerationTime;
      x = x0 + xp0 * Tempo + 0.5 * Decceleration * pow(Tempo, 2);
      y = y0;
      z = z0;
      psi = psi0;

      xp = xp0 + Decceleration * Tempo;
      yp = 0;
      zp = 0;
      psip = 0;

      xpp = Decceleration;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }
    else if (Tempo > (Turns * T + DeccelerationTime) && Tempo <= (Turns * T + DeccelerationTime + LandingTime))
    {
      Tempo = Tempo - (Turns * T + DeccelerationTime);
      Traj = TrajetoriaReferenciaHelicoidal(Turns * T + DeccelerationTime);

      double x0 = Traj(0);
      double y0 = Traj(1);
      double z0 = Traj(2);
      double psi0 = Traj(3);

      x = x0;
      y = y0;
      z = z0 + LandingVelocity * Tempo;
      psi = psi0;

      xp = 0;
      yp = 0;
      zp = LandingVelocity;
      psip = 0;

      xpp = 0;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }
    else
    {
      Traj = TrajetoriaReferenciaHelicoidal(Turns * T + DeccelerationTime + LandingTime);

      double x0 = Traj(0);
      double y0 = Traj(1);
      double z0 = Traj(2);
      double psi0 = Traj(3);

      x = x0;
      y = y0;
      z = z0;
      psi = psi0;

      xp = 0;
      yp = 0;
      zp = 0;
      psip = 0;

      xpp = 0;
      ypp = 0;
      zpp = 0;
      psipp = 0;
    }
    Traj << x, y, z, psi, xp, yp, zp, psip, xpp, ypp, zpp, psipp;

    return Traj;
  }

private:
  Eigen::VectorXd MakeTrajetoriaReferenciaHovering(double Tempo)
  {
    // Make trajectory to be used by the control strategy
    Eigen::VectorXd Traj(12), Trajectory(24);

    Traj = TrajetoriaReferenciaHelicoidal(Tempo);
    Eigen::VectorXd out(2);
    out = ForwardVelocitytoAlpha(Traj(5), Traj(8));

    double Alpha = out(0);
    double Alphap = out(1);
    double Theta = 0;
    double Thetap = 0;

    Eigen::VectorXd qr(8);
    Eigen::VectorXd qpr(8);
    Eigen::VectorXd qppr(8);
    qr << Alpha, Alpha, 0, Theta, Traj(3), Traj(0), Traj(1), Traj(2) + 1.2;
    qpr << Alphap, Alphap, 0, Thetap, Traj(7), Traj(4), Traj(5), Traj(6);
    qppr << 0, 0, 0, 0, Traj(11), Traj(8), Traj(9), Traj(10);

    Trajectory << qr, qpr, qppr;  // ar al phi the psi x y z

    return Trajectory;
  }

  std::vector<double> execute(simulator_msgs::SensorArray arraymsg)
  {
    simulator_msgs::Sensor msg;
    bool found = false;
    for (int i = 0; i < arraymsg.values.size(); i++)
    {
      if (arraymsg.values.at(i).name == "Estados")
      {
        msg = arraymsg.values.at(i);
        found = true;
        break;
      }
    }

    if (!found)
    {
      // In case of error, report the problem in a ROS_LOG, and returns an
      // empty array
      ROS_FATAL("[adaptivelqrVant5Aerod] State vector not found.");
      std::vector<double> out(Input.data(), Input.data() + Input.size());
      return out;
    }

    // msg = [x y z phi theta psi ar al xp yp zp phip thetap psip arp alp]

    static double Tempo = 0;  // Initialize time variable
    Tempo = Tempo + T;        // Compute the instant of time

    Trajectory = TrajectoryForwardAccelerationCircularDeceleration(Tempo);

    // Anothe examples of trajectory
    // Trajectory = MakeTrajetoriaReferenciaHovering(Tempo);
    // Trajectory = TrajetoriaReferenciaTransition(Tempo);

    // For better performance the controller must known the environment wind magnitudes
    EnvironmentWindI << 0, 0, 0;  /// Environment wind

    // If an environment wind is generated in the aerodynamic plugin, this environment wind must be provided here
    // EnvironmentWindI << 0, 3, 0; //Helicopter
    // EnvironmentWindI << -3, 0, 0; //Full-flight envelope

    double phi = msg.values.at(3);
    double theta = msg.values.at(4);
    double psi = msg.values.at(5);

    if (psi - Trajectory(4) < -pi)
    {
      psi = psi + 2 * pi;
    }
    else if (psi - Trajectory(4) > pi)
    {
      psi = psi - 2 * pi;
    }

    Eigen::VectorXd XYZB(3), XYZBnew(3), XpYpZpB(3), XpYpZpBnew(3), PhipThetapPsip(3), q(8), qp(8);
    Eigen::MatrixXd RIB(3, 3), Wn(3, 3);

    RIB << (cos(psi) * cos(theta)), (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)),
        (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)), (cos(theta) * sin(psi)),
        (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)),
        (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)), (-sin(theta)), (cos(theta) * sin(phi)),
        (cos(phi) * cos(theta));

    Wn << 1.0, 0.0, -sin(theta), 0.0, cos(phi), cos(theta) * sin(phi), 0.0, -sin(phi), cos(phi) * cos(theta);

    PhipThetapPsip << msg.values.at(11), msg.values.at(12), msg.values.at(13);

    XYZB << msg.values.at(0), msg.values.at(1), msg.values.at(2);
    XYZBnew = XYZB + RIB * TranslateBodyFrame;
    XpYpZpB << msg.values.at(8), msg.values.at(9), msg.values.at(10);
    XpYpZpBnew = XpYpZpB - RIB * SkewSymmetricMatrix(TranslateBodyFrame) * Wn * PhipThetapPsip;

    q << msg.values.at(6),  // ar
        msg.values.at(7),   // al
        msg.values.at(3),   // phi
        msg.values.at(4),   // theta
        psi,                // psi
        XYZBnew;            // xyz

    qp << msg.values.at(14),  // arp
        msg.values.at(15),    // alp
        msg.values.at(11),    // phip
        msg.values.at(12),    // thetap
        msg.values.at(13),    // psip
        XpYpZpBnew;           // xyzp
    SystemStates << q, qp;

    qref << Trajectory(0), Trajectory(1), Trajectory(2), Trajectory(3), Trajectory(4), Trajectory(5), Trajectory(6),
        Trajectory(7);
    qrefp << Trajectory(8), Trajectory(9), Trajectory(10), Trajectory(11), Trajectory(12), Trajectory(13),
        Trajectory(14), Trajectory(15);
    qrefpp << Trajectory(16), Trajectory(17), Trajectory(18), Trajectory(19), Trajectory(20), Trajectory(21),
        Trajectory(22), Trajectory(23);

    // integrators variables
    static double psiint = 0, psi_ant = 0;
    static double xint = 0, x_ant = 0;
    static double yint = 0, y_ant = 0;
    static double zint = 0, z_ant = 0;

    double psi_error = q(4) - qref(4);
    psiint = psiint + (T / 2.0) * (psi_error + psi_ant);
    psi_ant = psi_error;

    double x_error = q(5) - qref(5);
    xint = xint + (T / 2.0) * (x_error + x_ant);
    x_ant = x_error;

    double y_error = q(6) - qref(6);
    yint = yint + (T / 2.0) * (y_error + y_ant);
    y_ant = y_error;

    double z_error = q(7) - qref(7);
    zint = zint + (T / 2.0) * (z_error + z_ant);
    z_ant = z_error;

    M = InertiaMatrix(q);
    C = coriolisMatrix(q, qp);
    G = GravitationVector(q);
    double VPR, phiPR, VPL, phiPL, ub;
    ;
    Bo = InputCouplingMatrix(q, qp, RPSA, &VPR, &phiPR, &VPL, &phiPL, EnvironmentWindI);
    BAero = InputCouplingMatrixAero(q, qp, &ub, EnvironmentWindI);
    Ba << BAero.col(0), BAero.col(1), BAero.col(2), BAero.col(3);
    FUa = BAero.col(4);
    ForcesFriction << 0.0, 0.0, 0.005 * msg.values.at(14), 0.005 * msg.values.at(15), 0.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd qsp(2), qrtilp(3), qctilp(3), qrtil(3), qctil(3), intqctil(3);
    qsp << qp(2), qp(3);
    qrtilp << qp(0) - qrefp(0), qp(1) - qrefp(1), qp(4) - qrefp(4);
    qctilp << qp(5) - qrefp(5), qp(6) - qrefp(6), qp(7) - qrefp(7);
    qrtil << q(0) - qref(0), q(1) - qref(1), q(4) - qref(4);
    qctil << q(5) - qref(5), q(6) - qref(6), q(7) - qref(7);
    intqctil << xint, yint, zint;

    X << qsp, qrtilp, qctilp, qrtil, qctil, intqctil;  // zp

    qpnew << qp(2), qp(3), qp(0), qp(1), qp(4), qp(5), qp(6), qp(7);

    qdpp << qrefpp(2), qrefpp(3), qrefpp(0), qrefpp(1), qrefpp(4), qrefpp(5), qrefpp(6), qrefpp(7);
    M2 << Eigen::VectorXd::Zero(2), qdpp(2), qdpp(3), qdpp(4), qdpp(5), qdpp(6), qdpp(7);

    Uopt = -M * M1.transpose() * V * X + G + M * M2 + C * qpnew - FUa + ForcesFriction;

    // input coupling matrix
    B << Bo, Ba;

    Eigen::MatrixXd H(8, 8);
    Eigen::MatrixXd f(8, 1);
    Eigen::MatrixXd A(0, 8);
    Eigen::MatrixXd b(0, 1);
    Eigen::MatrixXd Aeq(0, 8);
    Eigen::MatrixXd beq(0, 1);
    Eigen::MatrixXd LB(8, 1);
    Eigen::MatrixXd UB(8, 1);

    H = B.transpose() * B;
    f = -B.transpose() * Uopt;

    double AngleDeflection = 0.436332312998582;
    double Transition = 12;
    UB << 90, 90, 2, 2, AngleDeflection * sigmf(ub, 2.0, Transition), AngleDeflection * sigmf(ub, 2.0, Transition),
        AngleDeflection * sigmf(ub, 2.0, Transition),
        AngleDeflection * sigmf(ub, 2.0, Transition);  // 0.52,  0.52,  0.52,  0.52;
    LB << 10, 10, -2, -2, -AngleDeflection * sigmf(ub, 2.0, Transition), -AngleDeflection * sigmf(ub, 2.0, Transition),
        -AngleDeflection * sigmf(ub, 2.0, Transition),
        -AngleDeflection * sigmf(ub, 2.0, Transition);  //-0.52, -0.52, -0.52, -0.52;

    double objvalue;
    Eigen::MatrixXd xstar(8, 1);
    int exitflag;

    cplexutils::qp(H, f, A, b, Aeq, beq, LB, UB, objvalue, xstar, exitflag);

    Input(0) = xstar(0);
    Input(1) = xstar(1);
    Input(2) = xstar(2);
    Input(3) = xstar(3);
    Input(4) = xstar(4);
    Input(5) = xstar(5);
    Input(6) = xstar(6);
    Input(7) = xstar(7);

    Input(0) = GetRPS(VPR, phiPR, Input(0), 1, 100);  // Converte propeller's force to RPS
    Input(1) = GetRPS(VPL, phiPL, Input(1), 1, 100);  // Converte propeller's force to RPS
    RPSA(0) = Input(0);                               // Save old RPS to use as initial condition to the solver
    RPSA(1) = Input(1);                               // Save old RPS to use as initial condition to the solver

    std::vector<double> out(Input.data(), Input.data() + Input.rows() * Input.cols());
    return out;
  }

  std::vector<double> Reference()
  {
    std::vector<double> out(Trajectory.data(), Trajectory.data() + Trajectory.rows() * Trajectory.cols());
    return out;
  }

  std::vector<double> Error()
  {
    std::vector<double> out(EnvironmentWindI.data(),
                            EnvironmentWindI.data() + EnvironmentWindI.rows() * EnvironmentWindI.cols());
    return out;
  }

  std::vector<double> State()
  {
    std::vector<double> out(SystemStates.data(), SystemStates.data() + SystemStates.rows() * SystemStates.cols());
    return out;
  }
};

extern "C" {
Icontroller* create(void)
{
  return new Vant42Winfnonlinear;
}
void destroy(Icontroller* p)
{
  delete p;
}
}
