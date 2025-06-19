/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @brief This file contains the declaration for the function that returns the
 * inertia matrix for the UAV 4.0.
 * @author Daniel Cardoso
 */

#include <eigen3/Eigen/Eigen>
//#include "TiltrotorParameters.h"
#include <cmath>

double GetCt(double J, double phi)
{
  const double J_min = 0.0;
  const double J_max = 0.75;
  const double phi_min = -pi / 2.0;
  const double phi_max = pi / 2.0;
  double val;

  if (J < J_min)
  {
    J = J_min;
  }
  else if (J > J_max)
  {
    J = J_max;
  }

  if (phi < phi_min)
  {
    phi = phi_min;
  }
  else if (phi > phi_max)
  {
    phi = phi_max;
  }

  phi = abs(phi);
  Eigen::VectorXd Ct(15);
  Eigen::VectorXd Coef(15);
  Ct << 1, J, phi, pow(J, 2), J * phi, pow(phi, 2), pow(J, 3), pow(J, 2) * phi, J * pow(phi, 2), pow(phi, 3), pow(J, 4),
      pow(J, 3) * phi, pow(J, 2) * pow(phi, 2), J * pow(phi, 3), pow(phi, 4);
  Coef << 0.0829, 0.0066, 0.0, -0.2197, -0.1504, 0.0, 0.0630, 0.2064, 0.2217, 0.0, 0.0302, -0.1115, 0.0173, -0.0856,
      0.0;
  val = Ct.transpose() * Coef;
  return val;
}

double GetCq(double J, double phi)
{
  const double J_min = 0.0;
  const double J_max = 0.75;
  const double phi_min = -pi / 2.0;
  const double phi_max = pi / 2.0;
  double val;

  if (J < J_min)
  {
    J = J_min;
  }
  else if (J > J_max)
  {
    J = J_max;
  }

  if (phi < phi_min)
  {
    phi = phi_min;
  }
  else if (phi > phi_max)
  {
    phi = phi_max;
  }

  phi = abs(phi);
  Eigen::VectorXd Cq(15);
  Eigen::VectorXd Coef(15);

  Cq << 1, J, phi, pow(J, 2), J * phi, pow(phi, 2), pow(J, 3), pow(J, 2) * phi, J * pow(phi, 2), pow(phi, 3), pow(J, 4),
      pow(J, 3) * phi, pow(J, 2) * pow(phi, 2), J * pow(phi, 3), pow(phi, 4);
  Coef << 0.0053, 0.0074, 0.0, -0.0078, -0.0190, 0.0, -0.0215, 0.0254, 0.0183, 0.0, 0.0112, -0.0024, -0.0049, -0.0068,
      0.0;
  val = Cq.transpose() * Coef;
  return val;
}

Eigen::MatrixXd InputCouplingMatrix(Eigen::VectorXd q, Eigen::VectorXd qp, Eigen::VectorXd RPSA, double* VR,
                                    double* phiR, double* VL, double* phiL, Eigen::VectorXd EnvironmentWindI)
{
  // Variáveis de Configuração
  double AlphaR = q(0);
  double AlphaL = q(1);
  double Phi = q(2);
  double Theta = q(3);
  double Psi = q(4);
  double X = q(5);
  double Y = q(6);
  double Z = q(7);

  Eigen::MatrixXd Bo(8, 4);
  Eigen::MatrixXd Bonew(8, 4);
  Eigen::MatrixXd Wn(3, 3);
  Eigen::VectorXd az(3), ay(3), ayR(3), ayL(3);
  Eigen::MatrixXd RB2(3, 3), RB3(3, 3), RIB(3, 3), RICR(3, 3), RICL(3, 3);
  Eigen::MatrixXd JvRp(3, 8), JvLp(3, 8), JwRp(3, 8), JwLp(3, 8), JwCm(3, 8);
  Eigen::VectorXd dPIpR(3), UVWpR(3), dPIpL(3), UVWpL(3);
  // Configuration Matrices
  Wn << 1, 0, -sin(Theta), 0, cos(Phi), sin(Phi) * cos(Theta), 0, -sin(Phi), cos(Phi) * cos(Theta);

  RIB = Rotz(Psi) * Roty(Theta) * Rotx(Phi);
  RB2 = Roty(AlphaR) * Rotx(-B);
  RB3 = Roty(AlphaL) * Rotx(B);
  RICR = RIB * RB2;
  RICL = RIB * RB3;

  ay << 0, 1, 0;
  az << 0, 0, 1;

  ayR = Rotx(-B).transpose() * ay;
  ayL = Rotx(B).transpose() * ay;

  JvRp << -RIB * RB2 * SkewSymmetricMatrix(DA2P2) * ayR, Eigen::VectorXd::Zero(3),
      -RIB * SkewSymmetricMatrix(DNBA2 + RB2 * DA2P2) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JvLp << Eigen::VectorXd::Zero(3), -RIB * RB3 * SkewSymmetricMatrix(DA3P3) * ayL,
      -RIB * SkewSymmetricMatrix(DNBA3 + RB3 * DA3P3) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JwRp << RICR * ayR, Eigen::VectorXd::Zero(3), RIB * Wn, Eigen::MatrixXd::Zero(3, 3);
  JwLp << Eigen::VectorXd::Zero(3), RICL * ayL, RIB * Wn, Eigen::MatrixXd::Zero(3, 3);
  JwCm << Eigen::VectorXd::Zero(3), Eigen::VectorXd::Zero(3), RIB * Wn, Eigen::MatrixXd::Zero(3, 3);

  dPIpR = JvRp * qp;
  UVWpR = (RICR.transpose()) * (dPIpR - EnvironmentWindI);
  double VpRxz = pow(pow(UVWpR(2), 2) + pow(UVWpR(0), 2), 0.5);  // problema aqui
  double AlphapR = atan2(UVWpR(0), UVWpR(2));
  *VR = VpRxz;
  *phiR = AlphapR;

  dPIpL = JvLp * qp;
  UVWpL = (RICL.transpose()) * (dPIpL - EnvironmentWindI);
  double VpLxz = pow(pow(UVWpL(2), 2) + pow(UVWpL(0), 2), 0.5);  // problema aqui
  double AlphapL = atan2(UVWpL(0), UVWpL(2));
  *VL = VpLxz;
  *phiL = AlphapL;

  Bo << (JvRp.transpose() + JwRp.transpose() * DiameterPropeller *
                                (GetCq(VpRxz / (DiameterPropeller * RPSA(0)), AlphapR) /
                                 GetCt(VpRxz / (DiameterPropeller * RPSA(0)), AlphapR))) *
            RICR * az,
      (JvLp.transpose() - JwLp.transpose() * DiameterPropeller *
                              (GetCq(VpLxz / (DiameterPropeller * RPSA(1)), AlphapL) /
                               GetCt(VpLxz / (DiameterPropeller * RPSA(1)), AlphapL))) *
          RICL * az,
      (JwRp - JwCm).transpose() * RICR * ay, (JwLp - JwCm).transpose() * RICL * ay;

  Bonew << Bo(2, 0), Bo(2, 1), Bo(2, 2), Bo(2, 3), Bo(3, 0), Bo(3, 1), Bo(3, 2), Bo(3, 3), Bo(0, 0), Bo(0, 1), Bo(0, 2),
      Bo(0, 3), Bo(1, 0), Bo(1, 1), Bo(1, 2), Bo(1, 3), Bo(4, 0), Bo(4, 1), Bo(4, 2), Bo(4, 3), Bo(5, 0), Bo(5, 1),
      Bo(5, 2), Bo(5, 3), Bo(6, 0), Bo(6, 1), Bo(6, 2), Bo(6, 3), Bo(7, 0), Bo(7, 1), Bo(7, 2), Bo(7, 3);
  return Bonew;
}

#define NUM_BISECTION_ITERATIONS (50)
#define NUM_INTERVALS (50)

double f(double V, double phi, double T_des, double n)
{
  double Ti_eng = GetCt(V / (n * DiameterPropeller), phi) * Ho * pow(n, 2) * pow(DiameterPropeller, 4);
  return (T_des - Ti_eng);
}

double bisection_method(double V, double phi, double T_des, double x0, double x1)
{
  for (unsigned int i = 0; i < NUM_BISECTION_ITERATIONS; i++)
  {
    double midpoint = 0.5 * x0 + 0.5 * x1;
    f(V, phi, T_des, x0) * f(V, phi, T_des, midpoint) < 0 ? x1 = midpoint : x0 = midpoint;
  }
  return 0.5 * x0 + 0.5 * x1;
}

double GetRPS(double V, double phi, double T_des, double Min, double Max)
{
  double x0, x1;
  for (unsigned int i = 0; i < NUM_INTERVALS - 1; i++)
  {
    x0 = Min + (-(Min - Max) / NUM_INTERVALS) * (i);
    x1 = Min + (-(Min - Max) / NUM_INTERVALS) * (i + 1);
    if (f(V, phi, T_des, x0) * f(V, phi, T_des, x1) < 0)
    {
      break;
      // std::cout << "Optimal Inverval: " << "[" << x0 << x1 << "]" << std::endl;
    }
  }
  // std::cout << "Optimal value: " << bisection_method(V, phi, T_des, x0, x1) << std::endl;
  return bisection_method(V, phi, T_des, x0, x1);
}

double Saturation(double Val, double Valmin, double Valmax)
{
  if (Val < Valmin)
  {
    return Valmin;
  }
  else if (Val > Valmax)
  {
    return Valmax;
  }
  else
  {
    return Val;
  }
}

double asin2(double x, double y)
{
  if (x == 0)
  {
    return 0;
  }
  else
  {
    return asin(x / y);
  }
}

Eigen::MatrixXd InputCouplingMatrixAero(Eigen::VectorXd q, Eigen::VectorXd qp, double* ub,
                                        Eigen::VectorXd EnvironmentWindI)
{
  // Variáveis de Configuração
  double AlphaR = q(0);
  double AlphaL = q(1);
  double Phi = q(2);
  double Theta = q(3);
  double Psi = q(4);
  double X = q(5);
  double Y = q(6);
  double Z = q(7);

  // Variáveis de Configuração
  double AlphaRp = qp(0);
  double AlphaLp = qp(1);
  double Phip = qp(2);
  double Thetap = qp(3);
  double Psip = qp(4);
  double Xp = qp(5);
  double Yp = qp(6);
  double Zp = qp(7);

  Eigen::MatrixXd RI_B(3, 3);
  Eigen::MatrixXd Wn(3, 3);
  Eigen::MatrixXd RBFw(3, 3);
  Eigen::MatrixXd RBwRw(3, 3);
  Eigen::MatrixXd RBwLw(3, 3);
  Eigen::MatrixXd RBtRw(3, 3);
  Eigen::MatrixXd RBtLw(3, 3);
  Eigen::MatrixXd RBcGw(3, 3);
  Eigen::MatrixXd Jvf(3, 8);
  Eigen::MatrixXd JvwR(3, 8);
  Eigen::MatrixXd JvwL(3, 8);
  Eigen::MatrixXd JvtR(3, 8);
  Eigen::MatrixXd JvtL(3, 8);
  Eigen::MatrixXd Jvcg(3, 8);
  Eigen::MatrixXd Jw(3, 8);

  // computing transformation matrices
  RI_B << Rotz(Psi) * Roty(Theta) * Rotx(Phi);
  // std::cout << "Phi" << Phi << std::endl;
  // std::cout << "Theta" << Theta << std::endl;
  // std::cout << "Psi" << Psi << std::endl;
  // std::cout << "RI_B" << RI_B << std::endl;
  // std::cout << "Rotz(Psi)" << Rotz(Psi) << std::endl;

  Wn << 1.0, 0.0, -sin(Theta), 0.0, cos(Phi), cos(Theta) * sin(Phi), 0.0, -sin(Phi), cos(Phi) * cos(Theta);

  // Velocities Aerodinamic centers
  Eigen::VectorXd dPI_CG(3);
  Eigen::VectorXd dPI_f(3);
  Eigen::VectorXd dPI_wr(3);
  Eigen::VectorXd dPI_wl(3);
  Eigen::VectorXd dPI_tr(3);
  Eigen::VectorXd dPI_tl(3);
  // aerodynamic center velocities and wind velocities
  Eigen::VectorXd UVWf(3);
  Eigen::VectorXd UVWCG(3);
  Eigen::VectorXd UVWwr(3);
  Eigen::VectorXd UVWwl(3);
  Eigen::VectorXd UVWtr(3);
  Eigen::VectorXd UVWtl(3);
  Eigen::VectorXd FUa(8);
  Eigen::VectorXd FUaFuselage(8);
  Eigen::VectorXd FUaWings(8);
  Eigen::VectorXd FUaTail(8);
  Eigen::VectorXd FUaInterference(8);
  Eigen::VectorXd FUaAeroMoments(8);
  Eigen::MatrixXd BAero(8, 4);
  Eigen::MatrixXd BAeroFUa(8, 5);
  Eigen::VectorXd ay(3);
  Eigen::VectorXd az(3);
  ay << 0, 1, 0;
  az << 0, 0, 1;

  // std::cout << "PosCG" << PosCGNr << std::endl;
  // std::cout << "DBfNr" << DBfNr << std::endl;
  // std::cout << "DBwrNr" << DBwrNr << std::endl;
  // std::cout << "DBwlNr" << DBwlNr << std::endl;
  // std::cout << "DBtrNr" << DBtrNr << std::endl;
  // std::cout << "DBtlNr" << DBtlNr << std::endl;

  Eigen::VectorXd PhipThetapPsip(3), XpYpZp(3);

  PhipThetapPsip << Phip, Thetap, Psip;
  //-----------Computing [Xdot Ydot Zdot]-------------------------------%
  XpYpZp << Xp, Yp, Zp;

  Jvf << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(DNBf) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JvwR << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(DNBwr) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JvwL << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(DNBwl) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JvtR << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(DNBtr) * Wn, Eigen::MatrixXd::Identity(3, 3);
  JvtL << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(DNBtl) * Wn, Eigen::MatrixXd::Identity(3, 3);
  Jvcg << Eigen::MatrixXd::Zero(3, 2), -RI_B * SkewSymmetricMatrix(PosNBCG) * Wn, Eigen::MatrixXd::Identity(3, 3);
  Jw << Eigen::MatrixXd::Zero(3, 2), RI_B * Wn, Eigen::MatrixXd::Zero(3, 3);

  // Compute the velocity of the aerodynamic centers expressed in the Inertial frame
  dPI_CG << Jvcg * qp;  // velocity of the aerodynamic center of fuselage w.r.t I expressed in I
  dPI_f << Jvf * qp;    // velocity of the aerodynamic center of fuselage w.r.t I expressed in I
  dPI_wr << JvwR * qp;  // velocity of the aerodynamic center of wing right w.r.t I expressed in I
  dPI_wl << JvwL * qp;  // velocity of the aerodynamic center of wing left w.r.t I expressed in I
  dPI_tr << JvtR * qp;  // velocity of the aerodynamic center of tail right w.r.t I expressed in I
  dPI_tl << JvtL * qp;  // velocity of the aerodynamic center of tail left w.r.t I expressed in I

  //----------Computing Properties of Relative wind for center of gravity---------%

  UVWCG = RI_B.transpose() * (dPI_CG - EnvironmentWindI);  // Express the relative airspeed velocity on the frame
                                                           // positioned at the aerodynamic center of fuselage
  double VCG = pow(pow(UVWCG(2), 2) + pow(UVWCG(1), 2) + pow(UVWCG(0), 2), 0.5);  // Magnitude x-z axis
  *ub = VCG;
  double AlphaCG = atan2(UVWCG(2), UVWCG(0));  // Orientation - Angle of attack CG
  double BetaCG = asin2(UVWCG(1), VCG);        // Orientation - Side slip angle CG
  // std::cout << "VCG" << VCG << std::endl;
  // std::cout << "AlphaCG" << AlphaCG << std::endl;
  // std::cout << "BetaCG" << BetaCG << std::endl;
  //----------Computing Properties of Relative wind for fuselage---------%

  UVWf = RI_B.transpose() * (dPI_f - EnvironmentWindI);  // Express the relative airspeed velocity on the frame
                                                         // positioned at the aerodynamic center of fuselage
  double Vf = pow(pow(UVWf(2), 2) + pow(UVWf(1), 2) + pow(UVWf(0), 2), 0.5);  // Magnitude x-z axis
  double Alphaf = atan2(UVWf(2), UVWf(0));                                    // Orientation - Angle of attack fuselage
  double Betaf = asin2(UVWf(1), Vf);                                          // Orientation - Side slip angle fuselage
  // std::cout << "Vf" << Vf << std::endl;
  // std::cout << "Alphaf" << Alphaf << std::endl;
  // std::cout << "Betaf" << Betaf << std::endl;

  //----------Computing Properties of Relative wind for wings---------%
  // Wing right
  UVWwr = RI_B.transpose() * (dPI_wr - EnvironmentWindI);  // Express the velocity on the frame positioned at the
                                                           // aerodynamic center of wing R
  double Vwr = pow(pow(UVWwr(2), 2) + pow(UVWwr(1), 2) + pow(UVWwr(0), 2), 0.5);  // compute Magnitude
  double Alphawr = atan2(UVWwr(2), UVWwr(0));                                     // compute Orientation
  double Betawr = asin2(UVWwr(1), Vwr);
  // std::cout << "Vwr" << Vwr << std::endl;
  // std::cout << "Alphawr" << Alphawr << std::endl;
  // std::cout << "Betawr" << Betawr << std::endl;

  // Wing left
  UVWwl = RI_B.transpose() * (dPI_wl - EnvironmentWindI);  // Express the velocity on the frame positioned at the
                                                           // aerodynamic center of wing L
  double Vwl = pow(pow(UVWwl(2), 2) + pow(UVWwl(1), 2) + pow(UVWwl(0), 2), 0.5);  // compute Magnitude
  double Alphawl = atan2(UVWwl(2), UVWwl(0));                                     // compute Orientation
  double Betawl = asin2(UVWwl(1), Vwl);
  // std::cout << "Vwl" << Vwl << std::endl;
  // std::cout << "Alphawl" << Alphawl << std::endl;
  // std::cout << "Betawl" << Betawl << std::endl;

  //----------Computing Properties of Relative wind for V-tail---------%

  // Tail right
  UVWtr = RI_B.transpose() * (dPI_tr - EnvironmentWindI);  // Express the velocity on the frame positioned at the
                                                           // aerodynamic center of Tail R
  double Vtr = pow(pow(UVWtr(2), 2) + pow(UVWtr(1), 2) + pow(UVWtr(0), 2), 0.5);  // compute Magnitude
  double Alphatr = atan2(UVWtr(2), UVWtr(0));                                     // compute Orientation
  double Betatr = asin2(UVWtr(1), Vtr);
  // std::cout << "Vtr" << Vtr << std::endl;
  // std::cout << "Alphatr" << Alphatr << std::endl;
  // std::cout << "Betatr" << Betatr << std::endl;

  // Tail left
  UVWtl = RI_B.transpose() * (dPI_tl - EnvironmentWindI);  // Express the velocity on the frame positioned at the
                                                           // aerodynamic center of Tail L
  double Vtl = pow(pow(UVWtl(2), 2) + pow(UVWtl(1), 2) + pow(UVWtl(0), 2), 0.5);  // compute Magnitude
  double Alphatl = atan2(UVWtl(2), UVWtl(0));                                     // compute Orientation
  double Betatl = asin2(UVWtl(1), Vtl);
  // std::cout << "Vtl" << Vtl << std::endl;
  // std::cout << "Alphatl" << Alphatl << std::endl;
  // std::cout << "Betatl" << Betatl << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Alphaf = Saturation(-Alphaf, -1.55, 1.55);
  Betaf = Saturation(-Betaf, -1.55, 1.55);
  Alphawr = Saturation(-Alphawr, -1.55, 1.55);
  Betawr = Saturation(-Betawr, -1.55, 1.55);
  Alphawl = Saturation(-Alphawl, -1.55, 1.55);
  Betawl = Saturation(-Betawl, -1.55, 1.55);
  Alphatr = Saturation(-Alphatr, -1.55, 1.55);
  Betatr = Saturation(-Betatr, -1.55, 1.55);
  Alphatl = Saturation(-Alphatl, -1.55, 1.55);
  Betatl = Saturation(-Betatl, -1.55, 1.55);
  AlphaCG = Saturation(-AlphaCG, -1.55, 1.55);
  BetaCG = Saturation(-BetaCG, -1.55, 1.55);

  // The vector of aerodynamic coefficients are given in the range -90 to 90 deg with data for each 0.008726646259972
  // rad/s

  double DeltaCoefRad = 0.008726646259972;

  // calculo dos coef. da fuselagem
  double Vala = ((Alphaf + pi / 2.0) / DeltaCoefRad);
  int Indexa = floor(Vala);
  double Proporcaoa = Vala - Indexa;

  double Vals = ((Betaf + pi / 2.0) / DeltaCoefRad);
  int Indexs = floor(Vals);
  double Proporcaos = Vals - Indexs;
  double CDfxz = VetCDfxz[Indexa] + Proporcaoa * (VetCDfxz[Indexa + 1] - VetCDfxz[Indexa]);
  double CLfxz = VetCLf[Indexa] + Proporcaoa * (VetCLf[Indexa + 1] - VetCLf[Indexa]);
  double CMfxz = VetCMf[Indexa] + Proporcaoa * (VetCMf[Indexa + 1] - VetCMf[Indexa]);

  double CDfxy = VetCDfxy[Indexs] + Proporcaos * (VetCDfxy[Indexs + 1] - VetCDfxy[Indexs]);
  double CYfxy = VetCYf[Indexs] + Proporcaos * (VetCYf[Indexs + 1] - VetCYf[Indexs]);
  double CNfxy = VetCNf[Indexs] + Proporcaos * (VetCNf[Indexs + 1] - VetCNf[Indexs]);

  // Interpolating data from table | Computing the right wing coeficients
  Vala = ((Alphawr + pi / 2.0) / DeltaCoefRad);
  Indexa = floor(Vala);
  Proporcaoa = Vala - Indexa;

  Vals = ((Betawr + pi / 2.0) / DeltaCoefRad);
  Indexs = floor(Vals);
  Proporcaos = Vals - Indexs;

  double CDwRxz = VetCDwxz[Indexa] + Proporcaoa * (VetCDwxz[Indexa + 1] - VetCDwxz[Indexa]);
  double CLwRxz = VetCLw[Indexa] + Proporcaoa * (VetCLw[Indexa + 1] - VetCLw[Indexa]);
  double CMwRxz = VetCMw[Indexa] + Proporcaoa * (VetCMw[Indexa + 1] - VetCMw[Indexa]);

  double CDwRxy = VetCDwxy[Indexs] + Proporcaos * (VetCDwxy[Indexs + 1] - VetCDwxy[Indexs]);
  double CYwRxy = VetCYw[Indexs] + Proporcaos * (VetCYw[Indexs + 1] - VetCYw[Indexs]);
  double CNwRxy = VetCNw[Indexs] + Proporcaos * (VetCNw[Indexs + 1] - VetCNw[Indexs]);

  // Interpolating data from table | Computing the left wing coeficients
  Vala = ((Alphawl + pi / 2.0) / DeltaCoefRad);
  Indexa = floor(Vala);
  Proporcaoa = Vala - Indexa;

  Vals = ((Betawl + pi / 2.0) / DeltaCoefRad);
  Indexs = floor(Vals);
  Proporcaos = Vals - Indexs;

  double CDwLxz = VetCDwxz[Indexa] + Proporcaoa * (VetCDwxz[Indexa + 1] - VetCDwxz[Indexa]);
  double CLwLxz = VetCLw[Indexa] + Proporcaoa * (VetCLw[Indexa + 1] - VetCLw[Indexa]);
  double CMwLxz = VetCMw[Indexa] + Proporcaoa * (VetCMw[Indexa + 1] - VetCMw[Indexa]);

  double CDwLxy = VetCDwxy[Indexs] + Proporcaos * (VetCDwxy[Indexs + 1] - VetCDwxy[Indexs]);
  double CYwLxy = VetCYw[Indexs] + Proporcaos * (VetCYw[Indexs + 1] - VetCYw[Indexs]);
  double CNwLxy = VetCNw[Indexs] + Proporcaos * (VetCNw[Indexs + 1] - VetCNw[Indexs]);

  // Interpolating data from table | Computing the right tail coeficients
  Vala = ((Alphatr + pi / 2.0) / DeltaCoefRad);
  Indexa = floor(Vala);
  Proporcaoa = Vala - Indexa;

  Vals = ((Betatr + pi / 2.0) / DeltaCoefRad);
  Indexs = floor(Vals);
  Proporcaos = Vals - Indexs;

  double CDtRxz = VetCDtxz[Indexa] + Proporcaoa * (VetCDtxz[Indexa + 1] - VetCDtxz[Indexa]);
  double CLtRxz = VetCLt[Indexa] + Proporcaoa * (VetCLt[Indexa + 1] - VetCLt[Indexa]);
  double CMtRxz = VetCMt[Indexa] + Proporcaoa * (VetCMt[Indexa + 1] - VetCMt[Indexa]);

  double CDtRxy = VetCDtxy[Indexs] + Proporcaos * (VetCDtxy[Indexs + 1] - VetCDtxy[Indexs]);
  double CYtRxy = VetCYt[Indexs] + Proporcaos * (VetCYt[Indexs + 1] - VetCYt[Indexs]);
  double CNtRxy = VetCNt[Indexs] + Proporcaos * (VetCNt[Indexs + 1] - VetCNt[Indexs]);

  // Interpolating data from table | Computing the left tail coeficients
  Vala = ((Alphatl + pi / 2.0) / DeltaCoefRad);
  Indexa = floor(Vala);
  Proporcaoa = Vala - Indexa;

  Vals = ((Betatl + pi / 2.0) / DeltaCoefRad);
  Indexs = floor(Vals);
  Proporcaos = Vals - Indexs;

  double CDtLxz = VetCDtxz[Indexa] + Proporcaoa * (VetCDtxz[Indexa + 1] - VetCDtxz[Indexa]);
  double CLtLxz = VetCLt[Indexa] + Proporcaoa * (VetCLt[Indexa + 1] - VetCLt[Indexa]);
  double CMtLxz = VetCMt[Indexa] + Proporcaoa * (VetCMt[Indexa + 1] - VetCMt[Indexa]);

  double CDtLxy = VetCDtxy[Indexs] + Proporcaos * (VetCDtxy[Indexs + 1] - VetCDtxy[Indexs]);
  double CYtLxy = VetCYt[Indexs] + Proporcaos * (VetCYt[Indexs + 1] - VetCYt[Indexs]);
  double CNtLxy = VetCNt[Indexs] + Proporcaos * (VetCNt[Indexs + 1] - VetCNt[Indexs]);

  // Interpolating data from table | Computing the interference coeficients

  Vala = ((AlphaCG + pi / 2.0) / DeltaCoefRad);
  Indexa = floor(Vala);
  Proporcaoa = Vala - Indexa;

  Vals = ((BetaCG + pi / 2.0) / DeltaCoefRad);
  Indexs = floor(Vals);
  Proporcaos = Vals - Indexs;

  double CDixz = VetCDixz[Indexa] + Proporcaoa * (VetCDixz[Indexa + 1] - VetCDixz[Indexa]);
  double CLixz = VetCLi[Indexa] + Proporcaoa * (VetCLi[Indexa + 1] - VetCLi[Indexa]);
  double CMixz = VetCMi[Indexa] + Proporcaoa * (VetCMi[Indexa + 1] - VetCMi[Indexa]);

  double CDixy = VetCDixy[Indexs] + Proporcaos * (VetCDixy[Indexs + 1] - VetCDixy[Indexs]);
  double CYixy = VetCYi[Indexs] + Proporcaos * (VetCYi[Indexs + 1] - VetCYi[Indexs]);
  double CNixy = VetCNi[Indexs] + Proporcaos * (VetCNi[Indexs + 1] - VetCNi[Indexs]);

  double Cla = 0.7104;
  double Cle = 0.3759;
  double Cyr = -0.0730;

  Alphaf = -Alphaf;
  Betaf = -Betaf;
  Alphawr = -Alphawr;
  Betawr = -Betawr;
  Alphawl = -Alphawl;
  Betawl = -Betawl;
  Alphatr = -Alphatr;
  Betatr = -Betatr;
  Alphatl = -Alphatl;
  Betatl = -Betatl;
  AlphaCG = -AlphaCG;
  BetaCG = -BetaCG;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Rotation matrix Wind frame to aerodynamic center
  RBFw = Roty(-Alphaf) * Rotz(-Betaf);
  RBwRw = Roty(-Alphawr) * Rotz(-Betawr);
  RBwLw = Roty(-Alphawl) * Rotz(-Betawl);
  RBtRw = Roty(-Alphatr) * Rotz(-Betatr);
  RBtLw = Roty(-Alphatl) * Rotz(-Betatl);
  RBcGw = Roty(-AlphaCG) * Rotz(-BetaCG);
  // std::cout << "Alphaf " << Alphaf << "Betaf" << Betaf << std::endl;
  // std::cout << RBFw << std::endl << std::endl;
  // std::cout << RBwRw << std::endl << std::endl;
  // std::cout << RBwLw << std::endl << std::endl;
  // std::cout << RBtRw << std::endl << std::endl;
  // std::cout << RBtLw << std::endl << std::endl;
  // std::cout << RBcGw << std::endl << std::endl;

  // dynamic pressure
  double kf = 0.5 * Ho * sw * pow(Vf, 2);
  double kwR = 0.5 * Ho * sw * pow(Vwr, 2);
  double kwL = 0.5 * Ho * sw * pow(Vwl, 2);
  double ktR = 0.5 * Ho * sw * pow(Vtr, 2);
  double ktL = 0.5 * Ho * sw * pow(Vtl, 2);
  double kcG = 0.5 * Ho * sw * pow(VCG, 2);
  // std::cout << " kf " << kf  << " kwR " << kwR << " kwL " << kwL << " ktR " << ktR << " ktL " << ktL << " kcG " <<
  // kcG << std::endl;

  Eigen::VectorXd FFuselage(3), FWingR(3), FWingL(3), FVTailR(3), FVTailL(3), FInterference(3), Fuselage(3), WingR(3),
      WingL(3), VTailR(3), VTailL(3), AerodynamicMoments(3), Interference(3);

  // std::cout << " CDfxy " << CDfxy  << " CDfxz " << CDfxz << " CYfxy " << CYfxy << " CLfxz " << CLfxz << std::endl;

  FFuselage << -kf * (CDfxy + CDfxz), kf * CYfxy, kf * CLfxz;
  Fuselage = RBFw * FFuselage;

  // std::cout << " CDwRxy " << CDwRxy  << " CDwRxz " << CDwRxz << " CYwRxy " << CYwRxy << " CLwRxz " << CLwRxz <<
  // std::endl;

  FWingR << -kwR * (CDwRxy + CDwRxz), kwR * CYwRxy, kwR * CLwRxz;
  WingR = RBwRw * FWingR;

  FWingL << -kwL * (CDwLxy + CDwLxz), kwL * CYwLxy, kwL * CLwLxz;
  WingL = RBwLw * FWingL;

  FVTailR << -ktR * (CDtRxy + CDtRxz), ktR * CYtRxy, ktR * CLtRxz;
  VTailR = RBtRw * FVTailR;

  FVTailL << -ktL * (CDtLxy + CDtLxz), ktL * CYtLxy, ktL * CLtLxz;
  VTailL = RBtLw * FVTailL;

  AerodynamicMoments << 0,
      -(kf * CMfxz + kwR * CMwRxz + kwL * CMwLxz + ktR * CMtRxz + ktL * CMtLxz) * MeanGeometricChord,
      -(kf * CNfxy + kwR * CNwRxy + kwL * CNwLxy + ktR * CNtRxy + ktL * CNtLxy) * WingSpan;

  FInterference << -kcG * (CDixy + CDixz), kcG * CYixy, kcG * CLixz;
  Interference = RBcGw * FInterference;

  FUaFuselage = Jvf.transpose() * RI_B * Fuselage;
  // std::cout << Jvf << std::endl << std::endl;
  // std::cout << RI_B << std::endl << std::endl;
  // std::cout << RBFw << std::endl << std::endl;
  // std::cout << Fuselage << std::endl << std::endl;

  FUaWings = JvwR.transpose() * RI_B * WingR + JvwL.transpose() * RI_B * WingL;
  FUaTail = JvtR.transpose() * RI_B * VTailR + JvtL.transpose() * RI_B * VTailL;
  FUaInterference = Jvcg.transpose() * RI_B * Interference;
  FUaAeroMoments = Jw.transpose() * RI_B * AerodynamicMoments;
  FUa = FUaFuselage + FUaWings + FUaTail + FUaInterference + FUaAeroMoments;
  BAero << JvwR.transpose() * RI_B * RBwRw * kwR * az * Cla, JvwL.transpose() * RI_B * RBwLw * kwL * az * Cla,
      (JvtR.transpose() * RI_B * RBtRw * ktR + JvtL.transpose() * RI_B * RBtLw * ktL) * az * Cle,
      (JvtR.transpose() * RI_B * RBtRw * ktR + JvtL.transpose() * RI_B * RBtLw * ktL) * ay * Cyr;

  BAeroFUa << BAero(2, 0), BAero(2, 1), BAero(2, 2), BAero(2, 3), FUa(2), BAero(3, 0), BAero(3, 1), BAero(3, 2),
      BAero(3, 3), FUa(3), BAero(0, 0), BAero(0, 1), BAero(0, 2), BAero(0, 3), FUa(0), BAero(1, 0), BAero(1, 1),
      BAero(1, 2), BAero(1, 3), FUa(1), BAero(4, 0), BAero(4, 1), BAero(4, 2), BAero(4, 3), FUa(4), BAero(5, 0),
      BAero(5, 1), BAero(5, 2), BAero(5, 3), FUa(5), BAero(6, 0), BAero(6, 1), BAero(6, 2), BAero(6, 3), FUa(6),
      BAero(7, 0), BAero(7, 1), BAero(7, 2), BAero(7, 3), FUa(7);
  return BAeroFUa;
}
