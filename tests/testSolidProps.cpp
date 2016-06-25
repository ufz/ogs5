//**
// * \file testSolidProps.cpp
// * 2011-01-31 TN Initial implementation
// *
// * Tests in the context of the Burgers model implementation
// */

// ** INCLUDES **
#include "gtest.h"

#include "fem_ele_vec.h"
#include "rf_msp_new.h"
#include "Eigen/Dense"
#include "makros.h"
#include <cfloat>
#include "burgers.h"
#include "minkley.h"
#include "invariants.h"
//#include "stdio.h"
#include <fstream>

TEST(SolidProps, EffectiveStress)
{
	SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;

		
	//stress vectors
	Eigen::Matrix<double,6,1> sig, sigd;
	sigd.setZero(6); sig.setZero(6);

    //stress values
    double *stress;
    stress = new double[6];
    stress[0] = 0.01; stress[1] = 0.02; stress[2] = 0.03;
    stress[3] = 0.04; stress[4] = 0.05; stress[5] = 0.06;
    sig  = solid.Voigt_to_Kelvin_Stress(stress);//Kelvin mapping
    delete [] stress;
    sigd = smath.P_dev*sig;

    const double s_eff(smath.CalEffectiveStress(sig));
    const double sd_eff(smath.CalEffectiveStress(sigd));

	ASSERT_NEAR(s_eff,0.158745078664,1.e-10);
	ASSERT_NEAR(sd_eff,0.152970585408,1.e-10);
}

TEST(SolidProps, LodeAngle)
{
    SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;

    //stress vectors
    Eigen::Matrix<double,6,1> sig, sigd;
    sigd.setZero(6); sig.setZero(6);

    //uniaxial tension
    sig(0) = 2.;
    sigd = smath.P_dev*sig;

    double theta(smath.CalLodeAngle(sigd));
    ASSERT_NEAR(theta,-30.*PI/180.,1.e-8);

    //uniaxial compression
    sig(0) = -2.;
    sigd = smath.P_dev*sig;

    theta = smath.CalLodeAngle(sigd);
    ASSERT_NEAR(theta,30.*PI/180.,1.e-8);

    //equibiaxial tension
    sig(0) = 2.;
    sig(2) = 2.;
    sigd = smath.P_dev*sig;

    theta = smath.CalLodeAngle(sigd);
    ASSERT_NEAR(theta,30.*PI/180.,1.e-8);

    //biaxial tension-compression
    sig(0) = -2.;
    sig(2) = 2.;
    sigd = smath.P_dev*sig;

    theta = smath.CalLodeAngle(sigd);
    ASSERT_NEAR(theta,0.,1.e-8);

    //hydrostatic loading
    sig(0) = 2.;
    sig(1) = 2.;
    sig(2) = 2.;
    sigd = smath.P_dev*sig;

    theta = smath.CalLodeAngle(sigd);
    ASSERT_NEAR(theta,0.,1.e-8);

}

//TEST(SolidProps, Lubby2Residual)
//{
//    SolidProp::CSolidProperties solid;
//    SolidMath::Invariants smath;
//    Math_Group::Matrix* data;
//    data = new Math_Group::Matrix(13);
//
//    //set Constants
//    (*data)(0) = 2.0; //Kelvin shear modulus
//    (*data)(1) = 10.0; //dependency parameter for "
//    (*data)(2) = 2.0; //Kelvin viscosity
//    (*data)(3) = 10.0; //dependency parameter for "
//    (*data)(4) = 1.7; //Maxwell shear modulus
//    (*data)(5) = 2.0; //Maxwell bulk modulus
//    (*data)(6) = 10.0; //Maxwell viscosity
//    (*data)(7) = 10.0; //dependency parameter for "
//    (*data)(8) = 0.; // slope of elesticity temperature dependence
//    (*data)(9) = 0.; // slope of elesticity temperature dependence
//    (*data)(10) = 273.; // reference temperature dependency parameter for "
//    (*data)(11) = 1.; // constant factor for Arrhenius term
//    (*data)(12) = 0.; // activation energy in Arrhenius term
//
//    Burgers::SolidBurgers material(data);
//
//    //state and trial variables
//    Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j;
//    Eigen::Matrix<double,6,1> sigd_t, sigd_j, epsd_i, epsd_t;
//    eps_i.setZero(6);
//    eps_K_t.setZero(6); eps_K_j.setZero(6);
//    sigd_t.setZero(6); sigd_j.setZero(6);
//    epsd_t.setZero(6); epsd_i.setZero(6);
//
//    //Residual vector
//    Eigen::Matrix<double,12,1> residual;
//
//    //set a strain increment and perform Kelvin mapping
//    double *strain_bc;
//    strain_bc = new double[6];
//    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
//    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
//    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
//
//    delete [] strain_bc;
//
//    epsd_i = smath.P_dev*eps_i;
//
//    //guess stress increment (dimensionless)
//    sigd_j = 2.0 * epsd_i;
//
//    //Update Material parameters
//    material.UpdateBurgersProperties(smath.CalEffectiveStress(sigd_j* material.GM0),273.);
//
//    //Calculate residual for time step
//    const double dt(0.01);
//    material.CalResidualBurgers(epsd_i,epsd_t,sigd_j,sigd_t,eps_K_j,eps_K_t,residual);
//
//    Eigen::Matrix<double,12,1> residual_Python;
//
//    //These values are taken from a constituive model test script coded in Python
//    residual_Python(0) = -1.1242544333e-04; residual_Python(1) = 3.9005411168e-20;
//    residual_Python(2) = 1.1242544333e-04; residual_Python(3) = 6.3597434684e-04;
//    residual_Python(4) = 7.9496793355e-04; residual_Python(5) = 9.5396152026e-04;
//    residual_Python(6) = 4.6843934720e-05; residual_Python(7) =  -1.6252254653e-20;
//    residual_Python(8) = -4.6843934720e-05 ; residual_Python(9) = -2.6498931118e-04;
//    residual_Python(10) = -3.3123663898e-04; residual_Python(11) = -3.9748396677e-04;
//
//    for (size_t i=0; i<12; i++)
//            ASSERT_NEAR(residual(i),residual_Python(i),1.e-10);
//}

TEST(SolidProps, Lubby2JacobianNumeric)
{
    SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;
    Math_Group::Matrix* data;
    data = new Math_Group::Matrix(13);

    //set Constants
    (*data)(0) = 2.0; //Kelvin shear modulus
    (*data)(1) = 10.0; //dependency parameter for "
    (*data)(2) = 2.0; //Kelvin viscosity
    (*data)(3) = 10.0; //dependency parameter for "
    (*data)(4) = 1.7; //Maxwell shear modulus
    (*data)(5) = 2.0; //Maxwell bulk modulus
    (*data)(6) = 10.0; //Maxwell viscosity
    (*data)(7) = 10.0; //dependency parameter for "
    (*data)(8) = 0.; // slope of elesticity temperature dependence
    (*data)(9) = 0.; // slope of elesticity temperature dependence
    (*data)(10) = 273.; // reference temperature dependency parameter for "
    (*data)(11) = 1.; // constant factor for Arrhenius term
    (*data)(12) = 0.; // activation energy in Arrhenius term
    Burgers::SolidBurgers material(data);
	
	//state and trial variables
	Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j, eps_M_t, eps_M_j;
	Eigen::Matrix<double,6,1> sigd_t, sigd_j, epsd_i, epsd_t;
	eps_i.setZero(6);
	eps_K_t.setZero(6); eps_K_j.setZero(6);
	eps_M_t.setZero(6); eps_M_j.setZero(6);
	sigd_t.setZero(6); sigd_j.setZero(6);
	epsd_t.setZero(6); epsd_i.setZero(6);

    //set a strain increment and perform Kelvin mapping
    double *strain_bc;
    strain_bc = new double[6];
    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
    delete [] strain_bc;
	
    epsd_i = smath.P_dev*eps_i;

    //guess stress increment (dimensionless)
    sigd_j = 2.0 * epsd_i;

    //Update Material parameters
    material.UpdateBurgersProperties(smath.CalEffectiveStress(sigd_j* material.GM0),273.);

	//set nontrivial internal variables
	eps_K_j = 0.1 * epsd_i;
	eps_M_j = 0.15 * epsd_i;

    Eigen::Matrix<double,18,18> Jacobian;
    Jacobian.setZero(18,18);

    //Calculate Jacobian for time step
    const double dt(0.01);
    material.CalJacobianBurgers(dt,Jacobian,smath.CalEffectiveStress(sigd_j),sigd_j,eps_K_j);

    //Residual vector
    Eigen::Matrix<double,18,1> residual;

    //Numerically calculate Jacobian
    Eigen::Matrix<double,18,18> Jac_num;
    Jac_num.setZero(18,18);

    Eigen::Matrix<double,6,1> upper, lower;
    const double pertub(1.e-6);
    double up, low;

    for (size_t i=0; i<18; i++)
        for (size_t j=0; j<18; j++)
        {
            if (j < 6)
            {
                upper = sigd_j;
                lower = sigd_j;
                upper(j) += pertub;
                lower(j) -= pertub;
                material.UpdateBurgersProperties(smath.CalEffectiveStress(upper* material.GM0),273.);
                material.CalResidualBurgers(dt,epsd_i,upper,eps_K_j,eps_K_t,eps_M_j,eps_M_t,residual);
                up = residual(i);
                material.UpdateBurgersProperties(smath.CalEffectiveStress(lower* material.GM0),273.);
                material.CalResidualBurgers(dt,epsd_i,lower,eps_K_j,eps_K_t,eps_M_j,eps_M_t,residual);
                low = residual(i);
            }
            else if (j < 12)
            {
                upper = eps_K_j;
                lower = eps_K_j;
                upper(j-6) += pertub;
                lower(j-6) -= pertub;
                //Update Material parameters
                material.UpdateBurgersProperties(smath.CalEffectiveStress(sigd_j* material.GM0),273.);
                material.CalResidualBurgers(dt,epsd_i,sigd_j,upper,eps_K_t,eps_M_j,eps_M_t,residual);
                up = residual(i);
                material.CalResidualBurgers(dt,epsd_i,sigd_j,lower,eps_K_t,eps_M_j,eps_M_t,residual);
                low = residual(i);
            }
			else
            {
                upper = eps_M_j;
                lower = eps_M_j;
                upper(j-12) += pertub;
                lower(j-12) -= pertub;
                //Update Material parameters
                material.UpdateBurgersProperties(smath.CalEffectiveStress(sigd_j* material.GM0),273.);
                material.CalResidualBurgers(dt,epsd_i,sigd_j,eps_K_j,eps_K_t,upper,eps_M_t,residual);
                up = residual(i);
                material.CalResidualBurgers(dt,epsd_i,sigd_j,eps_K_j,eps_K_t,lower,eps_M_t,residual);
                low = residual(i);
            }
        Jac_num(i,j) = (up-low)/(2.*pertub);
        ASSERT_NEAR(Jacobian(i,j),Jac_num(i,j),1.e-8);
        }
}

TEST(SolidProps, YieldMohrCoulomb)
{
    Math_Group::Matrix* data;
    data = new Math_Group::Matrix(15);

    //set Constants
    for (int i(0); i<15; i++)
        (*data)(i) = (double)i;

    (*data)(7) = 2.; //cohesion
    (*data)(9) = 30.; //friction angle
    const double phi((*data)(9)*PI/180.);
    (*data)(11) = 29.9; //transition angle

    Minkley::SolidMinkley material(data);

    //state and trial variables
    Eigen::Matrix<double,6,1> sig_i;
    sig_i.setZero(6,1);

    //MC formulation
    sig_i(0) = -2.; //sigma_1
    sig_i(2) = sig_i(0) * (1.+std::sin(phi))/(1.-std::sin(phi)) - 2. * (*data)(7) * std::cos(phi)/(1.-std::sin(phi)); //sigma_3
    sig_i(1) = 0.5*(sig_i(0) + sig_i(2)); //sigma_2

    ASSERT_NEAR(0.,material.YieldMohrCoulomb(sig_i),1.e-12);

    //Corner smoothed region
    sig_i(0) = -2.; //sigma_1
    sig_i(2) = sig_i(0) * (1.+std::sin(phi))/(1.-std::sin(phi)) - 2. * (*data)(7) * std::cos(phi)/(1.-std::sin(phi)); //sigma_3
    sig_i(1) = sig_i(0); //sigma_2

    //Large tolerance wanted -- yield surfaces deviate in this region.
    ASSERT_NEAR(0.,material.YieldMohrCoulomb(sig_i),1.e-2);
}

//TEST(SolidProps, MinkleyCreepResidual)
//{
//    SolidProp::CSolidProperties solid;
//    SolidMath::Invariants smath;
//    Math_Group::Matrix* data;
//    data = new Math_Group::Matrix(15);

//    //set Constants
//    (*data)(0) = 63.e3; //Kelvin shear modulus
//    (*data)(1) = 14.e6; //Kelvin viscosity
//    (*data)(2) = 12.e3; //Maxwell shear modulus
//    (*data)(3) = 18.e3; //Maxwell bulk modulus
//    (*data)(4) = 10.e10; //Maxwell viscosity
//    (*data)(5) = 4.9; //dependency parameter for " (m)
//    (*data)(6) = 0.33; //dependency parameter for " (n)
//    (*data)(7) = 2.; //cohesion
//    (*data)(8) = 1.; //hardening
//    (*data)(9) = 30.; //friction angle
//    (*data)(10) = 10.; //dilatancy angle
//    (*data)(11) = 28.; //transition angle
//    (*data)(12) = 0.1; //regularisation
//    (*data)(13) = 0.; // temperature parameter for Maxwell viscosity
//    (*data)(14) = 0.; // reference temperature for Maxwell viscosity
//    Minkley::SolidMinkley material(data);

//    //state and trial variables
//    Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j;
//    Eigen::Matrix<double,6,1> sig_t, sigd_j, sig_j, epsd_i, epsd_t;
//    eps_i.setZero(6);
//    eps_K_t.setZero(6); eps_K_j.setZero(6);
//    sig_t.setZero(6); sigd_j.setZero(6); sig_j.setZero(6);
//    epsd_t.setZero(6); epsd_i.setZero(6);

//    //Residual vector
//    Eigen::Matrix<double,12,1> residual;

//    //set a strain increment and perform Kelvin mapping
//    double *strain_bc;
//    strain_bc = new double[6];
//    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
//    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
//    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
//    eps_i *= 0.01;

//    delete [] strain_bc;

//    epsd_i = smath.P_dev*eps_i;

//    const double e_i(smath.CalI1(eps_i)), e_t(0.);

//    //guess stress increment (dimensionless)
//    sigd_j = 2.0 * epsd_i;
//    sig_j = sigd_j + material.KM0/material.GM0 * e_i * smath.ivec;

//    //Update Material parameters
//    material.UpdateMinkleyProperties(smath.CalEffectiveStress(sigd_j * material.GM0),0.,0.);

//    //Calculate residual for time step
//    const double dt(1.5);

//    material.CalViscoelasticResidual(dt,epsd_i,epsd_t,e_i,e_t,sig_j,sig_t,eps_K_j,eps_K_t,residual);

//    Eigen::Matrix<double,12,1> residual_Python;

//    //These values are taken from a constituive model test script coded in Python
//    residual_Python(0) = -1.1680676716e-04; residual_Python(1) = 1.0842021725e-19;
//    residual_Python(2) = 1.1680676716e-04; residual_Python(3) = 6.6075885719e-04;
//    residual_Python(4) = 8.2594857149e-04; residual_Python(5) = 9.9113828579e-04;
//    residual_Python(6) = 8.5714285714e-08; residual_Python(7) =  -6.9698711088e-23;
//    residual_Python(8) = -8.5714285714e-08 ; residual_Python(9) = -4.8487322139e-07;
//    residual_Python(10) = -6.0609152673e-07; residual_Python(11) = -7.2730983208e-07;

//    for (size_t i=0; i<12; i++)
//        ASSERT_NEAR(residual(i),residual_Python(i),1.e-10);
//}

TEST(SolidProps, MinkleyJacobianNumeric)
{
    SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;
    Math_Group::Matrix* data;
    data = new Math_Group::Matrix(15);

    //set Constants
    (*data)(0) = 63.e3; //Kelvin shear modulus
    (*data)(1) = 14.e6; //Kelvin viscosity
    (*data)(2) = 12.e3; //Maxwell shear modulus
    (*data)(3) = 18.e3; //Maxwell bulk modulus
    (*data)(4) = 10.e10; //Maxwell viscosity
    (*data)(5) = 4.9; //dependency parameter for " (m)
    (*data)(6) = 0.33; //dependency parameter for " (n)
    (*data)(7) = 2.; //cohesion
    (*data)(8) = 1.; //hardening
    (*data)(9) = 30.; //friction angle
    (*data)(10) = 10.; //dilatancy angle
    (*data)(11) = 28.; //transition angle
    (*data)(12) = 0.1; //regularisation
    (*data)(13) = 0.; // temperature parameter for Maxwell viscosity
    (*data)(14) = 313.; // reference temperature for Maxwell viscosity
    Minkley::SolidMinkley material(data);

    //state and trial variables
    Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j, eps_M_t, eps_M_j, eps_pl_j, eps_pl_t;
    Eigen::Matrix<double,6,1> sig_t, sigd_j, sig_j, epsd_i, epsd_t;
    eps_i.setZero(6);
    eps_K_t.setZero(6); eps_K_j.setZero(6);
    eps_M_t.setZero(6); eps_M_j.setZero(6);
    eps_pl_t.setZero(6); eps_pl_j.setZero(6);
    sig_t.setZero(6); sigd_j.setZero(6); sig_j.setZero(6);
    epsd_t.setZero(6); epsd_i.setZero(6);

    //set a strain increment and perform Kelvin mapping
    double *strain_bc;
    strain_bc = new double[6];
    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
    eps_i *= 0.01;

    delete [] strain_bc;

    epsd_i = smath.P_dev*eps_i;

    //set nontrivial internal variables
    eps_K_j = 0.1 * epsd_i;
    eps_M_j = 0.15 * epsd_i;
    eps_pl_j = 0.1 * epsd_i;

    const double e_i(smath.CalI1(eps_i));
    const double e_pv_j(0.1*e_i);

    //guess stress increment (dimensionless)
    sigd_j = 2.0 * epsd_i;
    sig_j = sigd_j + material.KM0/material.GM0 * e_i * smath.ivec;

    //Update Material parameters
    material.UpdateMinkleyProperties(smath.CalEffectiveStress(sigd_j * material.GM0),0.,0.);

    Eigen::Matrix<double,18,18> Jacobian;
    Jacobian.setZero(18,18);

    //Calculate Jacobian for time step
    const double dt(1.5);
    material.CalViscoelasticJacobian(dt,sig_j,smath.CalEffectiveStress(sigd_j),Jacobian);

    //Residual vector
    Eigen::Matrix<double,18,1> residual;

    //Numerically calculate Jacobian
    Eigen::Matrix<double,18,18> Jac_num;
    Jac_num.setZero(18,18);

    Eigen::Matrix<double,6,1> upper, lower;
    const double pertub(1.e-8);
    double up, low;

    for (size_t i=0; i<18; i++)
        for (size_t j=0; j<18; j++)
        {
            if (j < 6)
            {
                upper = sig_j;
                lower = sig_j;
                upper(j) += pertub;
                lower(j) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*upper * material.GM0),0.,0.);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,upper,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,residual);
                up = residual(i);
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*lower * material.GM0),0.,0.);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,lower,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,residual);
                low = residual(i);
            }
            else if (j < 12)
            {
                upper = eps_K_j;
                lower = eps_K_j;
                upper(j-6) += pertub;
                lower(j-6) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),0.,0.);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,sig_j,upper,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,residual);
                up = residual(i);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,sig_j,lower,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,residual);
                low = residual(i);
            }
            else
            {
                upper = eps_M_j;
                lower = eps_M_j;
                upper(j-12) += pertub;
                lower(j-12) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),0.,0.);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,sig_j,eps_K_j,eps_K_t,upper,eps_M_t,eps_pl_j,residual);
                up = residual(i);
                material.CalViscoelasticResidual(dt,epsd_i,e_i,e_pv_j,sig_j,eps_K_j,eps_K_t,lower,eps_M_t,eps_pl_j,residual);
                low = residual(i);
            }
            Jac_num(i,j) = (up-low)/(2.*pertub);
            ASSERT_NEAR(Jacobian(i,j),Jac_num(i,j),1.e-8);
        }
}

TEST(SolidProps, TensorInversion)
{
    SolidMath::Invariants smath;

    //state and trial variables
    Eigen::Matrix<double,6,1> vec, inv, inv_t;
    vec(0) = 1.7; vec(1) = -1.3; vec(2) = 0.8;
    vec(3) = 1.1*std::sqrt(2.); vec(4) = -2.0*std::sqrt(2.); vec(5) = 0.3*std::sqrt(2.);

    //Python result
    inv(0) = 0.4693174411;
    inv(1) = -0.1182605457;
    inv(2) = 0.3184654065;
    inv(3) = 0.1378154391*std::sqrt(2.);
    inv(4) = -0.3473321538*std::sqrt(2.);
    inv(5) = 0.1685445572*std::sqrt(2.);

    //OGS result
    inv_t = smath.InvertVector(vec);

    for (size_t i=0; i<6; i++)
        ASSERT_NEAR(inv_t(i),inv(i),1.e-10);
}

TEST(SolidProps, MinkleyFullResidual)
{
    SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;
    Math_Group::Matrix* data;
    data = new Math_Group::Matrix(15);

    //set Constants
    (*data)(0) = 63.e3; //Kelvin shear modulus
    (*data)(1) = 8.37e4; //Kelvin viscosity
    (*data)(2) = 12.e3; //Maxwell shear modulus
    (*data)(3) = 18.e3; //Maxwell bulk modulus
    (*data)(4) = 4.03e7; //Maxwell viscosity
    (*data)(5) = 4.9; //dependency parameter for " (m)
    (*data)(6) = 0.33; //dependency parameter for " (n)
    (*data)(7) = 2.; //cohesion
    (*data)(8) = 100.; //hardening
    (*data)(9) = 30.; //friction angle
    (*data)(10) = 10.; //dilatancy angle
    (*data)(11) = 29.; //transition angle
    (*data)(12) = 0.01; //regularisation
    (*data)(13) = 0.; // temperature parameter for Maxwell viscosity
    (*data)(14) = 0.; // reference temperature for Maxwell viscosity
    Minkley::SolidMinkley material(data);

    //state and trial variables
    Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j, eps_M_t, eps_M_j, eps_pl_t, eps_pl_j;
    Eigen::Matrix<double,6,1> sig_t, sigd_j, sig_j, epsd_i, epsd_t;
    eps_i.setZero(6);
    eps_K_t.setZero(6); eps_K_j.setZero(6);
    eps_M_t.setZero(6); eps_M_j.setZero(6);
    eps_pl_t.setZero(6); eps_pl_j.setZero(6);
    sig_t.setZero(6); sigd_j.setZero(6); sig_j.setZero(6);
    epsd_t.setZero(6); epsd_i.setZero(6);

    //Residual vector
    Eigen::Matrix<double,27,1> residual;

    //set a strain increment and perform Kelvin mapping
    double *strain_bc;
    strain_bc = new double[6];
    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
    eps_i *= -0.01;

    delete [] strain_bc;

    epsd_i = smath.P_dev*eps_i;

    const double e_i(smath.CalI1(eps_i));

    //guess stress increment (dimensionless)
    sigd_j = 2.0 * epsd_i;
    sig_j = sigd_j + material.KM0/material.GM0 * e_i * smath.ivec;

    //Calculate residual for time step
    const double dt(10.);
    //Modify plastic arrays to ensure full check
    //do get nonzero flow driven residuals
    const double lam_j(0.01);
    eps_K_j = 0.1 * epsd_i;
    eps_M_j = 0.15 * epsd_i;
    eps_pl_j = 0.2 * epsd_i;
    const double e_v_i = 0.2*e_i;
    const double e_eff_i = 1.2*e_i;

    //Update Material parameters
    material.UpdateMinkleyProperties(smath.CalEffectiveStress(sigd_j * material.GM0),e_eff_i,0.);

    material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);

    Eigen::Matrix<double,27,1> residual_Python;

    //These values are taken from a constituive model test script coded in Python
    residual_Python( 0 ) =  -9e-05 ;
    residual_Python( 1 ) =  -0.00018 ;
    residual_Python( 2 ) =  -0.00027 ;
    residual_Python( 3 ) =  -0.000509116882454 ;
    residual_Python( 4 ) =  -0.000636396103068 ;
    residual_Python( 5 ) =  -0.000763675323681 ;
    residual_Python( 6 ) =  -5.81003584229e-06 ;
    residual_Python( 7 ) =  8.19126399401e-21 ;
    residual_Python( 8 ) =  5.81003584229e-06 ;
    residual_Python( 9 ) =  3.28665259442e-05 ;
    residual_Python( 10 ) =  4.10831574302e-05 ;
    residual_Python( 11 ) =  4.92997889163e-05 ;
    residual_Python( 12 ) =  -0.144920047348 ;
    residual_Python( 13 ) =  1.17842582492e-16 ;
    residual_Python( 14 ) =  0.144920047348 ;
    residual_Python( 15 ) =  0.819791585675 ;
    residual_Python( 16 ) =  1.02473948209 ;
    residual_Python( 17 ) =  1.22968737851 ;
    residual_Python( 18 ) =  -0.00140743973352 ;
    residual_Python( 19 ) =  0.00173295588137 ;
    residual_Python( 20 ) =  -0.000325516147857 ;
    residual_Python( 21 ) =  0.00127776404536 ;
    residual_Python( 22 ) =  0.00237699541911 ;
    residual_Python( 23 ) =  0.00614979374969 ;
    residual_Python( 24 ) =  -0.00174848177667 ;
    residual_Python( 25 ) =  -0.00587444554726 ;
    residual_Python( 26 ) =  0.000736664715603 ;

    for (size_t i=0; i<27; i++){
        ASSERT_NEAR(residual(i),residual_Python(i),1.e-10);
    }
}

TEST(SolidProps, MinkleyFullJacobian)
{
    SolidProp::CSolidProperties solid;
    SolidMath::Invariants smath;
    Math_Group::Matrix* data;
    data = new Math_Group::Matrix(15);

    //set Constants
    (*data)(0) = 63.e3; //Kelvin shear modulus
    (*data)(1) = 8.37e4; //Kelvin viscosity
    (*data)(2) = 12.e3; //Maxwell shear modulus
    (*data)(3) = 18.e3; //Maxwell bulk modulus
	(*data)(4) = 4.03e7; //Maxwell viscosity
    (*data)(5) = 4.9; //dependency parameter for " (m)
    (*data)(6) = 0.33; //dependency parameter for " (n)
    (*data)(7) = 2.; //cohesion
    (*data)(8) = 100.; //hardening
    (*data)(9) = 30.; //friction angle
    (*data)(10) = 10.; //dilatancy angle
    (*data)(11) = 29.; //transition angle
    (*data)(12) = 0.01; //regularisation
    (*data)(13) = 0.; // temperature parameter for Maxwell viscosity
    (*data)(14) = 0.; // reference temperature for Maxwell viscosity
    Minkley::SolidMinkley material(data);

    //state and trial variables
    Eigen::Matrix<double,6,1> eps_i, eps_K_t, eps_K_j, eps_M_t, eps_M_j, eps_pl_t, eps_pl_j;
    Eigen::Matrix<double,6,1> sig_t, sigd_j, sig_j, epsd_i, epsd_t;
    eps_i.setZero(6);
    eps_K_t.setZero(6); eps_K_j.setZero(6);
    eps_M_t.setZero(6); eps_M_j.setZero(6);
    eps_pl_t.setZero(6); eps_pl_j.setZero(6);
    sig_t.setZero(6); sigd_j.setZero(6); sig_j.setZero(6);
    epsd_t.setZero(6); epsd_i.setZero(6);

    //Residual vector
    Eigen::Matrix<double,27,1> residual;

    //set a strain increment and perform Kelvin mapping
    double *strain_bc;
    strain_bc = new double[6];
    strain_bc[0] = 0.01; strain_bc[1] = 0.02; strain_bc[2] = 0.03;
    strain_bc[3] = 0.04*2.; strain_bc[4] = 0.05*2.; strain_bc[5] = 0.06*2.;
    eps_i = solid.Voigt_to_Kelvin_Strain(strain_bc);
    eps_i *= -0.01;

    delete [] strain_bc;

    epsd_i = smath.P_dev*eps_i;

    const double e_i(smath.CalI1(eps_i));

    //guess stress increment (dimensionless)
	sigd_j = 2.0 * epsd_i;
	sig_j = sigd_j + material.KM0/material.GM0 * e_i * smath.ivec;

    //Calculate residual for time step
    const double dt(10.);
    //Modify plastic arrays to ensure full check
    //do get nonzero flow driven residuals
    const double lam_j(0.01);
    eps_K_j = 0.1 * epsd_i;
    eps_M_j = 0.15 * epsd_i;
    eps_pl_j = 0.2 * epsd_i;
    const double e_v_i = 0.2*e_i;
    const double e_eff_i = 1.2*e_i;

    //Update Material parameters
    material.UpdateMinkleyProperties(smath.CalEffectiveStress(sigd_j * material.GM0),e_eff_i,0.);

    //Calculate Jacobian for time step
    Eigen::Matrix<double,27,27> Jacobian;
    material.CalViscoplasticJacobian(dt,sig_j,smath.CalEffectiveStress(sigd_j),lam_j,Jacobian);

    //Numerically calculate Jacobian
    Eigen::Matrix<double,27,27> Jac_num;
    Jac_num.setZero(27,27);

    Eigen::Matrix<double,6,1> upper, lower;
    const double pertub(1.e-8);
    double up, low;

    for (size_t i=0; i<27; i++)
        for (size_t j=0; j<27; j++)
        {
            if (j < 6)
            {
                upper = sig_j;
                lower = sig_j;
                upper(j) += pertub;
                lower(j) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*upper * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,upper,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                up = residual(i);
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*lower * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,lower,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 12)
            {
                upper = eps_K_j;
                lower = eps_K_j;
                upper(j-6) += pertub;
                lower(j-6) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,upper,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                up = residual(i);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,lower,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 18)
            {
                upper = eps_M_j;
                lower = eps_M_j;
                upper(j-12) += pertub;
                lower(j-12) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,upper,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                up = residual(i);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,lower,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 24)
            {
                upper = eps_pl_j;
                lower = eps_pl_j;
                upper(j-18) += pertub;
                lower(j-18) -= pertub;
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,upper,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                up = residual(i);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,lower,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 25)
            {
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i+pertub,0.,e_eff_i,0.,lam_j,residual);
                up = residual(i);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i-pertub,0.,e_eff_i,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 26)
            {
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i+pertub,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i+pertub,0.,lam_j,residual);
                up = residual(i);
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i-pertub,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i-pertub,0.,lam_j,residual);
                low = residual(i);
            }
            else if (j < 27)
            {
                //Update Material parameters
                material.UpdateMinkleyProperties(smath.CalEffectiveStress(smath.P_dev*sig_j * material.GM0),e_eff_i,0.);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j+pertub,residual);
                up = residual(i);
                material.CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t,e_v_i,0.,e_eff_i,0.,lam_j-pertub,residual);
                low = residual(i);
            }
            Jac_num(i,j) = (up-low)/(2.*pertub);
			//std::cout << i << " " << j << " "  << Jac_num(i,j) << " "  << Jacobian(i,j) << std::endl;
            ASSERT_NEAR(Jacobian(i,j),Jac_num(i,j),1.e-6);
        }
}
