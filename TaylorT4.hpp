/*
 * TaylorT4.hpp
 *
 *  Created on: Oct 23, 2018
 *      Author: blakemoore
 */

#ifndef TAYLORT4_HPP_
#define TAYLORT4_HPP_

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <fftw3.h>
#include <boost/array.hpp>
#include <boost/ref.hpp>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include <boost/numeric/odeint/config.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4_dense_output.hpp>
#include <functional>
#include "Amplitudes.hpp"

using namespace std;
using namespace boost::numeric::odeint;
namespace pl = std::placeholders;
typedef boost::array< double , 4 > state_type;


class TaylorT4 {
public:
	TaylorT4();
	TaylorT4(double p_c, double e_c, double M_c, double eta_c, double DL_c, double thet_c, double psi_c, double phi_c, double iota_c, double beta_c, int sr, int orb_pn_set);
	virtual ~TaylorT4();

	double edot(double y, double e);
	double ydot(double y, double e);
	double ldot(double y, double e);
	double lamdot(double y, double e);
	double t_to_coal();
	void rhs(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_1pn(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_2pn(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_7(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_8(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_6(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_10(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_12(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_11(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_exact_ratio(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_l_6(const state_type &x,  state_type &dxdt,  const double t);
	void rhs_e_l_8(const state_type &x,  state_type &dxdt,  const double t);
	void solorb();
	void sol_orb_e_param_6();
	void sol_orb_e_param_10();
	void sol_orb_e_param_8();
	void sol_orb_e_param_11();
	void sol_orb_e_param_12();
	void sol_orb_e_exact_ratio();
	void sol_orb_e_l_6();
	void sol_orb_e_l_8();
	double get_u3PN(double l, double e, double y);
	double f4t(double, double);
	double f6t(double, double);
	double g4t(double, double);
	double g6t(double, double);
	double i6t(double, double);
	double h6t(double, double);
	double f4p(double, double);
	double f6p(double, double);
	double g4p(double, double);
	double g6p(double, double);
	double i6p(double, double);
	double h6p(double, double);
	double get_v(double y, double e, double u);
	double get_W(double l, double e, double y, double u);
	void fill_all();
	vector<double> get_esol();
	vector<double> get_ysol();
	vector<double> get_lsol();
	vector<double> get_lamsol();
	vector<double> get_tsol();
	vector<double> get_hsol();
	vector<double> get_rsol();
	vector<double> get_phisol();
	double P_0_C_2_C_2(double, double);
	double P_0_C_2_S_2(double, double);
	double P_0(double, double);
	double X_0_C_2_C_2(double, double);
	double X_0_C_2_S_2(double, double);
	double H_P_0(double e, double u, double c2w, double s2w, double c2lam, double s2lam);
	double H_X_0(double e, double u, double c2w, double s2w, double c2lam, double s2lam);
	double F_p(double thet, double phi, double psi);
	double F_c(double thet, double phi, double psi);
	double h_0(double x, double e, double u, double w, double lam);
	double h_05(double x, double e, double u, double w, double lam);
	double h_1(double x, double e, double u, double w, double lam);
	double h_15(double x, double e, double u, double w, double lam);
	double h_2(double x, double e, double u, double w, double lam);
	void sol_h();
	vector<double> pad_h();
	void sol_h_e();
	vector<double> pad_h_e();
	void make_T4();
	void make_T4_e();
	vector<vector<double> > get_hfsol_e();
	void interp_map_solns_e();
	vector<vector<double> > get_hfsol();
	double H_harmdecomp(double y, double e, double l, double lam);
	vector<vector<double> > make_T4_decomp();
	vector<vector<double> > make_T4_decomp(int);
	void check_keplers();
	double H_harmdecomp_j(double y, double e, double l, double lam, int j);
	double amplookup_j(double y, double e, double eta, int j);
	vector<double> get_fnsol();
	vector<double> get_fwsol();
	double get_df();
	double get_F_lambda_end();
	double get_F_lambda_begin();
	double get_F_l_end();
	double get_F_l_begin();
	double H_harmdecomp_n(double y, double e, double l, int n);
	double amplookup_n(double e, int n);
	vector<vector<double> > make_T4_decomp_n(int n);
	double get_p0();
	double get_e0();
	double get_M();
	double get_eta();
	double get_DL();
	double get_theta();
	double get_psi();
	double get_phi();
	double get_beta();
	double get_iota();
	double get_last_f();
	void sol_orb_no_RR(int ncyc);
	void sol_h(int ncyc);
	void sol_h_pn(int pn);
	vector<double> pad_h(int ncyc);
	vector<double> pad_h_pn(int pn);
	void make_T4(int ncyc);
	void make_T4_pn(int pn);
	vector<vector<double> > get_hfsol(int ncyc);
	vector<vector<double> > get_hfsol_pn(int pn);
	vector<double> pad_h_hann();
	void sol_orb_e_param();
	double dyde(double, double);
	double dtde(double, double);
	double dlde(double, double);
	double dlamde(double, double);
	double dyde_6(double, double);
	double dtde_6(double, double);
	double dlde_6(double, double);
	double dlamde_6(double, double);
	vector<vector<double> > get_hfsol_downsamp(double fend, double df_lim);
	double get_eref();
	vector<vector<double> > get_hfsol_downsamp_e(double fend, double df_lim);

//////////////////////////////////////////////////////
// Member Variables
//////////////////////////////////////////////////////
private:
	double p0;
	double e0;
	double M;
	double eta;
	double DL;
	double thet;
	double psi;
	double phi;
	double iota;
	double beta;
	int orb_pn;
	int sr;
	vector<double> lsol;
	vector<double> lamsol;
	vector<double> esol;
	vector<double> ysol;
	vector<double> tsol;
	vector<double> usol;
	vector<double> wsol;
	vector<double> hsol;
	vector<vector<double> > h_T4_f;
};

#endif
