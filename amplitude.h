#ifndef BREIT_WIGNERS_WEIT_BRIGNERS
#define BREIT_WIGNERS_WEIT_BRIGNERS
#include<string>
#include<complex>
#include<iostream>
#include<limits>
#include"amplitude_functions.h"

//// Definitions of Breit wigner functions

////////	////////	constant_function 		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef amplitude constant_function; // Amplitude base class is already the constant function


////////	////////	breit_wigner 			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class breit_wigner : public amplitude{

	public:
		breit_wigner();

		std::string type()											const		{return "breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

breit_wigner::breit_wigner():amplitude(1,2,0,0){

	_name = "unnamed_breit_wigner";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "width";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_width";
};

template <typename xdouble>
std::complex<xdouble> breit_wigner::template_eval(const double* var, const xdouble* par, const xdouble* con)		const{
	std::complex<xdouble> denominator = std::complex<xdouble>(par[0]*par[0]-var[0]*var[0],-par[0]*par[1]);
	return std::complex<xdouble>(par[0]*par[1])/denominator;
};

////////	////////	mass_dep_breit_wigner		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class mass_dep_breit_wigner : public amplitude{

	public:
		mass_dep_breit_wigner();

		std::string type()											const		{return "mass_dep_breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

mass_dep_breit_wigner::mass_dep_breit_wigner():amplitude(1,2,3,1){

	_name = "unnamed_mass_dep_breit_wigner";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "width";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_width";

	_con_types[0] = "m_Pi";
	_con_types[1] = "m_Iso";
	_con_types[2] = "angular_momentum";

	_con_names[0] = "m_Pi";
	_con_names[1] = "m_Iso";
	_con_names[2] = "L";

	
};

template <typename xdouble>
std::complex<xdouble> mass_dep_breit_wigner::template_eval(const double* var, const xdouble* par, const xdouble* con)	const{
	double m   = var[0];

	xdouble m0  = par[0];
	xdouble G0  = par[1];

	xdouble mPi = con[0];
	xdouble mIso= con[1];

	xdouble q0 = breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso*mIso);
	xdouble q  = breakupMomentumReal<xdouble>(m* m ,mPi*mPi,mIso*mIso);
	xdouble Fl = barrierFactor<xdouble>(q,con[2]);
	xdouble Fl0= barrierFactor<xdouble>(q0,con[2]);

	xdouble G  = G0* m0/m * q*Fl*Fl/q0/Fl0/Fl0; //G0 * m0/m q*Fl^2/(q0*Fl0^2)
	std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
	return std::complex<xdouble>(m0*G0,0.)/denominator;	
};

////////	////////	two_channel_breit_wigner	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class two_channel_breit_wigner : public amplitude{

	public:
		two_channel_breit_wigner();

		std::string type()											const		{return "two_channel_breit_wigner";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

two_channel_breit_wigner::two_channel_breit_wigner():amplitude(1,2,6,2){

	_name = "unnamed_two_channel_breit_wigner";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "width";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_width";

	_con_types[0] = "m_Pi";
	_con_types[1] = "m_Iso1";
	_con_types[2] = "m_Iso2";
	_con_types[3] = "branching";
	_con_types[4] = "L1";
	_con_types[5] = "L2";


	_con_names[0] = "m_Pi";
	_con_names[1] = "m_Iso1";
	_con_names[2] = "m_Iso2";
	_con_names[3] = "branching";
	_con_names[4] = "L1";
	_con_names[5] = "L2";
};

template <typename xdouble>
std::complex<xdouble> two_channel_breit_wigner::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble m0    = par[0];
	xdouble G0    = par[1];

	xdouble mPi   = con[0];
	xdouble mIso1 = con[1];
	xdouble mIso2 = con[2];
	xdouble X  = con[3];

	xdouble L1 = con[4];
	xdouble L2 = con[5];

	xdouble R = 5.;

	xdouble psl1 = psl<xdouble>(m, mPi, mIso1, R, L1);
	xdouble psl2 = psl<xdouble>(m, mPi, mIso2, R, L2);

	xdouble psl10= psl<xdouble>(m0, mPi, mIso1, R, L1);
	xdouble psl20= psl<xdouble>(m0, mPi, mIso2, R, L2);

	xdouble G = G0 * m0/m * ((1-X) * psl1/psl10 + X * psl2/psl20);

	std::complex<xdouble> denominator  = std::complex<xdouble>(m0*m0-m*m,-m0*G);
	return std::complex<xdouble>(m0*G0,0)/denominator;
};

////////	////////	vandermeulen_phase_space	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class vandermeulen_phase_space : public amplitude{

	public:
		vandermeulen_phase_space();

		std::string type()											const		{return "vandermeulen_phase_space";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

vandermeulen_phase_space::vandermeulen_phase_space():amplitude(1,1,2,3){

	_name = "unnamed_vandermeulen_phase_space";

	_var_types[0] = "m";

	_par_types[0] = "alpha";

	_par_names[0] = "unnamed_alpha";

	_con_types[0] = "m_Pi";
	_con_types[1] = "m_Iso";

	_con_names[0] = "m_Pi";
	_con_names[1] = "m_Iso";
};

template <typename xdouble>
std::complex<xdouble> vandermeulen_phase_space::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble alpha = par[0];

	xdouble mPi   = con[0];
	xdouble mIso  = con[1];

	xdouble ampor = mPi + mIso;
	std::complex<xdouble> value;

	if ( m > ampor){
		xdouble S = m*m;
		xdouble E = (S + mPi * mPi - mIso*mIso)/(2*m);
		xdouble PSQ = E*E - mPi*mPi;
		value = std::complex<xdouble>(exp(alpha*PSQ),0.);
	}else{
		value = std::complex<xdouble>(1.,0.);			
	};
	return value;
};

////////	////////	valera_dorofeev_background	///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class valera_dorofeev_background : public amplitude{

	public:
		valera_dorofeev_background();

		std::string type()											const		{return "valera_dorofeev_background";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

valera_dorofeev_background::valera_dorofeev_background():amplitude(1,2,1,4){

	_name = "unnamed_valera_dorofeev_background";

	_var_types[0] = "m";

	_par_types[0] = "alpha";
	_par_types[1] = "beta";

	_par_names[0] = "unnamed_alpha";
	_par_names[1] = "unnamed_beta";

	_con_types[0] = "m_0";

	_con_names[0] = "m_0";
};

template <typename xdouble>
std::complex<xdouble> valera_dorofeev_background::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble alpha = par[0];
	xdouble beta  = par[1];

	xdouble m0    = con[0];		
	return std::complex<xdouble>(pow((m-m0)/0.5,alpha)*exp(-beta*(m-m0-0.5)),0);
};

////////	////////	bowler_parametrization		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class bowler_parametrization : public amplitude{

	public:
		bowler_parametrization();

		std::string type()											const		{return "bowler_parametrization";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

bowler_parametrization::bowler_parametrization():amplitude(1,2,0,5){

	_name = "unnamed_bowler_parametrization";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "width";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_width";
};

template <typename xdouble>
std::complex<xdouble> bowler_parametrization::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble m0    = par[0];
	xdouble G0    = par[1];

	xdouble G = G0*	bowler_integral_table<xdouble>(m)/bowler_integral_table<xdouble>(m0) * m0/m;
	std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);

	return std::complex<xdouble>(sqrt(m0*G0),0)/denominator;
};
////////	////////	flatte				///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class flatte : public amplitude{

	public:
		flatte();

		std::string type()											const		{return "flatte";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

flatte::flatte():amplitude(1,3,2,6){

	_name = "unnamed_flatte";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "g_1";
	_par_types[2] = "g_2";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_g_1";
	_par_names[2] = "unnamed_g_2";

	_con_types[0] = "m_Pi";
	_con_types[1] = "m_K";

	_con_names[0] = "m_Pi";
	_con_names[1] = "m_K";
};

template <typename xdouble>
std::complex<xdouble> flatte::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble m0    = par[0];
	xdouble g1    = par[1];
	xdouble g2    = par[2];

	xdouble mPi   = con[0];
	xdouble mK    = con[1];
	
	xdouble qpp= breakupMomentumReal<xdouble>(m*m,mPi*mPi,mPi*mPi);
	xdouble qKK= breakupMomentumReal<xdouble>(m*m,mK*mK,mK*mK);

	std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-(g1*qpp*qpp + g2*qKK*qKK));
	return std::complex<xdouble>(1,0)/denominator;
};

////////	////////	gaus				///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class gaus : public amplitude{

	public:
		gaus();

		std::string type()											const		{return "gaus";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

gaus::gaus():amplitude(1,2,0,9){

	_name = "unnamed_gaus";

	_var_types[0] = "x";

	_par_types[0] = "mu";
	_par_types[1] = "sigma";

	_par_names[0] = "unnamed_mu";
	_par_names[1] = "unnamed_sigma";
};

template <typename xdouble>
std::complex<xdouble> gaus::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble m0    = par[0];
	xdouble sig   = par[1];

	return std::complex<xdouble>(exp(-(m-m0)*(m-m0)/2/sig/sig),0);
};
////////	////////	polynomial			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class polynomial : public amplitude{

	public:
		polynomial();

		std::string type()											const		{return "polynomial";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

polynomial::polynomial():amplitude(1,5,1,10){

	_name = "unnamed_polynomial";

	_var_types[0] = "x";

	_par_types[0] = "c0";
	_par_types[1] = "c1";
	_par_types[2] = "c2";
	_par_types[3] = "c3";
	_par_types[4] = "c4";

	_par_names[0] = "unnamed_c0";
	_par_names[1] = "unnamed_c1";
	_par_names[2] = "unnamed_c2";
	_par_names[3] = "unnamed_c3";
	_par_names[4] = "unnamed_c4";

	_con_types[0] = "dummy_con";
	
	_con_names[0] = "unnamed_dummy";
};

template <typename xdouble>
std::complex<xdouble> polynomial::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
		double m     = var[0];

		xdouble c0    = par[0];
		xdouble c1    = par[1];
		xdouble c2    = par[2];
		xdouble c3    = par[3];
		xdouble c4    = par[4];

		xdouble ret = c4*m*m*m*m + c3*m*m*m + c2*m*m + c1*m + c0;

		ret += con[0];

		return std::complex<xdouble>(ret,0);
};

////////	////////	mass_dep_bw_2			///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class mass_dep_bw_2 : public amplitude{

	public:
		mass_dep_bw_2();

		std::string type()											const		{return "mass_dep_bw_2";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

mass_dep_bw_2::mass_dep_bw_2():amplitude(1,2,6,22){

	_name = "unnamed_mass_dep_bw_2";

	_var_types[0] = "m";

	_par_types[0] = "mass";
	_par_types[1] = "width";

	_par_names[0] = "unnamed_mass";
	_par_names[1] = "unnamed_width";

	_con_types[0] = "m_Pi";
	_con_types[1] = "m_Iso1";
	_con_types[2] = "m_Iso2";
	_con_types[3] = "branching";
	_con_types[4] = "angular_momentum_1";
	_con_types[5] = "angular_momentum_2";

	_con_names[0] = "m_Pi";
	_con_names[1] = "m_Iso1";
	_con_names[2] = "m_Iso2";
	_con_names[3] = "branching";
	_con_names[4] = "L1";
	_con_names[5] = "L2";

};

template <typename xdouble>
std::complex<xdouble> mass_dep_bw_2::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];

	xdouble m0    = par[0];
	xdouble G0    = par[1];

	xdouble mPi   = con[0];
	xdouble mIso1 = con[1];
	xdouble mIso2 = con[2];
	xdouble X     = con[3];

	xdouble q1 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso1*mIso1);
	xdouble q10= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso1*mIso1);
	xdouble q2 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso2*mIso2);
	xdouble q20= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso2*mIso2);
	xdouble Fl1= barrierFactor<xdouble>(q1,con[4]);
	xdouble Fl10=barrierFactor<xdouble>(q10,con[4]);
	xdouble Fl2= barrierFactor<xdouble>(q2,con[5]);
	xdouble Fl20=barrierFactor<xdouble>(q20,con[5]);	

	xdouble G  = G0 * m0/m* ((1-X)*q1*Fl1*Fl1/q10/Fl10/Fl10 + X* q2*Fl2*Fl2/q20/Fl20/Fl20);
	std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
	return std::complex<xdouble>(m0*G0,0.)/denominator;	
};

////////	////////	t_dependent_background		///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class t_dependent_background : public amplitude{

	public:
		t_dependent_background();

		std::string type()											const		{return "t_dependent_background";};	

		template<typename xdouble>
		std::complex<xdouble> template_eval(const double* var, const xdouble* par, const xdouble* con)		const;

		std::complex<double> Eval(const double* var, const double* par, const double* con)			const		{return template_eval(var,par,con);};
#ifdef ADOL_ON
		std::complex<adouble> Eval(const double* var, const adouble* par, const adouble* con)			const		{return template_eval(var,par,con);};
#endif//ADOL_ON	
};

t_dependent_background::t_dependent_background():amplitude(2,4,3,101){

	_name = "unnamed_t_dependent_background";

	_var_types[0] = "m";
	_var_types[1] = "t'";

	_par_types[0] = "b";
	_par_types[1] = "c0";
	_par_types[2] = "c1";
	_par_types[3] = "c2";

	_par_names[0] = "unnamed_b";
	_par_names[1] = "unnamed_c0";
	_par_names[2] = "unnamed_c1";
	_par_names[3] = "unnamed_c2";

	_con_types[0] = "m_0";
	_con_types[1] = "m_Pi";
	_con_types[2] = "m_Iso";

	_con_names[0] = "m_0";
	_con_names[1] = "m_Pi";
	_con_names[2] = "m_Iso";
};

template <typename xdouble>
std::complex<xdouble> t_dependent_background::template_eval(const double* var, const xdouble* par, const xdouble* con)const{
	double m     = var[0];
	double tPrime= var[1];

	xdouble b     = par[0];
	xdouble c0    = par[1];
	xdouble c1    = par[2];
	xdouble c2    = par[3];

	xdouble m0    = con[0];
	xdouble mPi   = con[1];
	xdouble mIso  = con[2];


	xdouble PSQ = 0.;
	xdouble mpor = mPi + mIso;		
	if (m > mpor){
		xdouble E = (m*m +mPi*mPi - mIso*mIso)/(2*m);
		PSQ = E*E - mPi*mPi;
	};
	return std::complex<xdouble>(pow(m-m0,b)*exp(PSQ*(c0+c1*tPrime+c2*tPrime*tPrime)),0.);
};
// End of special definitions

amplitude* get_amplitude(int id){
	amplitude* ret_amp;
	if (id==-1){
		ret_amp = new constant_function();
	}else if (id==0){
		ret_amp = new breit_wigner();
	}else if (id==1){
		ret_amp = new mass_dep_breit_wigner();
	}else if (id==2){
		ret_amp = new two_channel_breit_wigner();
	}else if(id==3){
		ret_amp = new vandermeulen_phase_space();
	}else if(id==4){
		ret_amp = new valera_dorofeev_background();
	}else if (id==5){
		ret_amp = new bowler_parametrization();
	}else if (id ==6){
		ret_amp = new flatte();
	}else if (id == 9){
		ret_amp = new gaus();
	}else if (id==10){
		ret_amp = new polynomial();
	}else if(id == 22){
		ret_amp = new mass_dep_bw_2();
	}else if(id == 101){
		ret_amp = new t_dependent_background();
	}else{
		std::cerr<<"amplitude id not defined"<<std::endl;
		ret_amp = new amplitude();
	};
	return ret_amp;
};
#endif//BREIT_WIGNERS_WEIT_BRIGNERS
