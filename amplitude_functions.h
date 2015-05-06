// This file contains a number of definitions for amplitude parameterizations. 
// Header only, because of templates.
#ifndef AMP_FUNCTIONS_CILLY_BO
#define AMP_FUNCTIONS_CILLY_BO

#include <cmath> 
#include <vector>
#include <complex>
#include <iostream>
#include <string>
// double mPi=1.3957018;///\pm0.00035MeV // Particle Data Booklet 2012

#ifdef ADOL_ON // Some function on std::complex<adtl::adouble>, needed for automatic differentiation.
#include <adolc/adouble.h>  
namespace complex_adtl{
	std::complex<adtl::adouble> log(std::complex<adtl::adouble> z){
		adtl::adouble re = std::real(z);
		adtl::adouble im = std::imag(z);

		adtl::adouble newReal = pow(re*re+im+im,.5);
		adtl::adouble newImag = atan2(im,re);

		return std::complex<adtl::adouble>(newReal,newImag);
	};
	std::complex<adtl::adouble> sqrt(std::complex<adtl::adouble> z){
		return pow(z,0.5);
	};
	adtl::adouble abs(std::complex<adtl::adouble> z){
		adtl::adouble squared = std::real(z*std::conj(z));
		return pow(squared,0.5);
	};
	adtl::adouble sqrt(adtl::adouble x){
		return pow(x,0.5);
	};
};
#endif//ADOL_ON

const double PION_MASS 	= 0.139;
const double PI  	= 3.141592653589793238463;

//////////////////////////  SOME COMMON DEFINITIONS  //////////////////////////////////////////////////////////////////////
template< typename xdouble> xdouble breakupMomentumReal(xdouble M2, xdouble m12, xdouble m22){ // Real breakup momentum: sqrt(lambda(M,m1,m2))/(2M) or 0.
	xdouble lambda= M2*M2 + m12*m12 + m22*m22 - 2*M2*m12 -2*M2*m22 - 2*m12*m22;
	if ( lambda >= 0. ){
		return sqrt(lambda)/(2*M2);
	}else{
		std::cerr << "amplitude_functions.h: Error: Found 0 > q^2("<<M2<<","<<m12<<","<<m22<<") = "<<lambda/(4*M2*M2)<<". Sub-threshold decay."<<std::endl;
		return 0.;
	};
};

template<typename xdouble> xdouble barrierFactor(xdouble q, xdouble L){
	double pr = 0.1973;
	xdouble z=q*q/pr/pr;
	xdouble res;
	if (L == 0){
		res=1.;
	}else if (L==1){
		res=sqrt(2*z/(z+1));
	}else if (L==2){
		res=sqrt(13*z*z/((z-3)*(z-3)+9*z));
	}else if (L==3){
		res=sqrt(277*z*z*z/(z*(z-15)*(z-15)+9*pow(2*z-5,2)));
	}else if (L==4){
		res=sqrt(12746*pow(z,4)/(pow(z*z-45*z+105,2)+25*z*pow(2*z-21,2)));
	} else {
		std::cerr << "amplitude_functions.h: Error: Barrier factors not defined for L =" << L <<std::endl;
		res =0.;
	};
	return res;
};

template<typename xdouble>  // Some different definition for Barrier factors... used by Dimas program.
xdouble fdl(xdouble P, xdouble R, xdouble L){
	xdouble X = P*R;
	if (L==0){
		return 1.;
	}else if(L==1){
		return 1. + X*X;
	}else if(L==2){
		return 9. + 3.*X*X + X*X*X*X;
	}else if(L==3){
		return 225. + 45.*X*X + 6.*X*X*X*X * X*X*X*X*X*X;	
	}else if(L==4){
		return 11025. + 1575.*X*X + 135.*X*X*X*X + 10*X*X*X*X*X*X + X*X*X*X*X*X*X*X;
	}else{
		std::cerr<<"amplitude_functions.h: Error: 'fdl(...)' not defined for 4 < L = "<< L << std::endl;
		return 1.;
	};
};

template<typename xdouble>
xdouble psl(xdouble m, xdouble m1, xdouble m2, xdouble R, xdouble L){
	xdouble ampor = m1+m2;
	if (m > ampor){
		xdouble E = (m*m + m1*m1 - m2 * m2)/(2*m);
		xdouble P = pow((E*E-m1*m1)*(E*E-m1*m1),.25);
		xdouble f = fdl<xdouble>(P,R,L);
		return pow(P,2*L+1)/f;
	}else{
		return 0.;
	};
};


template<typename xdouble>
xdouble bowler_integral_table(xdouble m){
	double points[]={
0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.1313E-09, 0.8940E-07, 0.6175E-06, 0.2020E-05, 0.4786E-05, 0.9458E-05, 0.1664E-04, 0.2701E-04, 0.4131E-04, 0.6037E-04, 0.8510E-04, 0.1165E-03, 0.1558E-03, 0.2040E-03, 0.2628E-03, 0.3335E-03, 0.4179E-03, 0.5178E-03, 0.6355E-03, 0.7732E-03, 0.9337E-03, 0.1120E-02, 0.1335E-02, 0.1583E-02, 0.1867E-02, 0.2194E-02, 0.2568E-02, 0.2995E-02, 0.3483E-02, 0.4039E-02, 0.4673E-02, 0.5396E-02, 0.6220E-02, 0.7160E-02, 0.8233E-02, 0.9458E-02, 0.1086E-01, 0.1246E-01, 0.1430E-01, 0.1641E-01, 0.1884E-01, 0.2163E-01, 0.2484E-01, 0.2853E-01, 0.3277E-01, 0.3759E-01, 0.4306E-01, 0.4917E-01, 0.5591E-01, 0.6322E-01, 0.7100E-01, 0.7913E-01, 0.8752E-01, 0.9604E-01, 0.1046E+00, 0.1132E+00, 0.1218E+00, 0.1302E+00, 0.1386E+00, 0.1469E+00, 0.1551E+00, 0.1631E+00, 0.1711E+00, 0.1790E+00, 0.1867E+00, 0.1944E+00, 0.2020E+00, 0.2095E+00, 0.2169E+00, 0.2243E+00, 0.2315E+00, 0.2387E+00, 0.2458E+00, 0.2529E+00, 0.2599E+00, 0.2668E+00, 0.2737E+00, 0.2805E+00, 0.2873E+00, 0.2940E+00, 0.3007E+00, 0.3073E+00, 0.3138E+00, 0.3204E+00, 0.3269E+00, 0.3333E+00, 0.3397E+00, 0.3461E+00, 0.3525E+00, 0.3587E+00, 0.3650E+00, 0.3713E+00, 0.3775E+00, 0.3837E+00, 0.3898E+00, 0.3959E+00, 0.4020E+00, 0.4081E+00, 0.4141E+00, 0.4201E+00, 0.4261E+00, 0.4320E+00, 0.4380E+00, 0.4439E+00, 0.4498E+00, 0.4556E+00, 0.4615E+00, 0.4673E+00, 0.4731E+00, 0.4790E+00, 0.4847E+00, 0.4905E+00, 0.4962E+00, 0.5019E+00, 0.5076E+00, 0.5134E+00, 0.5189E+00, 0.5246E+00, 0.5303E+00, 0.5359E+00, 0.5415E+00, 0.5471E+00, 0.5526E+00, 0.5582E+00, 0.5638E+00, 0.5693E+00, 0.5749E+00, 0.5804E+00, 0.5858E+00, 0.5914E+00, 0.5968E+00, 0.6023E+00, 0.6077E+00, 0.6132E+00, 0.6186E+00, 0.6241E+00, 0.6294E+00, 0.6348E+00, 0.6403E+00, 0.6456E+00, 0.6510E+00, 0.6563E+00, 0.6617E+00, 0.6671E+00, 0.6724E+00, 0.6777E+00, 0.6830E+00, 0.6882E+00, 0.6936E+00, 0.6990E+00, 0.7041E+00, 0.7095E+00, 0.7149E+00, 0.7199E+00, 0.7252E+00, 0.7305E+00, 0.7356E+00, 0.7410E+00, 0.7462E+00, 0.7514E+00, 0.7567E+00, 0.7619E+00, 0.7668E+00, 0.7723E+00, 0.7774E+00, 0.7826E+00, 0.7878E+00, 0.7930E+00, 0.7982E+00, 0.8033E+00, 0.8084E+00, 0.8135E+00, 0.8188E+00, 0.8239E+00, 0.8291E+00, 0.8340E+00, 0.8393E+00, 0.8444E+00, 0.8493E+00, 0.8547E+00, 0.8597E+00, 0.8649E+00, 0.8700E+00, 0.8750E+00, 0.8800E+00, 0.8851E+00, 0.8903E+00, 0.8953E+00, 0.9005E+00, 0.9054E+00, 0.9105E+00, 0.9156E+00, 0.9205E+00, 0.9256E+00, 0.9308E+00, 0.9358E+00, 0.9408E+00, 0.9458E+00, 0.9507E+00, 0.9560E+00, 0.9609E+00, 0.9659E+00, 0.9711E+00, 0.9760E+00, 0.9808E+00, 0.9860E+00, 0.9909E+00, 0.9960E+00, 0.1001E+01, 0.1006E+01, 0.1011E+01, 0.1016E+01, 0.1021E+01, 0.1026E+01, 0.1031E+01, 0.1036E+01, 0.1041E+01, 0.1046E+01, 0.1051E+01, 0.1056E+01, 0.1061E+01, 0.1066E+01, 0.1071E+01, 0.1076E+01, 0.1081E+01, 0.1085E+01, 0.1090E+01, 0.1096E+01, 0.1100E+01, 0.1105E+01, 0.1110E+01, 0.1115E+01, 0.1120E+01, 0.1125E+01, 0.1130E+01, 0.1135E+01, 0.1140E+01, 0.1145E+01, 0.1150E+01, 0.1154E+01, 0.1160E+01, 0.1164E+01, 0.1169E+01, 0.1174E+01, 0.1179E+01, 0.1184E+01, 0.1189E+01, 0.1194E+01, 0.1199E+01, 0.1204E+01, 0.1208E+01, 0.1214E+01, 0.1218E+01, 0.1223E+01, 0.1228E+01, 0.1233E+01, 0.1238E+01, 0.1243E+01, 0.1248E+01, 0.1253E+01, 0.1257E+01, 0.1262E+01, 0.1267E+01, 0.1272E+01, 0.1277E+01, 0.1282E+01, 0.1287E+01, 0.1292E+01, 0.1296E+01, 0.1301E+01, 0.1306E+01, 0.1311E+01, 0.1316E+01, 0.1321E+01, 0.1326E+01, 0.1330E+01, 0.1336E+01, 0.1340E+01, 0.1345E+01, 0.1350E+01, 0.1355E+01, 0.1359E+01, 0.1365E+01, 0.1369E+01, 0.1374E+01, 0.1379E+01, 0.1384E+01, 0.1389E+01, 0.1394E+01, 0.1398E+01, 0.1404E+01, 0.1408E+01, 0.1412E+01, 0.1418E+01, 0.1422E+01, 0.1427E+01, 0.1432E+01, 0.1437E+01, 0.1442E+01, 0.1447E+01, 0.1451E+01, 0.1457E+01, 0.1461E+01, 0.1466E+01, 0.1472E+01, 0.1475E+01, 0.1480E+01, 0.1486E+01, 0.1490E+01, 0.1495E+01, 0.1500E+01, 0.1504E+01, 0.1510E+01, 0.1514E+01, 0.1518E+01, 0.1524E+01, 0.1529E+01, 0.1534E+01, 0.1538E+01, 0.1542E+01, 0.1549E+01, 0.1552E+01, 0.1557E+01, 0.1562E+01, 0.1567E+01, 0.1573E+01, 0.1577E+01, 0.1581E+01, 0.1586E+01, 0.1592E+01, 0.1595E+01, 0.1601E+01, 0.1605E+01, 0.1610E+01, 0.1616E+01, 0.1619E+01, 0.1625E+01, 0.1630E+01, 0.1634E+01, 0.1639E+01, 0.1644E+01, 0.1648E+01, 0.1654E+01, 0.1658E+01, 0.1663E+01, 0.1668E+01, 0.1672E+01, 0.1678E+01, 0.1682E+01, 0.1687E+01, 0.1692E+01, 0.1696E+01, 0.1701E+01, 0.1708E+01, 0.1710E+01, 0.1716E+01, 0.1721E+01, 0.1724E+01, 0.1726E+01};
	double xmin = 0.;
	double xmax = 4.;
	double step = 0.01;
	if (m > xmin and m < xmax){
		int nStep = 0;
		while (xmin+(nStep+1)*step < m){ // nStep is the last integer, where the mass is smaller than the input mass.
			nStep+=1;
		};
		double upper = points[nStep+1];
		double lower = points[nStep];
		double mUpper = xmin + (nStep+1)*step;
		double mLower = xmin + nStep*step;
		xdouble x = (m - mLower)/(mUpper-mLower);
		return (1-x)*lower + x*upper;
	}else{
		return 0.;
	};
};

/// Basic class definition for the functor class //////////////////////////////////////////////////////////////

class amplitude {
	public:
		amplitude();
		amplitude(size_t nVar, size_t nPar, size_t nCon, int funcId);

		virtual std::string 	type()									const		{return "constant_function";};
		std::string		name()									const		{return _name;};

		std::complex<double> Eval(const double* var)							const		{return Eval(var, &_parameters[0], &_constants[0]);};
		std::complex<double> Eval(const double* var, const double* par)					const		{return Eval(var, par, &_constants[0]);};
		virtual std::complex<double> Eval(const double* var, const double* par, const double* con)	const		{return std::complex<double>(1.,0.);};

#ifdef ADOL_ON
		std::complex<adtl::adouble> Eval(const double* var, const adtl::adouble* par)				const		{return Eval(var, par, &_constants[0]);};
		virtual std::complex<adtl::adouble> Eval(const double* var, const adtl::adouble* par, const double* con)	const		{return std::complex<adtl::adouble>(1.,0.);};
#endif//ADOL_ON

		size_t 			nVar()									const		{return _nVar;};
		size_t 			nPar()									const		{return _nPar;};
		size_t 			nCon()									const		{return _nCon;};

		size_t			L()									const		{return _L;};
		void			setL(size_t L)										{_L=L;};

		bool 			setPar(size_t n, double val);
		bool 			setCon(size_t n, double val);

		bool			setPars(const double* vals);
		bool			setCons(const double* vals);

		const std::vector<double>*parameters()								const		{return &_parameters;};
		const std::vector<double>*constants()								const		{return &_constants;};
	
		double			getParameter(size_t n)							const;
		double			getConstant(size_t n)							const;

		std::string		getVarType(size_t n)							const;					
		std::string		getParType(size_t n)							const;
		std::string		getConType(size_t n)							const;
		std::string		getParName(size_t n)							const;
		std::string		getConName(size_t n)							const;

		int			getFuncId()								const		{return _funcId;};

		void			setFunctionName(std::string name)							{_name = name;};
		bool			setParName(size_t n, std::string name);
		bool			setConName(size_t n, std::string name);

		void			print()									const;

	protected:

		size_t 			_L;		// Spin, treat special here, since its essential

		std::vector<std::string>_var_types;	// Type of the variables // Convention: {m, t',...}
		std::vector<std::string>_par_types;	// Names of the parameters
		std::vector<std::string>_con_types;	// Names of the constants

		std::string		_name;		// Name of the functions
		std::vector<std::string>_par_names;	// Names of the parameters
		std::vector<std::string>_con_names;	// Names of the constants

		const size_t		_nVar;		// # of Variables
		const size_t		_nPar;		// # of Paramters
		const size_t		_nCon;		// # of Constants

		const int		_funcId;	// Id # of the function type, to be consistent with the old method

		std::vector<double>	_parameters;	// Paramters
		std::vector<double>	_constants;	// Constants

};

amplitude::amplitude():_nVar(0),_nPar(0),_nCon(0),_funcId(-1){ // Id for constant function is -1

		_name = "unnamed_constant_function";
};

amplitude::amplitude(size_t nVar, size_t nPar, size_t nCon, int funcId):_nVar(nVar), _nPar(nPar), _nCon(nCon), _funcId(funcId){

		_name = "unnamed_constant_function";

		_constants = std::vector<double>(nCon);
		_parameters= std::vector<double>(nPar);

		_var_types = std::vector<std::string>(nVar);

		_par_types = std::vector<std::string>(nPar);
		_par_names = std::vector<std::string>(nPar);

		_con_types = std::vector<std::string>(nCon);
		_con_names = std::vector<std::string>(nCon);
};

bool amplitude::setPar(size_t n, double val){
	if (n < _nPar){
		_parameters[n] = val;
		return true;
	}else{
		std::cerr << "Error: Can't set parameter #"<<n<<" for "<<type()<<std::endl;		
		return false;
	};
};

bool amplitude::setCon(size_t n, double val){
	if (n < _nCon){
		_constants[n] = val;
		return true;
	}else{
		std::cerr << "Error: Can't set constant #"<<n<<" for "<<type()<<std::endl;		
		return false;
	};
};

bool amplitude::setPars(const double* vals){
	bool ret = true;
	for (size_t iii =0;iii<_nPar;iii++){
		if (not setPar(iii,vals[iii])){
			ret = false;
		};
	};
	return ret;
};

bool amplitude::setCons(const double* vals){
	bool ret = true;
	for (size_t iii =0;iii<_nCon;iii++){
		if (not setCon(iii,vals[iii])){
			ret = false;
		};
	};
	return ret;
};

std::string amplitude::getVarType(size_t n)		const{
	if (n<_nVar){
		return _var_types[n];
	}else{
		std::cerr <<"Error: getVarType(...): _nVar  = "<<_nVar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getParType(size_t n)		const{
	if (n<_nPar){
		return _par_types[n];
	}else{
		std::cerr <<"Error: getParType(...): _nPar  = "<<_nPar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getConType(size_t n)		const{
	if (n<_nCon){
		return _con_types[n];
	}else{
		std::cerr <<"Error: getConType(...): _nCon  = "<<_nCon<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getParName(size_t n)		const{
	if (n<_nPar){
		return _par_names[n];
	}else{
		std::cerr <<"Error: getParName(...): _nPar  = "<<_nPar<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
std::string amplitude::getConName(size_t n)		const{
	if (n<_nCon){
		return _con_names[n];
	}else{
		std::cerr <<"Error: getConName(...): _nCon  = "<<_nCon<<" < "<<n<<std::endl;
		return "index_out_of_range";
	};
};
bool amplitude::setParName(size_t n, std::string name){
	if(n<_nPar){
		_par_names[n] = name;
		return true;
	}else{
		std::cerr<<"Error: Can't set _par_names["<<n<<"] for "<<type()<<std::endl;
		return false;
	};
};
bool amplitude::setConName(size_t n, std::string name){
	if(n<_nCon){
		_con_names[n] = name;
		return true;
	}else{
		std::cerr<<"Error: Can't set _con_names["<<n<<"] for "<<type()<<std::endl;
		return false;
	};
};

double amplitude::getParameter(size_t n)		const{
	if(n<_nPar){
		return _parameters[n];
	}else{
		std::cerr<<"Error: Can't get _parameters["<<n<<"] for "<<type()<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};

double amplitude::getConstant(size_t n)			const{
	if(n<_nCon){
		return _constants[n];
	}else{
		std::cerr<<"Error: Can't get _contants["<<n<<"] for "<<type()<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};

void amplitude::print()					const{
	std::cout<<type()<<": "<<name()<<" ("<<_funcId<<")"<<std::endl;
	if(_nVar>0){
		std::cout<<"  - variables: "<<std::endl;
		for (size_t i=0;i<_nVar;i++){
			std::cout<<"    - "<<_var_types[i]<<std::endl;
		};
	};
	if(_nPar>0){
		std::cout<<"  - parameters: "<<std::endl;
		for (size_t i=0;i<_nPar;i++){
			std::cout<<"    - "<<_par_types[i]<<": "<<_par_names[i]<<": "<<_parameters[i]<<std::endl;
		};
	};
	if(_nCon>0){
		std::cout<<"  - constants: "<<std::endl;
		for (size_t i=0;i<_nCon;i++){
			std::cout<<"    - "<<_con_types[i]<<": "<<_con_names[i]<<": "<<_constants[i]<<std::endl;
		};
	};
};

#endif//AMP_FUNCTIONS_CILLY_BO
