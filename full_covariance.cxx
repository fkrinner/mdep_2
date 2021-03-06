#include"full_covariance.h"
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<string>
#include<limits>
#include<stdexcept>

#include"matrix_utilities.h"
#include"string_utilities.h"

#ifdef ADOL_ON
namespace adtl{
	int ADOLC_numDir = NUMBER_DIRECTIONS;
};
#include "adolc/adouble.h"
std::complex<adtl::adouble> operator*(std::complex<adtl::adouble> a,double d){// Define std::complex<adtl::adouble> * double
	return std::complex<adtl::adouble>(a.real()*d,a.imag()*d);
};
#endif//ADOL_ON

full_covariance::full_covariance():
	full_SDM(){
	
	setTbinning((*_waveset.t_binning()));
};
//########################################################################################################################################################
#ifdef USE_YAML
///Constructror from YAML files
full_covariance::full_covariance(
							std::string 						card):

	full_SDM(card){
	YAML::Node Ycard   = YAML::LoadFile(card);
	std::string parametrizations = get_relative_path(Ycard["parametrization_file"].as<std::string>(),_waveset.YAML_file());
	YAML::Node Yparam  = YAML::LoadFile(parametrizations);
	update_n_cpls();
	setTbinning((*_waveset.t_binning()));
	update_definitions();
	update_min_max_bin();
	std::cout<<"full_covariance::full_covariance(...): Load full_covariance from YAML file\nLoad data and coma"<<std::endl;
	loadDataComa(Ycard);
	std::cout<<"full_covariance::full_covariance(...): Data and coma loaded\nLoad parameter values"<<std::endl;
	loadParameterValues(Ycard, Yparam);
	std::cout<<"full_covariance::full_covariance(...): Paramter values loaded"<<std::endl;
	std::cout<<"full_covariance::full_covariance(...): full_covariance loaded"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Evaluate Chi2 with the paramters xx
double full_covariance::mainEval(
							const double							*xx){


//	const std::complex<double>* cpl = (std::complex<double>*)xx; // This is forbidden by Charly, build complex variables by hand :(
	std::complex<double> cpl[_nCpl];
	for (size_t i=0;i<_nCpl;i++){
		cpl[i] = std::complex<double>(xx[2*i],xx[2*i+1]);
	};
	const double* par = xx+2*_nCpl;
//	const std::complex<double>* bra = (std::complex<double>*)(xx+2*_nCpl+_nPar);// This is forbidden by Charly, build complex variables by hand :(
	std::complex<double> bra[_nBra];
	for (size_t i=0;i<_nBra;i++){
		bra[i] = std::complex<double>(xx[2*_nCpl+_nPar+2*i],xx[2*_nCpl+_nPar+2*i+1]);
	};
	const double* iso_par = xx + 2*_nCpl + _nPar + 2*_nBra;

	// std::cout<<par[0]<<std::endl;
	double chi2;
	chi2 = EvalBranch(bra,cpl,par,iso_par);
	_count++;
	/* if (chi2 != chi2){ // Check for NaN
	std::cout<<"NaN paramters are:"<<std::endl;
	std::cout<<"bra:"<<std::endl;
	print_vector(bra);
	std::cout<<"cpl:"<<std::endl;
	print_vector(cpl);
	std::cout<<"par:"<<std::endl;
	print_vector(par);
	std::cout << xx << std::endl;
	throw;
	};*/
	if (0==_count%_nOut){ // Write every _nOut evaluation
		std::cout<<"#"<<_count<<": "<<chi2<<std::endl;
	};
	if (0==_count%_nOutFile){
		if (not _parameterFile.size() == 0){
			writeParameters(xx,_parameterFile);
		};
	};
	return chi2;
};
//########################################################################################################################################################
///Evaluates chi2 w/o branchings (EvalC(ouplings)P(arameters))
template<typename xdouble>
xdouble full_covariance::EvalCP(
							const std::complex<xdouble>					*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	xdouble chi2 = 0.;
	std::vector<std::complex<xdouble> > actCpl = std::vector<std::complex<xdouble> >(_waveset.nFtw());
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if((*_waveset.eval_tbin())[tbin]){
			for (size_t i=0;i<_waveset.nFtw();i++){
				actCpl[i] = cpl[i+tbin*_waveset.nFtw()];
			};
//			std::cout<<"VVVVVVVVVVVVVVVVVVVVVVVVV"<<std::endl;
//			print_vector(actCpl);
//			std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
			chi2+=EvalTbin(tbin,&actCpl[0],&par[0],&iso_par[0]);//[0]//
	//		std::cout<<tbin<<":E:"<<EvalTbin(tbin,actCpl,par)<<std::endl;
//			std::cout<<chi2<<std::endl;
		};
	};
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double full_covariance::EvalCP(const std::complex<double> *cpl,const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Evaluates chi2 with branchings
template<typename xdouble>
xdouble full_covariance::EvalBranch(
							const std::complex<xdouble>					*branch,
							const std::complex<xdouble>	 				*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	std::vector<std::complex<xdouble> > cplin(_waveset.nFtw()*_waveset.nTbin());
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if((*_waveset.eval_tbin())[tbin]){
			for (size_t i =0;i<_waveset.nFtw();i++){
				if ((*_waveset.n_branch())[i]==-1){
					cplin[i+tbin*_waveset.nFtw()] = cpl[(*_waveset.n_cpls())[i]+tbin*(_waveset.nBrCpl())];
				}else{
					cplin[i+tbin*_waveset.nFtw()] = cpl[(*_waveset.n_cpls())[i]+tbin*_waveset.nBrCpl()] * branch[(*_waveset.n_branch())[i]];
				};
			};
		};
	};
	xdouble chi2 = EvalCP(&cplin[0],par, iso_par);//&...[0] stays
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double full_covariance::EvalBranch(const std::complex<double> *branch, const std::complex<double> *cpl, const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets the chi2 for a single t' bin
double full_covariance::EvalTbin(
							int 								tbin,
							const std::complex<double>	 				*cpl,
							const double	 						*par,
							const double 							*iso_par)		const{

	double chi2 = 0.;
	std::vector<std::vector<std::complex<double> > > iso_eval = _waveset.iso_funcs(iso_par);
//Careful, multiprocessing messes up the output file!!! Bed for debugging
#pragma omp parallel for reduction(+:chi2)
	for (size_t bin=_waveset.minBin(); bin<_waveset.maxBin(); bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
#ifdef ADOL_ON // Has to be twice, otherwise the openmp pragma would not work. Did not fins another workaround
adtl::adouble full_covariance::EvalTbin(
							int 								tbin,
							const std::complex<adtl::adouble>	 				*cpl,
							const adtl::adouble	 						*par,
							const adtl::adouble 							*iso_par)		const{

	adtl::adouble chi2 = 0.;
	std::vector<std::vector<std::complex<adtl::adouble> > > iso_eval = _waveset.iso_funcs(iso_par);
	for (size_t bin=_waveset.minBin(); bin<_waveset.maxBin(); bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
#endif//ADOL_ON
//########################################################################################################################################################
///Gets the Chi2 for a single t' and m3pi bin
template<typename xdouble>
xdouble full_covariance::EvalBin(
							int 								tbin,
							int 								bin,
							const std::complex<xdouble>					*cpl,
							const xdouble	 						*par,
							std::vector<std::vector<std::complex<xdouble> > > 		&iso_eval)		const{

	double mass = _waveset.get_m(bin); // Eval at bin center.

//	std::vector<std::complex<xdouble> > ccpl(_waveset.nFtw(),std::complex<xdouble>(1.,0.));
//	std::cout<<ccpl.size()<<std::endl;
//	std::cout<<_nCpl<<std::endl;
//	for (size_t i=0; i<_nPar;++i){
//		std::cout<<i<<" "<<par[i]<<std::endl;
//	};


	std::vector<xdouble> deltas = delta(tbin,bin,mass, cpl, par,iso_eval);


	xdouble chi2 = 0.;
//	print_vector(deltas);
	for (size_t i=0;i<_waveset.nPoints()*_waveset.nPoints();++i){
//		int iWave1 = (*_waveset.point_to_wave())[(i+1)/2];
//		if (mass >= _lowerLims[iWave] and mass < _upperLims[iWave]){
//		if (_is_active[bin][i]){
		if (deltas[i] != 0.){

//			std::cout<<"le_addite:D"<<pow(deltas[i],2.)*_coma[tbin][bin][i][i]<<std::endl;
			if (_waveset.write_out()){
				*_waveset.outStream() <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[i]<<"        "<<deltas[i]<<std::endl;
				*_waveset.outStream() <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][i][i]<<std::endl;
				*_waveset.outStream() <<"  Ci^2 before      "<<chi2<<"      +     "<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<"       =    ";
			};

//			if (tbin == 3 and bin == 50){
//				std::cout<< i<<" "<<deltas[i]<<" "<<_coma[tbin][bin][i][i]<<std::endl;
//			};

			chi2+= deltas[i]*deltas[i]*_coma[tbin][bin][i][i];
//			std::cout<<deltas[i]<<"*"<<deltas[i]<<"*"<<_coma[tbin][bin][i][i]<<"="<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<std::endl;
			if(_waveset.write_out()){
				*_waveset.outStream() <<chi2<<std::endl;
			};

			for (size_t j=0;j<i;j++){
//				int jWave = (*_waveset.point_to_wave())[(j+1)/2];
//				if(mass >= _lowerLims[jWave] and mass < _upperLims[jWave]){
//				if(_is_active[bin][j]){
				if (deltas[j] != 0.){
//					std::cout<<"le_addite: "<<2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if (_waveset.write_out()){
						*_waveset.outStream() <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<j+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
						*_waveset.outStream() <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[j]<<"        "<<deltas[i]<<std::endl;
						*_waveset.outStream() <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][j][i]<<std::endl;
						*_waveset.outStream() <<"  Ci^2 before      "<<chi2<<"      +     "<<2.*deltas[j]*deltas[i]*_coma[tbin][bin][i][j]<<"       =    ";
					};
					chi2+=2.*deltas[i]*deltas[j]*_coma[tbin][bin][i][j]; // Factor 2. because _coma[t][m][i][j] is symmetric under i<->j
//					std::cout<<deltas[i]<<"*"<<deltas[j]<<"*"<<_coma[tbin][bin][i][j]<<"="<<deltas[i]*deltas[j]*_coma[tbin][bin][i][j]<<std::endl;
					if(_waveset.write_out()){
						*_waveset.outStream() <<chi2<<std::endl;
					};
				};
			};
		};
	};
//	std::cout<<tbin<<":"<<bin<<"::"<<chi2<<std::endl;
	return chi2;
};
template double full_covariance::EvalBin(int tbin,int bin,const std::complex<double> *cpl,const double *par,std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
///Returns f(m,...) - data[...] for each SDM entry in the fit
template<typename xdouble>
std::vector<xdouble> full_covariance::delta(
							int 								tbin,
							int 								bin,
							double 								mass,
							const std::complex<xdouble> 					*cpl,
							const xdouble	 						*par,
							std::vector<std::vector<std::complex<xdouble> > > 		&iso_eval)		const{

	size_t nPoints = _waveset.nPoints();
	std::vector<double> var = _waveset.getVar(mass,tbin);
	std::vector<std::complex<xdouble> > ampls = _waveset.amps(&var[0], cpl, par, iso_eval);


///	std::cout<<"''''''''''''''''''''''''''''''''''''''''"<<std::endl;
///	for (size_t i=0;i<ampls.size();++i){
///		std::cout<<i<<": "<<ampls[i]<<std::endl;
///	};
	std::vector<xdouble> del(nPoints*nPoints,0.);
//	del[0]=std::norm(ampls[0]) - _data[tbin][bin][0]; ////// HERE ADD if(_is_point_bin[bin][i]){...}; so that no deltas are aclculated for turend off bins
	for (size_t i = 0; i<_waveset.nPoints();i++){
//		if (not _is_active[bin][i]){

//			continue;
//		};
		if (ampls[i] == std::complex<xdouble>(0.,0.)){
			continue;
		};

		del[i*(nPoints+1)] = std::norm(ampls[i]) - _data[tbin][bin][i*(nPoints+1)];
		if(_waveset.write_out()){
			*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
			*_waveset.outStream() << " Re1 data    "<<_data[tbin][bin][i*(nPoints+1)]<<"     - theory   "<<std::norm(ampls[i])<<"       =    "<<del[i*(nPoints+1)]<<std::endl;
		};
		for (size_t j=0; j<i;++j){
//			if (not _is_active[bin][j]){
//				continue;
//			};
			
			if (ampls[j] == std::complex<xdouble>(0.,0.)){
				continue;
			};

			std::complex<xdouble> inter = ampls[i]*std::conj(ampls[j]);
			del[nPoints*i+j]=   imag(inter) - _data[tbin][bin][nPoints*i+j];    // imag part
			del[nPoints*j+i]=   real(inter) - _data[tbin][bin][nPoints*j+i];    // real part

			if(_waveset.write_out()){
				*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<nPoints*i+j<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() << " Re data    "<<_data[tbin][bin][nPoints*i+j]<<"     - theory   "<<real(inter)<<"       =    "<<del[nPoints*i+j]<<std::endl;
				*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<nPoints*j+i<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() << " Im data    "<<_data[tbin][bin][nPoints*j+i]<<"     - theory   "<<imag(inter)<<"       =    "<<del[nPoints*j+i]<<std::endl;
			};
		};
	};
//	print_vector(del);	

	return del;
};
template std::vector<double> full_covariance::delta(int tbin, int bin,double mass, const std::complex<double> *cpl, const double *par, std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
///Instantiate auto diff methods, if needed (Enable adtl::adouble operations, if the auto diff package is loaded)
#ifdef ADOL_ON
template adtl::adouble full_covariance::EvalCP(const std::complex<adtl::adouble> *cpl,const adtl::adouble *par, const adtl::adouble *iso_par) const;
template adtl::adouble full_covariance::EvalBranch(const std::complex<adtl::adouble> *branch, const std::complex<adtl::adouble> *cpl, const adtl::adouble *par, const adtl::adouble *iso_par) const;
template adtl::adouble full_covariance::EvalBin(int tbin,int bin,const std::complex<adtl::adouble> *cpl,const adtl::adouble *par, std::vector<std::vector<std::complex<adtl::adouble> > > &iso_par) const;
template std::vector<adtl::adouble> full_covariance::delta(int tbin, int bin,double mass, const std::complex<adtl::adouble> *cpl, const adtl::adouble *par, std::vector<std::vector<std::complex<adtl::adouble> > > &iso_eval) const;
//#######################################################################################################################################################
///Gets the gradient w.r.t. xx
std::vector<double> full_covariance::Diff(
                                                        const std::vector<double>                                       &xx)                            const{

	std::vector<adtl::adouble>aCpl_r(2*_nCpl);
	std::vector<adtl::adouble>aPar(_nPar);
	std::vector<adtl::adouble>aBra_r(2*_nBra);
	std::vector<adtl::adouble>aIso(_nIso);
	int count=0;
	for (size_t i=0;i<2*_nCpl;i++){
		aCpl_r[i] = xx[count];
		aCpl_r[i].setADValue(count,1.0);
		++count;
	};
	for (size_t i=0; i<_nPar;i++){
		aPar[i] = xx[count];
		aPar[i].setADValue(count,1.0);
		++count;
	};
	for (size_t i=0;i<2*_nBra;i++){
		aBra_r[i] = xx[count];
		aBra_r[i].setADValue(count,1.0);
		++count;
	};
	for (size_t i=0;i<_nIso;i++){
		aIso[i] = xx[count];
		aIso[i].setADValue(count,1.0);
		++count;
	};
	std::vector<std::complex<adtl::adouble> > aCpl_c(_nCpl);
	std::vector<std::complex<adtl::adouble> > aBra_c(_nBra);
	for (size_t i=0;i<_nCpl;i++){
		aCpl_c[i] = std::complex<adtl::adouble>(aCpl_r[2*i],aCpl_r[2*i+1]);
	};
	for (size_t i=0;i<_nBra;i++){
		aBra_c[i] = std::complex<adtl::adouble>(aBra_r[2*i],aBra_r[2*i+1]);
	};
	adtl::adouble aChi2;
	aChi2 = EvalBranch(&aBra_c[0],&aCpl_c[0],&aPar[0],&aIso[0]);//[0]//
	vector<double> gradient(_nTot);
	for (size_t i=0;i<_nTot;i++){
		gradient[i]=aChi2.getADValue(i);
	};
	return gradient;
};
#endif//ADOL_ON
//########################################################################################################################################################
///Set coma manually
bool full_covariance::set_coma(
							int 								tbin,
							int 								bin,
							std::vector<std::vector<double> > 				coma){

	//print_matrix(coma);
	// Ensure right number of bins
	if(_coma.size() != _waveset.nTbin()){
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
	if(_coma[tbin].size() != _waveset.nBins()){
		_coma[tbin] = std::vector<std::vector<std::vector<double> > >(_waveset.nBins());
	};
	_coma[tbin][bin] = coma;
	if (coma.size() == _waveset.nPoints()*_waveset.nPoints()){
		for (size_t i=0;i<coma.size();i++){
			if (coma[i].size() != coma.size()){
				std::cout<<"full_covariance::set_coma(...): Warning: Set coma is not quadratic."<<std::endl;
				return false;
			};
		};
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Loads the (inverted) covariance matrix for a t' bin from 'comaFile'
void full_covariance::loadComa(
							int 								tbin,
							const char* 							comaFile){

	_coma[tbin]=std::vector<std::vector<std::vector<double> > >();
	std::fstream data(comaFile,std::ios_base::in);
	double val;
	size_t nDat = _waveset.nPoints()*_waveset.nPoints();
	std::vector<std::vector<double> > coma_bin;
	std::vector<double> coma_line;
	while(data>>val){
		coma_line.push_back(val);
		if (coma_line.size() == nDat){
			coma_bin.push_back(coma_line);
			coma_line = std::vector<double>();
		};
		if (coma_bin.size() == nDat){
			_coma[tbin].push_back(coma_bin);
			coma_bin = std::vector<std::vector<double> >();
		};
	};
	if (_waveset.nBins() != _coma[tbin].size()){
		std::cout << "full_covariance::loadComa(...): Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _coma.size()="<<_coma[tbin].size() << std::endl;
	}else{
		std::cout<< "full_covariance::loadComa(...): File delivered the right size for _coma"<<std::endl;
	};
};
//########################################################################################################################################################
///Complex conjugates _data, changes _coma accordingly
void full_covariance::conjugate(){ //<<need>>

	size_t nPoints = _waveset.nPoints();
	for (size_t ii = 0;ii<nPoints*nPoints;++ii){
		size_t i = ii/nPoints;
		size_t j = ii-nPoints*i;
		if (j<i){
			for(size_t tbin=0;tbin<_waveset.nTbin();++tbin){
				for(size_t bin=0;bin<_waveset.nBins();++bin){
//					_data[tbin][bin][ii]*=-1;
					for (size_t jj=0;jj<nPoints*nPoints;++jj){
						_coma[tbin][bin][ii][jj]*=-1;
						_coma[tbin][bin][jj][ii]*=-1;
					};
				};
			};
		};
	};
	for (size_t i=0;i<nPoints;++i){
		for (size_t j=0;j<i;++j){
			for(size_t tbin=0;tbin<_waveset.nTbin();++tbin){
				for(size_t bin=0;bin<_waveset.nBins();++bin){
					_data[tbin][bin][i*nPoints+j]*=-1;
				};
			};
		};
	};

};
//########################################################################################################################################################
///Initializes _data and _coma, when setting the t' binning
void full_covariance::setTbinning(std::vector<std::vector<double> > binning){
	_waveset.setTbinning(binning);
	if (_data.size() != _waveset.nTbin()){
		if (_data.size() != 0){
			std::cout<<"full_covariance::setTbinning(...): Warning: _data.size() != _waveset.nTbin(), but nonzero. Previous set _data[][][] will be lost."<<std::endl;
		};
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_coma.size() != _waveset.nTbin()){
		if (_coma.size() != 0){
			std::cout<<"full_covariance::setTbinning(...): Warning: _coma.size() != _waveset.nTbin(), bun nonzero. Previous set _coma[][][][] will be lost."<<std::endl;
		};
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
};
//########################################################################################################################################################
///Writes plots for one t' bin, only diagonal elements at the moment
void full_covariance::write_plots(
							std::string						filename,
							int 							tbin,
							const std::vector<std::complex<double> >		&cpl,
							const std::vector<double>				&par,
							const std::vector<std::complex<double> >		&bra,
							const std::vector<double> 				&iso)					const{

	std::cout<<std::endl;
	std::vector<std::complex<double> > cpl_all = getAllCouplings(tbin,cpl,par,bra,iso);
	std::cout<<std::endl;
	std::vector<std::vector<std::complex<double> > > iso_eval;
	if(_waveset.has_isobars()){
		iso_eval = _waveset.iso_funcs(&iso[0]);
	};
	std::ofstream write_out;
	write_out.open(filename.c_str());
	std::cout<<"full_covariance::write_plots(...): write_plots(...): Chi2 for the used paramters is: "<<EvalTbin(tbin,&cpl_all[0],&par[0],&iso[0])<<std::endl;//[0]//
	for (size_t bin=0;bin<_waveset.nBins();bin++){
		double mass = _waveset.get_m(bin);
		std::vector<double> var = _waveset.getVar(mass,tbin);
		std::vector<std::complex<double> > amplitudes = _waveset.amps(&var[0],&cpl_all[0],&par[0],iso_eval);
		write_out<<mass<<" ";
		for (size_t i=0; i<_waveset.nPoints(); i++){
			for (size_t j=0; j<_waveset.nPoints(); j++){
				std::complex<double> inter = amplitudes[i]*conj(amplitudes[j]);
				if (i==j){
					write_out<<inter.real()<<" "<<_data[tbin][bin][(_waveset.nPoints()+1)*i]<<" "<<pow(_coma[tbin][bin][(_waveset.nPoints()+1)*i][(_waveset.nPoints()+1)*i],-.5)<<" ";
				};
				if (i<j){
					write_out<<inter.real()<<" "<<_data[tbin][bin][_waveset.nPoints()*i+j]<<" "<<pow(_coma[tbin][bin][_waveset.nPoints()*i+j][_waveset.nPoints()*i+j],-.5)<<" ";
				};
				if(i>j){
					write_out<<inter.imag()<<" "<<_data[tbin][bin][_waveset.nPoints()*i+j]<<" "<<pow(_coma[tbin][bin][_waveset.nPoints()*i+j][_waveset.nPoints()*i+j],-.5)<<" ";
				};
			};
		};
		write_out<<std::endl;
	};
	write_out.close();
};
//########################################################################################################################################################
