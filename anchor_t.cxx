#include"anchor_t.h"
#include<vector>
#include<complex>
#include<fstream>
#include<iostream>
#include<string>
#include<limits>

#include"matrix_utilities.h"

#ifdef ADOL_ON
#include "adolc/adolc.h"
std::complex<adouble> operator*(std::complex<adouble> a,double d){// Define std::complex<adouble> * double
	return std::complex<adouble>(a.real()*d,a.imag()*d);
};
#endif//ADOL_ON

anchor_t::anchor_t():

	method(),
	_is_ampl(false),
	_useBranch(true){
	setTbinning((*_waveset.t_binning()));
};
//########################################################################################################################################################
#ifdef USE_YAML
///Constructror from YAML files
anchor_t::anchor_t(
							std::string 						card):

	method(card),
	_is_ampl(false),
	_useBranch(true){
	YAML::Node Ycard   = YAML::LoadFile(card);
	std::string parametrizations = Ycard["parametrization_file"].as<std::string>();
	YAML::Node Yparam  = YAML::LoadFile(parametrizations);
	update_n_cpls();
	setTbinning((*_waveset.t_binning()));
	update_definitions();
	update_min_max_bin();
	std::cout<<"anchor_t::anchot_t(...): Load anchor_t from YAML file\nLoad data and coma"<<std::endl;
	loadDataComa(Ycard);
	std::cout<<"anchor_t::anchot_t(...): Data and coma loaded\nLoad parameter values"<<std::endl;
	loadParameterValues(Ycard, Yparam);
	std::cout<<"anchor_t::anchot_t(...): Paramter values loaded"<<std::endl;
	std::cout<<"anchor_t::anchot_t(...): anchor_t loaded"<<std::endl;
};
#endif//USE_YAML
//########################################################################################################################################################
///Evaluate Chi2 with the paramters xx
double anchor_t::mainEval(
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
	if (_useBranch){ // Evaluate
		chi2 = EvalAutoCplBranch(bra,cpl,par,iso_par);
	}else{
		chi2 = EvalAutoCpl(cpl,par,iso_par); // This only works, because of the Automatic coupling finding algorithm, switching off the branchings does not give additional couplings
	};
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
xdouble anchor_t::EvalCP(
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
template double anchor_t::EvalCP(const std::complex<double> *cpl,const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Evaluates chi2 with branchings
template<typename xdouble>
xdouble anchor_t::EvalBranch(
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
template double anchor_t::EvalBranch(const std::complex<double> *branch, const std::complex<double> *cpl, const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets Chi2 with automatic non anchor couplings (no branchings)
template<typename xdouble>
xdouble anchor_t::EvalAutoCpl(
							const std::complex<xdouble>	 				*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	xdouble chi2 = 0.;
	int nNon = (*_waveset.borders_waves())[0];
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if((*_waveset.eval_tbin())[tbin]){
			const std::complex<xdouble>* actCpl = cpl+tbin*nNon;
			chi2+=EvalAutoCplTbin(tbin,actCpl,par,iso_par);
	//		std::cout << tbin<<":A:"<< EvalAutoCplTbin(tbin,actCpl,par)<<std::endl;
		};
	};
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalAutoCpl(const std::complex<double> *cpl,const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets Chi2 for automatically calculated couplings with branchings (the BEST)
template<typename xdouble>
xdouble anchor_t::EvalAutoCplBranch(
							const std::complex<xdouble>					*bra,
							const std::complex<xdouble>					*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	xdouble chi2=0.;
	int nNon = _nBrCplAnc;
	for(size_t tbin=0;tbin<_waveset.nTbin();tbin++){
		if ((*_waveset.eval_tbin())[tbin]){
			const std::complex<xdouble>* actCpl = cpl+tbin*nNon;
			std::vector<std::complex<xdouble> > bestcpl = getMinimumCplBra(tbin,bra,actCpl,par,iso_par);
			std::vector<std::complex<xdouble> > best_cpl_std = std::vector<std::complex<xdouble> >(_waveset.nFtw());
			for (size_t i=0;i<_waveset.nFtw();i++){
				if(-1==(*_waveset.n_branch())[i]){
					best_cpl_std[i] = bestcpl[(*_waveset.n_cpls())[i]];
				}else{
					best_cpl_std[i] = bestcpl[(*_waveset.n_cpls())[i]] * bra[(*_waveset.n_branch())[i]];
				};
			};
			chi2+=EvalTbin(tbin,&best_cpl_std[0],par,iso_par);//&...[0] stays ...
		};
	};
	if (_waveset.write_out()){
		*_waveset.outStream() <<"  Ci^2 final      "<<chi2<<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalAutoCplBranch(const std::complex<double> *bra, const std::complex<double> *cpl, const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets the chi2 for a single t' bin
template<typename xdouble>
xdouble anchor_t::EvalTbin(
							int 								tbin,
							const std::complex<xdouble>	 				*cpl,
							const xdouble	 						*par,
							const xdouble 							*iso_par)		const{

	xdouble chi2 = 0.;
	std::vector<std::vector<std::complex<xdouble> > > iso_eval = _waveset.iso_funcs(iso_par);
	for (size_t bin=_waveset.minBin(); bin<_waveset.maxBin(); bin++){
		chi2+=EvalBin(tbin,bin,cpl,par,iso_eval);
//		std::cout<<"bin #"<<bin<<" chi2++"<<EvalBin(tbin,bin,cpl,par,iso_eval)<<std::endl;
//		std::cout << bin<<"   "<< EvalBin(tbin,bin,cpl,par) <<std::endl;
	};
	return chi2;
};
template double anchor_t::EvalTbin(int tbin, const std::complex<double> *cpl,const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets the chi2 for a certain t' bin, with automatically calculated couplings
template<typename xdouble>
xdouble anchor_t::EvalAutoCplTbin(
							int 								tbin,
							const std::complex<xdouble>	 				*cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{

	std::vector<std::complex<xdouble> > mincpl = getMinimumCpl(tbin,cpl,par,iso_par);
//	std::cout<<tbin<<": ";
//	print_vector(mincpl);
	xdouble chi2 = EvalTbin(tbin,&mincpl[0],par,iso_par); // &...[0] stays, since mincpl is built inside the function
//	std::cout<<chi2<<std::endl;
	return chi2;
};
template double anchor_t::EvalAutoCplTbin(int tbin, const std::complex<double> *cpl, const double *par, const double *iso_par)const;
//########################################################################################################################################################
///Gets the Chi2 for a single t' and m3pi bin
template<typename xdouble>
xdouble anchor_t::EvalBin(
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
	for (size_t i=0;i<2*_waveset.nPoints()-1;i++){
		int iWave = (*_waveset.point_to_wave())[(i+1)/2];
//		if (mass >= _lowerLims[iWave] and mass < _upperLims[iWave]){
		if (_is_active[tbin][bin][iWave]){
//			std::cout<<"le_addite:D"<<pow(deltas[i],2.)*_coma[tbin][bin][i][i]<<std::endl;
			if (_waveset.write_out()){
				*_waveset.outStream() <<" mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<i+1<<"  ipj=            "<<i+1<<"  isectd=           "<<tbin+1<<std::endl;
				*_waveset.outStream() <<"  delta1(imb, ipi, isectd) delta1(imb, ipj, isectd)     "<<deltas[i]<<"        "<<deltas[i]<<std::endl;
				*_waveset.outStream() <<"   cov_key(imb, ipi, ipj,  isectd) =     "<<_coma[tbin][bin][i][i]<<std::endl;
				*_waveset.outStream() <<"  Ci^2 before      "<<chi2<<"      +     "<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<"       =    ";
			};
			chi2+= deltas[i]*deltas[i]*_coma[tbin][bin][i][i];


//			std::cout<<deltas[i]<<"*"<<deltas[i]<<"*"<<_coma[tbin][bin][i][i]<<"="<<deltas[i]*deltas[i]*_coma[tbin][bin][i][i]<<std::endl;
			if(_waveset.write_out()){
				*_waveset.outStream() <<chi2<<std::endl;
			};
			for (size_t j=0;j<i;j++){
				int jWave = (*_waveset.point_to_wave())[(j+1)/2];
//				if(mass >= _lowerLims[jWave] and mass < _upperLims[jWave]){
				if(_is_active[tbin][bin][jWave]){
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
template double anchor_t::EvalBin(int tbin,int bin,const std::complex<double> *cpl,const double *par,std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
///Returns f(m,...) - data[...] for each SDM entry in the fit
template<typename xdouble>
std::vector<xdouble> anchor_t::delta(
							int 								tbin,
							int 								bin,
							double 								mass,
							const std::complex<xdouble> 					*cpl,
							const xdouble	 						*par,
							std::vector<std::vector<std::complex<xdouble> > > 		&iso_eval)		const{

	std::vector<double> var = _waveset.getVar(mass,tbin);
	std::vector<std::complex<xdouble> > ampls = _waveset.amps(&var[0], cpl, par, iso_eval);

///	std::cout<<"''''''''''''''''''''''''''''''''''''''''"<<std::endl;
///	for (size_t i=0;i<ampls.size();++i){
///		std::cout<<i<<": "<<ampls[i]<<std::endl;
///	};

	std::vector<xdouble> del = std::vector<xdouble> (2*_waveset.nPoints()-1);
	std::complex<xdouble> divider = std::complex<xdouble>(1.,0.); // When amplitudes are fitted, this is set to |ampls[0]|
	if(_is_ampl){
		divider = std::complex<xdouble>(pow(std::norm(ampls[0]),.5),0.);
	};
	if(_waveset.write_out()){
		*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           1  isectd=           "<<tbin+1<<std::endl;
		*_waveset.outStream() << " Re1 data    "<<_data[tbin][bin][0]<<"     - theory   "<<std::norm(ampls[0]/divider)<<"       =    "<<std::norm(ampls[0]/divider) - _data[tbin][bin][0]<<std::endl;
	};
	del[0]=std::norm(ampls[0]/divider) - _data[tbin][bin][0]; ////// HERE ADD if(_is_point_bin[bin][i]){...}; so that no deltas are aclculated for turend off bins
	for (size_t i = 1; i<_waveset.nPoints();i++){
		if (not _is_active[tbin][bin][i]){
			continue;
		};
		std::complex<xdouble> inter = ampls[0]*std::conj(ampls[i])/divider;
		del[2*i-1]=real(inter) - _data[tbin][bin][2*i-1]; // real part
		del[2*i  ]=imag(inter) - _data[tbin][bin][2*i];    // imag part

		if(_waveset.write_out()){
			*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<2*i<<"  isectd=           "<<tbin+1<<std::endl;
			*_waveset.outStream() << " Re data    "<<_data[tbin][bin][2*i-1]<<"     - theory   "<<real(inter)<<"       =    "<<real(inter) - _data[tbin][bin][2*i-1]<<std::endl;
			*_waveset.outStream() << " mass   "<<mass<<"      imb=           "<<bin<<"  ipi=           "<<2*i+1<<"  isectd=           "<<tbin+1<<std::endl;
			*_waveset.outStream() << " Im data    "<<_data[tbin][bin][2*i]<<"     - theory   "<<imag(inter)<<"       =    "<<imag(inter) - _data[tbin][bin][2*i]<<std::endl;
		};
	};
	return del;
};
template std::vector<double> anchor_t::delta(int tbin, int bin,double mass, const std::complex<double> *cpl, const double *par, std::vector<std::vector<std::complex<double> > > &iso_eval) const;
//########################################################################################################################################################
/*
		┌───────────────────────────────────────┐
		│	Es folgt eine Indexschlacht	│
		└───────────────────────────────────────┘
*/
//########################################################################################################################################################
///Gets the coefficients A and B for Chi2 = x^T*A*x + B*x + C, works with isobars
template<typename xdouble>
AandB<xdouble> anchor_t::get_AB(
							int 								tbin,
							const std::complex<xdouble>	 				*anchor_cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)		const{



	int nCplAnc = (*_waveset.borders_waves())[0]; // Number of couplings for the anchor wave
	std::vector<std::vector<std::complex<xdouble> > > iso_eval  =  _waveset.iso_funcs(iso_par);					//// ><>< Maybe just give pointer to already evaluated functions to save time???
	int nNon = _waveset.nFtw() - nCplAnc;
	AandB<xdouble> AB(2*nNon);
	std::vector<xdouble> lines = std::vector<xdouble>(2*_waveset.nPoints()-2,0.); 							//// >>>>  _waveset.nPoints()
	for(size_t bin=_waveset.minBin();bin<_waveset.maxBin();bin++){
		std::vector<double> var = _waveset.getVar(_waveset.get_m(bin),tbin);
		std::vector<std::complex<xdouble> > func = _waveset.funcs(&var[0],par);

		std::vector<double> phase = _waveset.phase_space(&var[0]);
		std::complex<xdouble> ampAnc = std::complex<xdouble>(0.,0.);
		for(int i=0;i<nCplAnc;i++){
			ampAnc+= anchor_cpl[i]*func[(*_waveset.funcs_to_waves())[i]]*phase[0];
		};
		if(_is_ampl){ // Divide ampAnc by |ampAnc| to get the amplitude with phase 0. // Need to be tested
			ampAnc/=pow(std::norm(ampAnc),.5);
		};
		if(ampAnc == std::complex<xdouble>(0.,0.)){
			continue;
		};
		xdouble RR = ampAnc.real();
		xdouble II = ampAnc.imag();
		for (size_t i =0;i<2*_waveset.nPoints()-2;i++){ 										//// >>>>  _waveset.nPoints()
			xdouble val = (RR*RR+II*II - _data[tbin][bin][0])*_coma[tbin][bin][i+1][0];
			for (size_t j=1;j<2*_waveset.nPoints()-1;j++){ 									//// >>>>  _waveset.nPoints()
				val+=-_data[tbin][bin][j]*_coma[tbin][bin][i+1][j];
			};
			lines[i]=val;
		};
		size_t wave1 = 1;	// Start with 1, leave anchor wave out								//// |||| For de-isobarred, this stays, because anchor wave must be isobarred
		size_t up1 = (*_waveset.borders_waves())[1];
		for (size_t nf1 = nCplAnc;nf1<_waveset.nFtw();nf1++){
			if (nf1==up1){
				wave1+=1;
				up1 = (*_waveset.borders_waves())[wave1];
			};
			int f1 = (*_waveset.funcs_to_waves())[nf1];
			int iso_f1 = (*_waveset.iso_to_waves())[nf1];									//// >>>> Int f isobar 1
			int iso_nBin1 = std::abs((*_waveset.wave_binning_pts())[wave1]);
			int ind1=nf1-nCplAnc;
			int wave_start1 = (*_waveset.point_borders_wave())[wave1-1];								//// >>>> To calculate the _data and _coma indices

			for (int i_iso_1=0;i_iso_1<iso_nBin1;i_iso_1++){	//// >>>> Here loop ofer isobar bins, if necessary (i_iso_1)(_wave_binning_pts[wave1)
				int point1 = wave_start1 + i_iso_1;	// Position of data and coma points
				if (not _is_active[tbin][bin][point1]){
					continue;
				};
				std::complex<xdouble> AB1;
				if (iso_f1 == -1){
					AB1 = ampAnc * conj(func[f1]) * phase[wave1];
				}else{
					AB1 = ampAnc * conj(func[f1]*iso_eval[iso_f1][i_iso_1]) * phase[wave1];				//// >>>> Here take isobar function evaluation
				};
				if (AB1 == std::complex<xdouble>(0.,0.)){
					continue;
				};

				size_t wave2=1;
				size_t up2=(*_waveset.borders_waves())[1];
				//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)

				AB.B[2*ind1  ] += AB1.real() * lines[2*point1-2] * 2;							//// >>>> now is Point
				AB.B[2*ind1  ] += AB1.imag() * lines[2*point1-1] * 2;							//// >>>> now is Point

				AB.B[2*ind1+1] +=-AB1.real() * lines[2*point1-1] * 2;							//// >>>> now is Point
				AB.B[2*ind1+1] += AB1.imag() * lines[2*point1-2] * 2;							//// >>>> now is Point

				for (size_t nf2= nCplAnc;nf2<_waveset.nFtw();nf2++){
					if(nf2==up2){
						wave2+=1;
						up2 = (*_waveset.borders_waves())[wave2];
					};
					int f2 = (*_waveset.funcs_to_waves())[nf2];
					int iso_f2 = (*_waveset.iso_to_waves())[nf2];								//// >>>> Int f isobar 2
					int iso_nBin2 = abs((*_waveset.wave_binning_pts())[wave2]);
					int wave_start2 = (*_waveset.point_borders_wave())[wave2-1];							//// >>>> To calculate the _data and _coma indices
					int ind2=nf2-nCplAnc;										//// |||| ind still applies!!!
					for (int i_iso_2=0;i_iso_2<iso_nBin2;i_iso_2++){	//// >>>> loop over second isobar bins (i_iso_2) (_wave_binning_pts[wave2)

						int point2 = wave_start2+i_iso_2;
						if (not _is_active[tbin][bin][point2]){
							continue;
						};
						std::complex<xdouble> AB2;
//						xdouble r2;
//						xdouble i2;
						if (iso_f2 == -1){
							AB2 = ampAnc * conj(func[f2]) * phase[wave2];
						}else{
							AB2 = ampAnc * conj(func[f2]*iso_eval[iso_f2][i_iso_2]) * phase[wave2];				//// >>>> Here take isobar function evaluation
						};
				

						//// (RR + j II)*(R - j I)*(r - j i) = (RR R r - RR I i - II I r + II R i)+j(II R r + II I i - RR r I - RR R i)

						AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1  ][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1  ][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2  ]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.imag();	//// >>>> are points1/2 now

						AB.A[2*ind1+1][2*ind2+1]+= AB1.real() * _coma[tbin][bin][2*point1  ][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+=-AB1.real() * _coma[tbin][bin][2*point1  ][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+=-AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2  ] * AB2.real();	//// >>>> are points1/2 now
						AB.A[2*ind1+1][2*ind2+1]+= AB1.imag() * _coma[tbin][bin][2*point1-1][2*point2-1] * AB2.imag();	//// >>>> are points1/2 now

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now

						std::cout << AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.real() <<" "<< _coma[tbin][bin][2*point1  ][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
						std::cout <<-AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2  ] <<" "<< AB2.real() <<std::endl;	//// >>>> are points1/2 now
						std::cout << AB1.imag() <<" "<< _coma[tbin][bin][2*point1-1][2*point2-1] <<" "<< AB2.imag() <<std::endl;	//// >>>> are points1/2 now
*/
					};
				};
			};
		};
	};
	return AB;
};
template AandB<double> anchor_t::get_AB(int tbin,const std::complex<double> *anchor_cpl, const double *par, const double *iso_par) const;
//########################################################################################################################################################
///Gets the best couplings without branchings
template<typename xdouble>
std::vector<std::complex<xdouble> > anchor_t::getMinimumCpl(
							int 								tbin,
							const std::complex<xdouble>					*anchor_cpl,
							const xdouble 							*par,
							const xdouble 							*iso_par)		const{
	int nCplAnc = (*_waveset.borders_waves())[0]; // Number of couplings for the anchor wave
	int nNon = _waveset.nFtw() - nCplAnc;
	std::vector<std::complex<xdouble> > cpl = std::vector<std::complex<xdouble> >(_waveset.nFtw());
	for (int i=0; i<nCplAnc; i++){
		cpl[i] = anchor_cpl[i];
	};
	AandB<xdouble>AB = get_AB(tbin,anchor_cpl,par,iso_par);

	std::vector<xdouble> D = cholesky::cholesky_solve(AB.A,AB.B);

	for (int i=0;i<nNon;i++){
		// Short explanation:
		// ┌────────────────────────────────────────────────────────────────────┐
		// │C^T * A * C + B * C + const. = Chi2 is a quadratic form		│
		// │To find the minimum: dChi2/dC_i = 0. for all i in {0,...,2*nNon-1}	│
		// │=> (A*C)_i + (C^T*A)_i + B_i = 0.					│
		// │Since in the Chi2, the assymmetric part of A cancels out anyway:	│
		// │=> 2*A*C = -B (**)							│
		// │D solves A*D = B (*)						│
		// │D = -2*C								│
		// │=> C solves (**)							│
		// └────────────────────────────────────────────────────────────────────┘
		cpl[nCplAnc+i]=std::complex<xdouble>(-D[2*i]/2.,-D[2*i+1]/2.);
	};
	return cpl;
};
template std::vector<std::complex<double> > anchor_t::getMinimumCpl(int tbin,const std::complex<double> *anchor_cpl, const double *par, const double *iso_par) const;
//########################################################################################################################################################
///Calculated all non anchor couplings (no need for fitting)
template<typename xdouble>
std::vector<std::complex<xdouble> > anchor_t::getMinimumCplBra(
							int 								tbin,
							const std::complex<xdouble>	 				*branch,
							const std::complex<xdouble>	 				*anchor_cpl,
							const xdouble	 						*par,
							const xdouble	 						*iso_par)			const{

	if (0==_waveset.nBranch()){ // Do not do the complicated stuff, when no branchings are used
		return getMinimumCpl(tbin,anchor_cpl,par, iso_par);
	}else{
		int nCplAnc = (*_waveset.borders_waves())[0];
		int nNon = _waveset.nFtw() - nCplAnc;
		std::vector<std::complex<xdouble> > cplAncBr = std::vector<std::complex<xdouble> >(nCplAnc);
		std::vector<std::complex<xdouble> > ccppll(_waveset.nBrCpl());
		for (size_t i=0;i<_nBrCplAnc;i++){
			ccppll[i] = anchor_cpl[i];
		};
		for (int i=0;i<nCplAnc;i++){
			int ncpl_act = (*_waveset.n_cpls())[i];
			int nbranch_act = (*_waveset.n_branch())[i];
			if (-1 == nbranch_act){
				cplAncBr[i] = anchor_cpl[ncpl_act];
			}else{
				cplAncBr[i] = anchor_cpl[ncpl_act] * branch[nbranch_act];
			};
		};
		AandB<xdouble>AB = get_AB(tbin,&cplAncBr[0],par,iso_par); // &...[0] stays here, since the vector is built inside the function
		AandB<xdouble> ABprime(2*(_waveset.nBrCpl() - _nBrCplAnc));
		for (int i =0;i<nNon;i++){ // Reshuffle the coefficients
			int i_tot = i+nCplAnc;
			int i_cpl_tot = (*_waveset.n_cpls())[i_tot];
			int i_cpl = i_cpl_tot - _nBrCplAnc;
			int i_br  = (*_waveset.n_branch())[i_tot];
	//		std::cout<<"────────────────────────────────────────────────"<<std::endl;
	//		std::cout<<i<<" "<<i_tot<<" "<<i_cpl_tot<<" "<<i_cpl<<" "<<i_br<<std::endl;
			if (-1 == i_br){ // If the funcion is uncoupled, just take the old coefficients
				ABprime.B[2*i_cpl  ]+= AB.B[2*i  ];
				ABprime.B[2*i_cpl+1]+= AB.B[2*i+1];
			}else if(i_cpl>-1){ // If the function is not coupled to the anchor wave, add it to the corresponding coupling
				ABprime.B[2*i_cpl  ]+= AB.B[2*i  ]*branch[i_br].real()   +   AB.B[2*i+1]*branch[i_br].imag();
				ABprime.B[2*i_cpl+1]+=-AB.B[2*i  ]*branch[i_br].imag()   +   AB.B[2*i+1]*branch[i_br].real();
			};
			for (int j=0; j<nNon;j++){
				int j_tot = j+nCplAnc;
				int j_cpl_tot = (*_waveset.n_cpls())[j_tot];
				int j_cpl = j_cpl_tot -_nBrCplAnc;
				int j_br  = (*_waveset.n_branch())[j_tot];
				if ((-1 == i_br or i_cpl>-1)and(-1 == j_br or j_cpl>-1)){
					std::complex<xdouble> ifak(1.,0.);
					std::complex<xdouble> jfak(1.,0.);
					if (i_br >-1){
						ifak = branch[i_br];
					};
					if (j_br >-1){
						jfak = branch[j_br];
					};
					ABprime.A[2*i_cpl  ][2*j_cpl  ]+=	 AB.A[2*i  ][2*j  ]*ifak.real()*jfak.real() +
										 AB.A[2*i  ][2*j+1]*ifak.real()*jfak.imag() +
										 AB.A[2*i+1][2*j  ]*ifak.imag()*jfak.real() +
										 AB.A[2*i+1][2*j+1]*ifak.imag()*jfak.imag() ;
					ABprime.A[2*i_cpl  ][2*j_cpl+1]+=	-AB.A[2*i  ][2*j  ]*ifak.real()*jfak.imag() +
										 AB.A[2*i  ][2*j+1]*ifak.real()*jfak.real() +
										-AB.A[2*i+1][2*j  ]*ifak.imag()*jfak.imag() +
										 AB.A[2*i+1][2*j+1]*ifak.imag()*jfak.real() ;
					ABprime.A[2*i_cpl+1][2*j_cpl  ]+=	-AB.A[2*i  ][2*j  ]*ifak.imag()*jfak.real() +
										-AB.A[2*i  ][2*j+1]*ifak.imag()*jfak.imag() +
										 AB.A[2*i+1][2*j  ]*ifak.real()*jfak.real() +
										 AB.A[2*i+1][2*j+1]*ifak.real()*jfak.imag() ;
					ABprime.A[2*i_cpl+1][2*j_cpl+1]+=	 AB.A[2*i  ][2*j  ]*ifak.imag()*jfak.imag() +
										-AB.A[2*i  ][2*j+1]*ifak.imag()*jfak.real() +
										-AB.A[2*i+1][2*j  ]*ifak.real()*jfak.imag() +
										 AB.A[2*i+1][2*j+1]*ifak.real()*jfak.real() ;
	//				std::cout<<ifak<<jfak<<std::endl;
				}else if(-1 == i_br or i_cpl>-1){
					std::complex<xdouble> ifak(1.,0.);
					if (i_br>-1){
						ifak = branch[i_br];
					};
					std::complex<xdouble> cpl_j = branch[j_br]*anchor_cpl[j_cpl_tot];
					ABprime.B[2*i_cpl  ]+=   AB.A[2*i  ][2*j  ]*ifak.real()*cpl_j.real() +
								 AB.A[2*i  ][2*j+1]*ifak.real()*cpl_j.imag() +
								 AB.A[2*i+1][2*j  ]*ifak.imag()*cpl_j.real() +
								 AB.A[2*i+1][2*j+1]*ifak.imag()*cpl_j.imag() ;
					ABprime.B[2*i_cpl+1]+=	-AB.A[2*i  ][2*j  ]*ifak.imag()*cpl_j.real() +
								-AB.A[2*i  ][2*j+1]*ifak.imag()*cpl_j.imag() +
								 AB.A[2*i+1][2*j  ]*ifak.real()*cpl_j.real() +
								 AB.A[2*i+1][2*j+1]*ifak.real()*cpl_j.imag() ;
	//				std::cout<<"Nirrr"<<std::endl;
				}else if(-1 == j_br or j_cpl>-1){
					std::complex<xdouble> jfak(1.,0.);
					if (j_br>-1){
						jfak = branch[j_br];
					};
					std::complex<xdouble> cpl_i = branch[i_br]*anchor_cpl[i_cpl_tot];
					ABprime.B[2*j_cpl  ]+=	 AB.A[2*i  ][2*j  ]*cpl_i.real()*jfak.real() +
								 AB.A[2*i  ][2*j+1]*cpl_i.real()*jfak.imag() +
								 AB.A[2*i+1][2*j  ]*cpl_i.imag()*jfak.real() +
								 AB.A[2*i+1][2*j+1]*cpl_i.imag()*jfak.imag() ;
					ABprime.B[2*j_cpl+1]+=	-AB.A[2*i  ][2*j  ]*cpl_i.real()*jfak.imag() +
								 AB.A[2*i  ][2*j+1]*cpl_i.real()*jfak.real() +
								-AB.A[2*i+1][2*j  ]*cpl_i.imag()*jfak.imag() +
								 AB.A[2*i+1][2*j+1]*cpl_i.imag()*jfak.real() ;
	//				std::cout<<"Njrrr"<<std::endl;
				};// If both functions are coupled to the anchor wave, they just add to the constant C, (Chi2 = x^TAx + B^Tx +C) and do not change the position of the minimum.
			};
		};
	//	std::vector<double> les_as;
	//	for (int i=0;i<2*(_waveset.nBrCpl() - _nBrCplAnc);i++){
	//		les_as.push_back(ABprime.A[i][i]);
	//	};
	//	print_vector(les_as);
		std::vector<xdouble> D = cholesky::cholesky_solve(ABprime.A,ABprime.B);
		for (size_t i=_nBrCplAnc;i<_waveset.nBrCpl();i++){
			int irel = i-_nBrCplAnc;
			ccppll[i] = std::complex<xdouble>(-D[2*irel]/2,-D[2*irel+1]/2);
		};
	//	std::cout<<"made it through"<<std::endl;
		return ccppll;
	};
};
template std::vector<std::complex<double> > anchor_t::getMinimumCplBra(int tbin, const std::complex<double> *branch, const std::complex<double> *anchor_cpl, const double *par, const double *iso_par) const;
//########################################################################################################################################################
///Instantiate auto diff methods, if needed (Enable adouble operations, if the auto diff package is loaded)
#ifdef ADOL_ON
template adouble anchor_t::EvalCP(const std::complex<adouble> *cpl,const adouble *par, const adouble *iso_par) const;
template adouble anchor_t::EvalBranch(const std::complex<adouble> *branch, const std::complex<adouble> *cpl, const adouble *par, const adouble *iso_par) const;
template adouble anchor_t::EvalTbin(int tbin, const std::complex<adouble> *cpl,const adouble *par, const adouble *iso_par) const;
template adouble anchor_t::EvalBin(int tbin,int bin,const std::complex<adouble> *cpl,const adouble *par, std::vector<std::vector<std::complex<adouble> > > &iso_par) const;
template std::vector<adouble> anchor_t::delta(int tbin, int bin,double mass, const std::complex<adouble> *cpl, const adouble *par, std::vector<std::vector<std::complex<adouble> > > &iso_eval) const;
template adouble anchor_t::EvalAutoCpl(const std::complex<adouble> *cpl,const adouble *par, const adouble *iso_par)const;
template adouble anchor_t::EvalAutoCplBranch(const std::complex<adouble> *bra, const std::complex<adouble> *cpl, const adouble *par, const adouble *iso_par)const;
template adouble anchor_t::EvalAutoCplTbin(int tbin, const std::complex<adouble> *cpl, const adouble *par, const adouble *iso_par)const;
template AandB<adouble> anchor_t::get_AB(int tbin,const std::complex<adouble> *anchor_cpl, const adouble *par, const adouble *iso_par) const;
template std::vector<std::complex<adouble> > anchor_t::getMinimumCpl(int tbin,const std::complex<adouble> *anchor_cpl, const adouble *par, const adouble *iso_par) const;
template std::vector<std::complex<adouble> > anchor_t::getMinimumCplBra(int tbin, const std::complex<adouble> *branch, const std::complex<adouble> *anchor_cpl, const adouble *par, const adouble *iso_par) const;
//#######################################################################################################################################################
///Gets the gradient w.r.t. xx
std::vector<double> anchor_t::Diff(
							std::vector<double> 				&xx)						const{

	return Diff(&xx[0]);
};
//#######################################################################################################################################################
///Gets the gradient w.r.t. xx
std::vector<double> anchor_t::Diff(
							const double 					*xx)						const{

	int nTape = 0;
	double x[_nTot];
	for (size_t i=0;i<_nTot;i++){
		x[i] = xx[i];
	};

	trace_on(nTape);
	std::vector<adouble>aCpl_r(2*_nCpl);
	std::vector<adouble>aPar(_nPar);
	std::vector<adouble>aBra_r(2*_nBra);
	std::vector<adouble>aIso(_nIso);
	int count=0;
	for (size_t i=0;i<2*_nCpl;i++){
		aCpl_r[i] <<= x[count];
		count++;
	};
	for (size_t i=0; i<_nPar;i++){
		aPar[i] <<= x[count];
		count++;
	};
	for (size_t i=0;i<2*_nBra;i++){
		aBra_r[i] <<= x[count];
		count++;
	};
	for (size_t i=0;i<_nIso;i++){
		aIso[i] <<= x[count];
		count++;
	};
	std::vector<std::complex<adouble> > aCpl_c(_nCpl);
	std::vector<std::complex<adouble> > aBra_c(_nBra);
	for (size_t i=0;i<_nCpl;i++){
		aCpl_c[i] = std::complex<adouble>(aCpl_r[2*i],aCpl_r[2*i+1]);
	};
	for (size_t i=0;i<_nBra;i++){
		aBra_c[i] = std::complex<adouble>(aBra_r[2*i],aBra_r[2*i+1]);
	};
	double Chi2;
	adouble aChi2;
	aChi2 = EvalAutoCplBranch(&aBra_c[0],&aCpl_c[0],&aPar[0],&aIso[0]);//[0]//
	aChi2 >>= Chi2;
	trace_off();
	double grad[_nTot];
	gradient(nTape,_nTot,x,grad);
	vector<double> gradient(_nTot);
	for (size_t i=0;i<_nTot;i++){
		gradient[i]=grad[i];
	};
	return gradient;
};
#endif//ADOL_ON
//########################################################################################################################################################
///Calculates some initial values for the branchigs
std::vector<std::complex<double> > anchor_t::get_branchings(
							const std::vector<std::complex<double> > 	&cpl,
							const std::vector<double> 			&par,
							const std::vector<double> 			&iso_par)					const{

	// Gets branchings as follows:
	//	- Take Anchor cpls to get couplings for each t-bin
	//	- Average over the tbins for all couplings
	//	- Determine the branchings from these averages
	std::vector<std::complex<double> > brchs = std::vector<std::complex<double> >(_waveset.nBranch());
	std::vector<std::complex<double> > all_cpls = std::vector<std::complex<double> >(_waveset.nFtw());
	int active_t_bins=0;
	int cpl_used=0;
	for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){ // Average couplings over all t' bins, after rotating to phase of cpl[0] to 0
		if ((*_waveset.eval_tbin())[tbin]){
			std::vector<std::complex<double> > cpl_t = std::vector<std::complex<double> >(_nBrCplAnc);
			for (size_t i=0;i<_nBrCplAnc;i++){
				cpl_t[i] = cpl[cpl_used];
				cpl_used++;
			};
			active_t_bins++;
			std::vector<std::complex<double> > all_cpls_t= getMinimumCpl(tbin,&cpl_t[0],&par[0], &iso_par[0]);//[0]//
			std::complex<double> phase = all_cpls_t[0]/abs(all_cpls_t[0]);
			for (size_t i=0;i<_waveset.nFtw();i++){
				all_cpls[i]+=all_cpls_t[i]/phase;
			};
		};
	};
	for (size_t i=0;i<_waveset.nFtw();i++){ // Divide by number of t' bins
		all_cpls[i]/=active_t_bins;
	};
	for (size_t j=0;j<(*_waveset.n_branch()).size();j++){ // Write coupling-average to corresponding branching
		if ((*_waveset.n_branch())[j] >= 0){
			brchs[(*_waveset.n_branch())[j]] = all_cpls[j];
		};
	};
	return brchs;
};
//########################################################################################################################################################
///Gets all couplings for one t' bin (.size() = _waveset.nFtw()) for different input formats
std::vector<std::complex<double> > anchor_t::getAllCouplings(
							int						tbin,
							const std::vector<std::complex<double> >	&cpl,
							const std::vector<double>			&par,
							const std::vector<std::complex<double> >	&bra,
							const std::vector<double>			&iso)						const{

	std::vector<std::complex<double> > cpl_all;
	if (cpl.size() == _nCpl){ 				// Anchor couplings for all t' bins
		std::vector<std::complex<double> > cpl_t;
		for (size_t i=0;i<_nBrCplAnc;i++){
			cpl_t.push_back(cpl[tbin*_nBrCplAnc+i]);
		};
		cpl_all = getMinimumCplBra(tbin,&bra[0],&cpl_t[0],&par[0],&iso[0]);//[0]//
		cpl_all = getUnbranchedCouplings(cpl_all,bra);
//		std::cout<<"getAllCouplings(...): Take couplings as anchor couplings for all t' bins"<<std::endl;
	}else if (cpl.size() == _nBrCplAnc){ 			// Anchor couplings for one t' bin
		cpl_all = getMinimumCplBra(tbin,&bra[0],&cpl[0],&par[0],&iso[0]);//[0]//
		cpl_all = getUnbranchedCouplings(cpl_all,bra);
//		std::cout<<"getAllCouplings(...): Take couplings as anchor couplings for one t' bin"<<std::endl;
	}else if (cpl.size() == _waveset.nBrCpl()*_waveset.nTbin()){ 		// Branched couplings for all t' bins
		std::vector<std::complex<double> > cpl_t;
		for (size_t i=0;i<_waveset.nBrCpl();i++){
			cpl_t.push_back(cpl[tbin*_waveset.nBrCpl()+i]);
		};
		cpl_all = getUnbranchedCouplings(cpl_t,bra);
//		std::cout<<"getAllCouplings(...): Take couplings as branched couplings for all t' bins"<<std::endl;
	}else if (cpl.size() == _waveset.nBrCpl()){ 			// Branched couplings for one t' bin
		cpl_all = getUnbranchedCouplings(cpl,bra);
//		std::cout<<"getAllCouplings(...): Take couplings as branched couplings for one t' bin"<<std::endl;
	}else if (cpl.size() == _waveset.nTbin()*_waveset.nFtw()){ 			// Simple couplings for all t' bins
		for (size_t i=0;i<_waveset.nFtw();i++){
			cpl_all.push_back(cpl[tbin*_waveset.nFtw()+i]);
		};
//		std::cout<<"getAllCouplings(...): Take couplings as simple couplings for all t' bins"<<std::endl;
	}else if (cpl.size() == _waveset.nFtw()){ 			// Simple couplings for one t' bin
		cpl_all = cpl;
//		std::cout<<"getAllCouplings(...): Take couplings as simple couplings for one t' bin"<<std::endl;
	}else{
		std::cerr<<"anchor_t::getAllCouplings(...): Error: Can't determine the format of the couplings"<<std::endl;
	};
	return cpl_all;
};
//########################################################################################################################################################
///Set the anchor coupligs to one, that are coupled to another function
void anchor_t::branchCouplingsToOne(){
	for (size_t i=0;i<(*_waveset.coupled()).size();i++){
		if((*_waveset.coupled())[i] >=0){
			size_t nCpl = (*_waveset.n_cpls())[i];
			if (nCpl <= _nBrCplAnc){
				for (size_t tbin =0;tbin<_waveset.nTbin();tbin++){
					setParameter(2*_nBrCplAnc*tbin+2*nCpl  ,1.);
					setParameter(2*_nBrCplAnc*tbin+2*nCpl+1,0.);
				};
			};
		};
	};
};
//########################################################################################################################################################
///Set data points manually
bool anchor_t::set_data(
							int								tbin,
							int 								bin,
							std::vector<double> 						data){

	if (_data.size() != _waveset.nTbin()){
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_data[tbin].size() != _waveset.nBins()){
		_data[tbin] = std::vector<std::vector<double> >(_waveset.nBins());
	};
	_data[tbin][bin] = data;
	update_is_active();

	if (data.size() == 2*_waveset.nPoints() -1){
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Set coma manually
bool anchor_t::set_coma(
							int 								tbin,
							int 								bin,
							std::vector<std::vector<double> > 				coma){

	//print_matrix(coma);
	if(_coma.size() != _waveset.nTbin()){
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
	if(_coma[tbin].size() != _waveset.nBins()){
		_coma[tbin] = std::vector<std::vector<std::vector<double> > >(_waveset.nBins());
	};
	_coma[tbin][bin] = coma;
	if (coma.size() == 2*_waveset.nPoints()-1){
		for (size_t i=0;i<coma.size();i++){
			if (coma[i].size() != coma.size()){
				std::cout<<"anchor_t::set_coma(...): Warning: Set coma is not quadratic."<<std::endl;
				return false;
			};
		};
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Loads the data for a t' bin from a file
void anchor_t::loadData(
							int 								tbin,
							const char* 							dataFile){

	_data[tbin] = std::vector<std::vector<double> >();
	std::fstream data(dataFile,std::ios_base::in);
	double val;
	size_t nDat = 2*_waveset.nPoints()-1;
	std::vector<double> data_bin;
	while(data >> val){
		data_bin.push_back(val);
		if (data_bin.size() == nDat){
			_data[tbin].push_back(data_bin);
			data_bin = std::vector<double>();
		};
	};
	update_is_active();
	if (_waveset.nBins() != _data[tbin].size()){
		std::cout << "anchor_t::loadData(...): Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _data.size()="<<_data[tbin].size()<<std::endl;
	}else{
		std::cout << "anchor_t::loadData(...): File delivered the right size for _data"<<std::endl;
	};
};
//########################################################################################################################################################
///Loads the (inverted) covariance matrix for a t' bin from 'comaFile'
void anchor_t::loadComa(
							int 								tbin,
							const char* 							comaFile){

	_coma[tbin]=std::vector<std::vector<std::vector<double> > >();
	std::fstream data(comaFile,std::ios_base::in);
	double val;
	size_t nDat = 2*_waveset.nPoints()-1;
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
		std::cout << "anchor_t::loadComa(...): Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _coma.size()="<<_coma[tbin].size() << std::endl;
	}else{
		std::cout<< "anchor_t::loadComa(...): File delivered the right size for _coma"<<std::endl;
	};
};
//########################################################################################################################################################
///Returns the number of anchor couplings
int anchor_t::getNanc(){
	return _nBrCplAnc * _waveset.nTbin();
};
//########################################################################################################################################################
///Switches on/off usage of branchings in the operator()
bool anchor_t::setUseBranch(bool in){
	_useBranch = in;
	return true;
};
//########################################################################################################################################################
///Gets the total number of paramters with only anchor couplings
int anchor_t::getNtotAnc(){

	return 2*getNanc() + _waveset.getNpar() + 2*_waveset.nBranch() + _waveset.getNiso();
};
//########################################################################################################################################################
///Sets all data to 0.
void anchor_t::nullify(){

	for (size_t tbin = 0; tbin < _waveset.nTbin();tbin++){
		for (size_t bin =0; bin<_waveset.nBins();bin++){
			int nDat = 2*_waveset.nPoints()-1;
			for (int i=0;i<nDat;i++){
				_data[tbin][bin][i]=0.;
			};
		};
	};
};
//########################################################################################################################################################
///Updates the number of couplings
void anchor_t::update_n_cpls(){

	_waveset.update_n_cpls();
	_nBrCplAnc = (*_waveset.n_cpls())[(*_waveset.borders_waves())[0]-1]+1;
};
//########################################################################################################################################################
/// Prints the internal status
void anchor_t::printStatus()																const{
	_waveset.printStatus();
	std::cout<<std::endl<<"ANCHOR:"<<std::endl;
	std::cout<<std::endl<<"PARAMETER NUMBERS:"<<std::endl;
	std::cout<<"_nTot: "<<_nTot<<"; _nPar: "<<_nPar<<"; _nCpl: "<<_nCpl<<std::endl;
	std::cout<<"_nBra: "<<_nBra<<"; _nIso: "<<_nIso<< "; _nBrCplAnc: "<<_nBrCplAnc<<std::endl;
	std::cout<<std::endl<<std::endl<<"PARAMETERS:"<<std::endl<<"_parameters"<<std::endl;
	print_vector(_parameters);
	std::cout<<std::endl<<"_upper_parameter_limits"<<std::endl;
	print_vector(_upper_parameter_limits);
	std::cout<<std::endl<<"_lower_parameter_limits"<<std::endl;
	print_vector(_lower_parameter_limits);
	std::cout<<std::endl<<"_parNames"<<std::endl;
	print_vector(_parNames);
	std::cout<<std::endl<<"Omit printing of _data and _coma"<<std::endl;
	std::cout<<std::endl<<"OTHER MEMBERS:"<<std::endl;
	std::cout<<"_is_ampl: "<<_is_ampl<<std::endl;
	std::cout<<std::endl<<"_useBranch: "<<_useBranch<<std::endl;
	std::cout<<std::endl<<"_nOut: "<<_nOut<<std::endl;
	std::cout<<std::endl<<"_count: "<<_count<<std::endl;

};
//########################################################################################################################################################
///Set the mode to amplitude fitting, instead of anchor wave
void anchor_t::set_is_ampl(
							bool 						is_ampl){

	_is_ampl = is_ampl;
};
//########################################################################################################################################################
///Complex conjugates _data, changes _coma accordingly
void anchor_t::conjugate(){

	for (size_t i=2;i<2*(_waveset.nPoints());i+=2){
		for (size_t tbin=0;tbin<_waveset.nTbin();tbin++){
			for(size_t bin=0;bin<_waveset.nBins();bin++){
				_data[tbin][bin][i]*=-1;
				for (size_t j=0;j<2*_waveset.nPoints()-1;j++){
					_coma[tbin][bin][i][j]*=-1;
					_coma[tbin][bin][j][i]*=-1;
				};
			};
		};
	};
};
//########################################################################################################################################################
///Initializes _data and _coma, when setting the t' binning
void anchor_t::setTbinning(std::vector<std::vector<double> > binning){
	_waveset.setTbinning(binning);
	if (_data.size() != _waveset.nTbin()){
		if (_data.size() != 0){
			std::cout<<"anchor_t::setTbinning(...): Warning: _data.size() != _waveset.nTbin(), but nonzero. Previous set _data[][][] will be lost."<<std::endl;
		};
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_coma.size() != _waveset.nTbin()){
		if (_coma.size() != 0){
			std::cout<<"anchor_t::setTbinning(...): Warning: _coma.size() != _waveset.nTbin(), bun nonzero. Previous set _coma[][][][] will be lost."<<std::endl;
		};
		_coma = std::vector<std::vector<std::vector<std::vector<double> > > >(_waveset.nTbin());
	};
};
//########################################################################################################################################################
///Updates internal definitions
void anchor_t::update_definitions(){

	_nTot = getNtotAnc();
	_nPar = _waveset.getNpar();
	_nCpl = getNanc();
	_nBra = _waveset.nBranch();
	_nIso = _waveset.getNiso();
	std::vector<std::string> names;
	std::stringstream str_count;
	for (size_t i=0;i<_nCpl;i++){
		str_count<<i;
		names.push_back("reC"+str_count.str());
		names.push_back("imC"+str_count.str());
		str_count.str("");
	};
	for (size_t i=0;i<_nPar;i++){
		names.push_back(_waveset.getParameterName(i));
	};
	for (size_t i=0;i<_nBra;i++){
		str_count<<i;
		names.push_back("reB"+str_count.str());
		names.push_back("imB"+str_count.str());
		str_count.str("");
	};
	for (size_t i=0;i<_nIso;i++){
		names.push_back(_waveset.getIsoParName(i));
	};
	_parNames = names;
};
//########################################################################################################################################################
///Updates, which data-point is actually active
void anchor_t::update_is_active(){

	_is_active = std::vector<std::vector<std::vector<bool> > >(_waveset.nTbin(),std::vector<std::vector<bool> >(_waveset.nBins(),std::vector<bool>(_waveset.nPoints(),true)));
	for (size_t tbin =0;tbin<_waveset.nTbin();tbin++){
		if (_data.size() < tbin){
			break;
		};
		for(size_t bin=0;bin<_waveset.nBins();bin++){
			if (_data[tbin].size() < bin-1){
				break;
			};
			if (_data[tbin][bin].size() < 2*_waveset.nPoints()-1){
				break;
			};
			if(_data[tbin][bin][0] == 0.){
				_is_active[tbin][bin][0] = false;
			};
			for(size_t point=1;point<_waveset.nPoints();point++){
				if (_data[tbin][bin][2*point-1] == 0. and _data[tbin][bin][2*point] ==0.){
					_is_active[tbin][bin][point] = false;
				};
			};
		};
	};
};
//########################################################################################################################################################
///Writes plots for one t' bin
void anchor_t::write_plots(
							std::string						filename,
							int 							tbin,
							const std::vector<std::complex<double> >		&cpl,
							const std::vector<double>				&par,
							const std::vector<std::complex<double> >		&bra,
							const std::vector<double> 				&iso)					const{

	std::vector<std::complex<double> > cpl_all = getAllCouplings(tbin,cpl,par,bra,iso);
	std::vector<std::vector<std::complex<double> > > iso_eval;
	if(_waveset.has_isobars()){
		iso_eval = _waveset.iso_funcs(&iso[0]);
	};
	std::ofstream write_out;
	write_out.open(filename.c_str());
	std::cout<<"anchor_t::write_plots(...): Chi2 for the used paramters is: "<<EvalTbin(tbin,&cpl_all[0],&par[0],&iso[0])<<std::endl;//[0]//
	for (size_t bin=0;bin<_waveset.nBins();bin++){
		double mass = _waveset.get_m(bin);
		std::vector<double> var = _waveset.getVar(mass,tbin);
		std::vector<std::complex<double> > amplitudes = _waveset.amps(&var[0],&cpl_all[0],&par[0],iso_eval);
		std::complex<double> ancAmp = amplitudes[0];
		write_out<<mass<<" "<<norm(ancAmp)<<" "<<_data[tbin][bin][0]<<" "<<1/sqrt(_coma[tbin][bin][0][0]);
		for (size_t i=1; i<_waveset.nPoints(); i++){
			std::complex<double> inter = ancAmp*conj(amplitudes[i]);
			write_out<<" "<<inter.real()<<" "<<_data[tbin][bin][2*i-1]<<" "<<1/sqrt(_coma[tbin][bin][2*i-1][2*i-1]);
			write_out<<" "<<inter.imag()<<" "<<_data[tbin][bin][2*i  ]<<" "<<1/sqrt(_coma[tbin][bin][2*i  ][2*i  ]);
		};
		write_out<<std::endl;
	};
	write_out.close();
};
//########################################################################################################################################################
