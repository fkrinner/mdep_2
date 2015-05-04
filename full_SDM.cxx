#include"full_SDM.h"
#include"matrix_utilities.h"


//########################################################################################################################################################
///Calculates some initial values for the branchigs
std::vector<std::complex<double> > full_SDM::get_branchings(
							const std::vector<std::complex<double> > 	&cpl,
							const std::vector<double> 			&par = std::vector<double>(),
							const std::vector<double> 			&iso_par = std::vector<double>())		const{

	


	throw; // Does not work at the moment

	return std::vector<std::complex<double> >(_nBra,std::complex<double>(1.,0.));;
};
//########################################################################################################################################################
///Gets all couplings for one t' bin (.size() = _waveset.nFtw()) for different input formats
std::vector<std::complex<double> > full_SDM::getAllCouplings(
							int						tbin,
							const std::vector<std::complex<double> >	&cpl,
							const std::vector<double>			&par,
							const std::vector<std::complex<double> >	&bra,
							const std::vector<double>			&iso)						const{


	size_t nTbin = _waveset.nTbin();
	size_t nFtw = _waveset.nFtw();
	size_t cpls = _nCpl/nTbin;
	std::vector<std::complex<double> > cpl_all(nFtw);
	if (cpl.size() == nFtw*nTbin){
//		std::cout<<"Got all couplings for all 't bins"<<std::endl;
		for (size_t i=0;i<nFtw;++i){
			cpl_all[i] = std::complex<double>(cpl[nFtw*tbin+i]);
		};	
	} else if (cpl.size() == _nCpl){
//		std::cout<< "Got branched couplings for all t' bins"<<std::endl;
		std::vector<std::complex<double> > tCpls(cpls);
		for (size_t i=0;i<cpls;++i){
			tCpls[i]=cpl[tbin*cpls+i];
		};
		cpl_all = getUnbranchedCouplings(tCpls,bra);
	} else if (cpl.size() == nFtw){
//		std::cout<<"Got all couplings for on t' bin"<<std::endl;
		cpl_all = cpl;
	} else if (cpl.size() == cpls){
//		std::cout<<"Got branched couplings for one t' bin"<<std::endl;
		cpl_all = getUnbranchedCouplings(cpl,bra);
	}else{
		throw std::invalid_argument("Unknown format for couplings");
	};
//	std::cout<<"SDM: ";
//	print_vector(cpl_all);
	return cpl_all;
};
//########################################################################################################################################################
///Set the anchor coupligs to one, that are coupled to another function
void full_SDM::branchCouplingsToOne(){
	for (size_t i=0;i<(*_waveset.coupled()).size();i++){
		if((*_waveset.coupled())[i] >=0){
			size_t nCpl = (*_waveset.n_cpls())[i];
			if (nCpl <= _nBrCpl){
				for (size_t tbin =0;tbin<_waveset.nTbin();tbin++){
					setParameter(2*_nBrCpl*tbin+2*nCpl  ,1.);
					setParameter(2*_nBrCpl*tbin+2*nCpl+1,0.);
				};
			};
		};
	};
};
//########################################################################################################################################################
///Set data points manually
bool full_SDM::set_data(
							int								tbin,
							int 								bin,
							std::vector<double> 						data){

	// Ensure right number of bins
	if (_data.size() != _waveset.nTbin()){
		_data = std::vector<std::vector<std::vector<double> > >(_waveset.nTbin());
	};
	if (_data[tbin].size() != _waveset.nBins()){
		_data[tbin] = std::vector<std::vector<double> >(_waveset.nBins());
	};
	_data[tbin][bin] = data;
	update_is_active();

	if (data.size() == _waveset.nPoints()*_waveset.nPoints()){
		return true;
	};
	return false;
};
//########################################################################################################################################################
///Loads the data for a t' bin from a file
void full_SDM::loadData(
							int 								tbin,
							const char* 							dataFile){

	_data[tbin] = std::vector<std::vector<double> >();
	std::fstream data(dataFile,std::ios_base::in);
	double val;
	size_t nDat = _waveset.nPoints()*_waveset.nPoints();
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
		std::cout << "'full_SDM.cxx' loadData(...): Warning: _waveset.nBins()="<<_waveset.nBins()<<" != _data.size()="<<_data[tbin].size()<<std::endl;
	}else{
		std::cout << "'full_SDM.cxx' loadData(...): File delivered the right size for _data"<<std::endl;
	};
};
//########################################################################################################################################################
///Returns the number of couplings
size_t full_SDM::getNcpl()																const{
	std::vector<size_t> nnn = *_waveset.n_cpls();
	size_t n = nnn.size();
	size_t maxx = 0;
	for (size_t i=0;i<n;++i){
		if (nnn[i] > maxx){
			maxx = nnn[i];
		};
	};
	maxx+=1;
	maxx *= _waveset.nTbin();
	return maxx;
};
//########################################################################################################################################################
/// Sets br and cpl in one t' bin so, that it corresponds tothe full cpls
std::vector<std::vector<std::complex<double> > > full_SDM::full_to_br_cpl(std::vector<std::complex<double> > &cpl){

	size_t nFtw  = _waveset.nFtw();
	size_t nTbin = _waveset.nTbin();
	size_t nCpl  = _nCpl;
	size_t cpls  = nCpl/nTbin;
	size_t nBra  = _waveset.nBranch();

	std::vector<std::vector<std::complex<double> > > ret;
	ret.push_back(std::vector<std::complex<double> >(cpls,std::complex<double>(0.,0.)));
	ret.push_back(std::vector<std::complex<double> >(nBra,std::complex<double>(0.,0.)));

	std::vector<int> n_branch= *_waveset.n_branch();
	std::vector<size_t> n_cpls = *_waveset.n_cpls();

	std::complex<double> phase = conj(cpl[0])/abs(cpl[0]);

	for (size_t func=0;func<nFtw;++func){
		int nBr = n_branch[func];
		size_t nCp = n_cpls[func];
		if (nBr ==-1){
			ret[0][nCp] = cpl[func]*phase;
		}else{
			if(ret[0][nCp] == std::complex<double>(0.,0.)){
				ret[0][nCp] = cpl[func]*phase;
				ret[1][nBr] = std::complex<double>(1.,0.);
			}else{
				ret[1][nBr] = cpl[func]*phase/ret[0][nCp];
			};
		};
	};
	return ret;
};
//########################################################################################################################################################
///Gets the total number of paramters with only anchor couplings
size_t full_SDM::getNtot()																const{

	return 2*getNcpl() + _waveset.getNpar() + 2*_waveset.nBranch() + _waveset.getNiso();
};
//########################################################################################################################################################
///Sets all data to 0.
void full_SDM::nullify(){

	for (size_t tbin = 0; tbin < _waveset.nTbin();tbin++){
		for (size_t bin =0; bin<_waveset.nBins();bin++){
			int nDat = _waveset.nPoints()*_waveset.nPoints();
			for (int i=0;i<nDat;i++){
				_data[tbin][bin][i]=0.;
			};
		};
	};
};
//########################################################################################################################################################
///Updates the number of couplings
void full_SDM::update_n_cpls(){

	std::vector<size_t> nnn = *_waveset.n_cpls();
	size_t n = nnn.size();
	size_t maxx = 0;
	for (size_t i=0;i<n;++i){
		if (nnn[i] > maxx){
			maxx = nnn[i];
		};
	};
	maxx+=1;
	_nBrCpl =  maxx;
};
//########################################################################################################################################################
/// Prints the internal status
void full_SDM::printStatus()																const{
	_waveset.printStatus();
	std::cout<<std::endl<<"FULL_SDM:"<<std::endl;
	std::cout<<std::endl<<"PARAMETER NUMBERS:"<<std::endl;
	std::cout<<"_nTot: "<<_nTot<<"; _nPar: "<<_nPar<<"; _nCpl: "<<_nCpl<<std::endl;
	std::cout<<"_nBra: "<<_nBra<<"; _nIso: "<<_nIso<< "; _nBrCpl: "<<_nBrCpl<<std::endl;
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
	std::cout<<std::endl<<"_nOut: "<<_nOut<<std::endl;
	std::cout<<std::endl<<"_count: "<<_count<<std::endl;
};
//########################################################################################################################################################
///Updates, which data-point is actually active
void full_SDM::update_is_active(){ //<<need>>

	_is_active = std::vector<std::vector<bool> >(_waveset.nBins(),std::vector<bool>(_waveset.nPoints(),true));
	size_t  nPoint = _waveset.nPoints();
	std::vector<int> point_to_wave = *(_waveset.point_to_wave());
	for(size_t bin=0;bin<_waveset.nBins();++bin){
		double mass = _waveset.get_m(bin);
		for (size_t point = 0;point<nPoint;++point){
			size_t wave = point_to_wave[point];
			double upperLim = (*_waveset.upperLims())[wave];
			double lowerLim = (*_waveset.lowerLims())[wave];
			if (mass > upperLim or mass < lowerLim){
				_is_active[bin][point] = false;
			};
		};
	};
};
//########################################################################################################################################################
