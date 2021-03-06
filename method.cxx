#include"method.h"
#include"matrix_utilities.h"
//#######################################################################################################################################################
///Set a single parameter by number
void method::setParameter(
							size_t 						i, 	// # of parameter
							double 						par){

	if(_parameters.size() < 2*(_nCpl+_nBra)){
		_parameters = std::vector<double>(2*(_nCpl+_nBra),0.);
	};
	if (i<2*_nCpl){
		_parameters[i] = par;
	}else if(i<2*_nCpl+_nPar){
		_waveset.setPar(i-2*_nCpl,par);
	}else if(i<2*_nCpl+_nPar+2*_nBra){
		_parameters[i-_nPar] = par;
	}else if(i<2*_nCpl+_nPar+2*_nBra+_nIso){
		_waveset.setIsoPar(i-2*_nCpl-_nPar-2*_nBra,par);
	}else{
		std::cerr<<"method::setParameter(...): Error: Can't set parameter #"<<i<<std::endl;
	};
};
//#######################################################################################################################################################
double method::getParameter(				size_t 						i)						const{

	if(i<2*_nCpl){
		return _parameters[i];
	}else if(i<2*_nCpl+_nPar){
		return _waveset.getPar(i-2*_nCpl);
	}else if(i<2*_nCpl+_nPar+2*_nBra){
		return _parameters[i-_nPar];
	}else if(i<2*_nCpl+_nPar+2*_nBra+_nIso){
		return _waveset.getIsoPar(i-2*_nCpl-_nPar-2*_nBra);
	}else{
		std::cerr<<"method::getParameter(...): Error: Can't get parameter #"<<i<<std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	};
};
//#######################################################################################################################################################
///Set internal parameters
void method::setParameters(
							std::vector<double> 				pars){

	if(pars.size()>=_nTot){
		for (size_t i=0;i<_nTot;i++){
			setParameter(i,pars[i]);
		};
	}else{
		std::cerr<<"method::setParameters(...): Error: Input pars.size() too small"<<std::endl;
	};
};
//#######################################################################################################################################################
///Get parameters
const std::vector<double> method::parameters()														const{

	std::vector<double> par(_nTot);
	for (size_t i=0;i<_nTot;i++){
		par[i] = getParameter(i);
	};
	return par;
};
//#######################################################################################################################################################
///Set a single parameter by name
bool method::setParameter(
							std::string 					name,
							double 						par){

	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "method::setParameter(...): Error: Parameter '"<<name<<"' not found"<<std::endl;
		return false;
	}else if ((int)_nTot <= number){
		std::cerr << "method::setParameter(...): Error: Parameter number too high"<<std::endl;
		return false;
	}else{
		setParameter(number,par);
		return true;
	};
};
//########################################################################################################################################################
///Gets the number for a given name, if the name does not exist, return -1
int method::getParNumber(
							std::string 					name)						const{

	for (size_t i=0;i<_parNames.size();i++){
		if (_parNames[i] == name){
			return i;
		};
	};
	return -1;
};
//########################################################################################################################################################
///Set Parameter limits. If upper < lower, no limits will be applied
void method::setParLimits(
							int 						par,
							double 						upper,
							double 						lower){

	if(_lower_parameter_limits.size() != _nTot){
		init_lower_limits();
	};
	if(_upper_parameter_limits.size() != _nTot){
		init_upper_limits();
	};
	if (upper < lower){
		std::cout<<"method::setParLimits(...): Warning: Trying to set limits for parameter number #"<<par<<": Upper limit smaller than lower limit( "<<upper<<" < "<<lower<<" )"<<std::endl;
	};
	_upper_parameter_limits[par]=upper;
	_lower_parameter_limits[par]=lower;
};
//########################################################################################################################################################
void method::setParLimits(
							std::string					name,
							double 						upper,
							double 						lower){

	int number = getParNumber(name);
	if (-1==number){
		std::cerr << "method::setParLimits: Error: Parameter '"<<name<<"' not found"<<std::endl;
	}else if ((int)_nTot <= number){
		std::cerr << "method::setParLimits: Error: Parameter number too high"<<std::endl;
	}else{
		setParLimits(number,upper,lower);
	};
};
//########################################################################################################################################################
///Inits all lower parateter limits to 1.
void method::init_lower_limits(int n){

	if (n==-1){
		n=_nTot;
	};
	_lower_parameter_limits = std::vector<double>(n,1.);
};
//########################################################################################################################################################
///Inits all upper parateter limits to 0.
void method::init_upper_limits(int n){

	if(n==-1){
		n=_nTot;
	};
	_upper_parameter_limits = std::vector<double>(n,0.);
};
//########################################################################################################################################################
/// Writes parametes to a YAML file (#ifdef USE_YAML
void method::writeParameters(				const double* 					param,
							std::string					fileName)					const{

#ifdef USE_YAML
	if (fileName.size() ==0){
		std::cout<<"method::writeParameters(...): Warning: No fileName given"<<std::endl;
	}else{
		if (not _parameterFile.size() == 0){
			YAML::Node YAML_parameters;
			for (size_t i=0;i<_nTot;++i){
				YAML_parameters[_parNames[i]] = parameters()[i];
			};
			std::ofstream fout(fileName.c_str());
			fout<<YAML_parameters;
			fout.close();	
			std::cout<<"method::writeParameters(...): Parameters written to '"<<fileName<<"'"<<std::endl;
		};
	};
#else//USE_YAML
	std::cout<<"method::writeParameters(): Error: YAML not used, do nothig"<<std::endl;
#endif//USE_YAML
};
//########################################################################################################################################################
/// Wirtes internal paramters to a YAML file
void method::writeParameters(				std::string					fileName)					const{

	writeParameters(&parameters()[0], fileName);
};
//########################################################################################################################################################
/// Reads parameters from a YAML file
void method::readParameters(				std::string					fileName){
#ifdef USE_YAML
	YAML::Node YAML_parameters = YAML::LoadFile(fileName);
	for( YAML::const_iterator it=YAML_parameters.begin(); it != YAML_parameters.end(); ++it){
		setParameter(it->first.as<std::string>(), it->second.as<double>());
	};
#else//USE_YAML
	std::cout<<"method::readParameters(): Error: YAML not used, do nothig"<<std::endl;
#endif//USE_YAML
};
//########################################################################################################################################################
///Gets the couplings-without-branchings from couplings-with-branchings and branchings for one t' bin
std::vector<std::complex<double> > method::getUnbranchedCouplings(
							const std::vector<std::complex<double> >	&cpl,
							const std::vector<std::complex<double> >	&bra)						const{
	
	std::vector<std::complex<double> > cpl_un(_waveset.nFtw());
	for (size_t i=0;i<_waveset.nFtw();i++){
		if(-1==(*_waveset.n_branch())[i]){
			cpl_un[i] = cpl[(*_waveset.n_cpls())[i]];
		}else{
			cpl_un[i] = cpl[(*_waveset.n_cpls())[i]] * bra[(*_waveset.n_branch())[i]];
		};
	};
	return cpl_un;
};
//########################################################################################################################################################
///Updates the internal ranges for evaluation
void method::update_min_max_bin(){

// Override the mthod in waveset, because the limits are given by the anchor wave in this method
	_waveset.setMinBin(0);
//	std::cout<<"Call update_min_max_bin()"<<std::endl;
	_waveset.setMaxBin((*_waveset.binning()).size());
	if ((*_waveset.lowerLims()).size() == 0 or (*_waveset.upperLims()).size() == 0){
		std::cout<< "method::update_min_max_bin(): Warning: Wave limits not set, omitting internal setting of _waveset.minBin() and _waveset.maxBin()"<<std::endl;
		std::cout<< "method::update_min_max_bin(): Warning: Setting _waveset.minBin() = 0; _waveset.maxBin() = _binning.size() = "<<(*_waveset.binning()).size()<<std::endl;
		_waveset.setMinBin(0);
		_waveset.setMaxBin((*_waveset.binning()).size());
		return;
	};
	double ancMin = (*_waveset.lowerLims())[0];
	double ancMax = (*_waveset.upperLims())[0];
//	std::cout<<ancMin<<"-"<<ancMax<<std::endl;
	if ((*_waveset.binning()).size() == 0){
		std::cout<< "method::update_min_max_bin(): Warning: No binning set, cannot determine _waveset.minBin(), _waveset.maxBin()" <<std::endl;
	};
	for (size_t i=0;i<(*_waveset.binning()).size()-1;i++){
		double up = (*_waveset.binning())[i+1];
		double low= (*_waveset.binning())[i];
		if (ancMin >= low and ancMin <up){
			_waveset.setMinBin(i);
		};
		if (ancMax > low and ancMax <=up){
			_waveset.setMaxBin(i+1);
		};
	};
	update_n_cpls();
//	std::cout<<"_waveset.minBin(): "<<_waveset.minBin()<<std::endl;
//	std::cout<<"_waveset.maxBin(): "<<_waveset.maxBin()<<std::endl;
};
//########################################################################################################################################################
///Updates internal definitions
void method::update_definitions(){
	_nTot = getNtot();
	_nPar = _waveset.getNpar();
	_nCpl = getNcpl();
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
	size_t nFtw = _waveset.nFtw();
	size_t nTbin = _waveset.nTbin();
	size_t wave = 0;
	size_t cpl_per_t = _nCpl/nTbin;
	for (size_t ftw = 0;ftw<nFtw;++ftw){
		if(_waveset.borders_waves()->at(wave) == ftw){
			++wave;
		};
		size_t n_cpl = _waveset.n_cpls()->at(ftw);
		size_t n_func = _waveset.funcs_to_waves()->at(ftw);
		int n_bra = _waveset.n_branch()->at(ftw);
		std::string addname = _waveset.amp_funcs()->at(n_func)->name()+"/"+_waveset.waveNames()->at(wave);

		for (size_t tbin=0;tbin<nTbin;++tbin){
			str_count<<tbin;
//			std::cout<<cpl_per_t<<" "<<n_cpl<<" "<<_nCpl<<" uiubai "<<2*(tbin*cpl_per_t+n_cpl)<<std::endl;
			if (n_cpl < cpl_per_t){ // To prevent auto cpl methods from renaming parameters they don't use
				names[2*(tbin*cpl_per_t+n_cpl)  ] = "reC/"+str_count.str()+"/"+addname;
				names[2*(tbin*cpl_per_t+n_cpl)+1] = "imC/"+str_count.str()+"/"+addname;
			};
			str_count.str("");
		};
		if (n_bra >-1){
			names[2*_nCpl+_nPar+2*n_bra  ] = "reB/"+addname;
			names[2*_nCpl+_nPar+2*n_bra+1] = "imB/"+addname;
		};



	};
//	for (size_t i=0;i<_nPar;i++){
//		names[2*_nCpl+i]=_waveset.getParameterName(i);
//	};
	_parNames = names;
};
//########################################################################################################################################################
#ifdef USE_YAML
//Loads _data and _coma from a YAML file
bool method::loadDataComa(
							YAML::Node					&waveset){


	if (waveset["data_files"].size() != _waveset.nTbin() ){
		std::cerr<<"method::loadDataComa(...): Error: Number of data-files does not match number of t' bins"<<std::endl;
	};
	if (waveset["coma_files"].size() != _waveset.nTbin()){
		std::cerr<<"method::loadDataComa(...): Error: Number of coma-files does not match number of t' bins"<<std::endl;
	};
	update_min_max_bin(); //Has to be called to set internal definitions right, since not all belong to the same class now
	update_definitions(); // Has to be called once with all functions set, to get internal handling of parameter numbers right
	if (waveset["data_files"].size() == _waveset.nTbin() and waveset["coma_files"].size() == _waveset.nTbin()){
		for(size_t i=0;i<_waveset.nTbin();i++){
			loadData(i,waveset["data_files"][i].as<std::string>().c_str());
			loadComa(i,waveset["coma_files"][i].as<std::string>().c_str());
		};
	}else{
		std::cerr<<"method::loadDataComa(...): Error: Wrong number of files, no _data or _coma loaded"<<std::endl;
		return false;
	};
	if (waveset["conjugate_fit_result"]){
		if (waveset["conjugate_fit_result"].as<bool>()){
			conjugate();
		};
	};
	return true;
};
//########################################################################################################################################################
///Loads paramters values and limits from a YAML file
bool method::loadParameterValues(
							YAML::Node					&waveset,
							YAML::Node					&param){

	std::map<std::string,int> fMap;
	int fCount = 1; // Starts as one, since map[nonExistingKey] = 0 -> use this to determine, whether a function exists.
	int iCount = 1;
	int pCount = 0;
	int ipCount = 0;
	int nWaves = waveset["waves"].size();
	bool ookk = true;
	for (int i=0;i<nWaves;i++){
		int nFunc = waveset["waves"][i]["parametrizations"].size();
		for (int j=0;j<nFunc;j++){
			std::string fName;
			std::string iName;
			bool isisobar = false;
			if (waveset["waves"][i]["parametrizations"][j].size()==0){
				fName = waveset["waves"][i]["parametrizations"][j].as<std::string>();
//				std::cout<<fName<<" alone"<<std::endl;
			}else if(waveset["waves"][i]["parametrizations"][j].size()==2) {
				isisobar = true;
				fName = waveset["waves"][i]["parametrizations"][j][0].as<std::string>();
				iName = waveset["waves"][i]["parametrizations"][j][1].as<std::string>();
//				std::cout<<iName<<" to "<<fName<<std::endl;
			};
			if(param[fName]){
				if (fMap[fName]==0){
					fMap[fName]=fCount;
					fCount++;
					int nPar = param[fName]["parameters"].size();
					for(int par=0;par<nPar;par++){
						double value = param[fName]["parameters"][par]["value"].as<double>();
						if(param[fName]["parameters"][par]["upper_limit"] and param[fName]["parameters"][par]["lower_limit"]){
							double upper = param[fName]["parameters"][par]["upper_limit"].as<double>();
							double lower = param[fName]["parameters"][par]["lower_limit"].as<double>();
							setParLimits(2*_nCpl+pCount,upper,lower);
						};
						setParameter(2*_nCpl+pCount,value);
						pCount++;
					};
				};
			}else{
				ookk = false;
				std::cerr << "method::loadParameterValues(...): Error: '"<<fName<<"' not defined in parametrization file, no values set"<<std::endl;
			};
			if (isisobar){
				if(param[iName]){
					if(fMap[iName]==0){
						fMap[iName]=iCount;
						iCount++;
						int nPar =  param[iName]["parameters"].size();
						for(int par=0;par<nPar;par++){
							double value = param[iName]["parameters"][par]["value"].as<double>();
							if(param[iName]["parameters"][par]["upper_limit"] and param[iName]["parameters"][par]["lower_limit"]){
								double upper = param[iName]["parameters"][par]["upper_limit"].as<double>();
								double lower = param[iName]["parameters"][par]["lower_limit"].as<double>();
								setParLimits(2*_nCpl+_nPar+2*_nBra+ipCount,upper,lower);
							};
							setParameter(2*_nCpl+_nPar+2*_nBra+ipCount,value);
							ipCount++;
						};
					};
				}else{
					ookk = false;
				};
			};
		};
	};
	// Load explicit definitions
	if (waveset["parameters"]){
		int nDef = waveset["parameters"].size();
		for (int i=0;i<nDef;i++){
			std::string name = waveset["parameters"][i]["name"].as<std::string>();
			int par = getParNumber(name);
			if (par==-1){
				std::cout<<"method::loadParameterValues(...): Warning: Unknown parameter in YAML file: "<<name<<std::endl;
				ookk = false;
			}else{
				if (waveset["parameters"][i]["value"]){
					setParameter(name,waveset["parameters"][i]["value"].as<double>());
				};
				if (waveset["parameters"][i]["upper_limit"] and waveset["parameters"][i]["lower_limit"]){
					setParLimits(par,waveset["parameters"][i]["upper_limit"].as<double>(),waveset["parameters"][i]["lower_limit"].as<double>());
				};
			};
		};
	};
	return ookk;
};
#endif//USE_YAML
//########################################################################################################################################################
///Writes plots with the internal _parameters
void method::write_plots(
							std::string						filename,
							int							tbin)					const{

	write_plots(filename,tbin,parameters());
};
//########################################################################################################################################################
///Write plots with the usual parameter ordering
void method::write_plots(
							std::string						filename,
							int							tbin,
							const std::vector<double>				&param)					const{

	std::vector<std::complex<double> > cpl(_nCpl);
	std::vector<double> par(_nPar);	
	std::vector<std::complex<double> > bra(_nBra);
	std::vector<double> iso(_nIso);
	int count =0;
	for (size_t i=0;i<_nCpl;i++){
		cpl[i] = std::complex<double>(param[count],param[count+1]);
		count+=2;
	};
	for (size_t i=0;i<_nPar;i++){
		par[i] = param[count];
		count++;
	};
	for (size_t i=0;i<_nBra;i++){
		bra[i] = std::complex<double>(param[count],param[count+1]);
		count+=2;
	};
	for (size_t i=0;i<_nIso;i++){
		iso[i]=param[count];
		count++;
	};
	write_plots(filename,tbin,cpl,par,bra,iso);
};
//########################################################################################################################################################
///Get the amplitudes for internal parameters
std::vector<std::complex<double> > method::amplitudes( 	double 							mass, 
							int							tbin,
							bool							ignore_limits)				const{
	return amplitudes(mass, tbin, parameters());
};
//########################################################################################################################################################

std::vector<std::complex<double> > method::amplitudes(	double							mass,
							int							tbin,
							const std::vector<double>				&param,
							bool							ignore_limits)				const{

	std::vector<double> par(_nPar);	
	std::vector<std::complex<double> > bra(_nBra);
	std::vector<double> iso(_nIso);
	size_t count =0;
	size_t vec_size = param.size();
	for (size_t i=1;i<=_nIso;++i){ // Di the inverse direction
		count+=1;
		iso[_nIso-i] = param[vec_size-count];

	};
	for (size_t i=1;i<=_nBra;++i){
		count+=1;
		bra[_nBra-i] = std::complex<double> (param[vec_size-count-1],param[vec_size-count]);
		count+=1;
	};
	for (size_t i=1;i<=_nPar;++i){
		count+=1;
		par[_nPar-i] =  param[vec_size-count];
		
	};
	size_t diff = vec_size - count;
	std::vector<std::complex<double> > cpl(diff/2);
	for (size_t i=0;i<diff/2;++i){
		cpl[i] = std::complex<double>(param[2*i],param[2*i+1]);
	};
	std::vector<std::complex<double> > cpl_all = getAllCouplings(tbin,cpl,par,bra,iso);
	std::vector<std::vector<std::complex<double> > > iso_eval;
	if(_waveset.has_isobars()){
		iso_eval = _waveset.iso_funcs(&iso[0]);
	};
	std::vector<double> var = _waveset.getVar(mass,tbin);
	std::vector<std::complex<double> > amplitudes = _waveset.amps(&var[0],&cpl_all[0],&par[0],iso_eval,ignore_limits);
	return amplitudes;
};
//########################################################################################################################################################
