#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>
#include"chi2container.h"
#include<string>
namespace bp = boost::python;

template<class T>
bp::list std_vector_to_py_list(const std::vector<T>& v){
	bp::object get_iter = bp::iterator<std::vector<T> >();
	bp::object iter = get_iter(v);
	bp::list l(iter);
	return l;
};


template< typename T >
std::vector<T> to_std_vector(const bp::object& iterable){
    return std::vector<T>(bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>( ));
};

template< typename T >
std::vector<std::vector<T> > to_std_vector_vector(const bp::object& iterable){
	std::vector<bp::object> intermediate = to_std_vector<bp::object>(iterable);
	std::vector<T> line;
	std::vector<std::vector<T> > final;
	for (unsigned int i=0;i<intermediate.size();i++){
		line=to_std_vector<T>(intermediate[i]);
		final.push_back(line);
	};
	return final;	
};


struct chi2py:public chi2container{

	chi2py(std::string card):
		chi2container(card){};

	size_t nTot()  {return _minimizeChi2->method()->nTot();};
	size_t nTbin() {return _minimizeChi2->method()->Waveset()->nTbin();};
	size_t nCpl()  {return _minimizeChi2->method()->nCpl();};
	size_t nBrCpl(){return _minimizeChi2->method()->Waveset()->nBrCpl();};
	size_t nPar()  {return _minimizeChi2->method()->Waveset()->getNpar();};
	size_t nBra()  {return _minimizeChi2->method()->Waveset()->nBranch();};
	size_t nFtw()  {return _minimizeChi2->method()->Waveset()->nFtw();};
	size_t nWaves(){return _minimizeChi2->method()->Waveset()->nWaves();};

	void relPar(std::string inin){
		_minimizeChi2->relPar(inin);
	};
	void fixPar(std::string inin){
		_minimizeChi2->fixPar(inin);
	};
	void printParameters(){
		_minimizeChi2->method()->Waveset()->printParameters();
	};
	std::string method(){
		return _minimizeChi2->method()->className();
	};

	double Eval(bp::object par){
		std::vector<double> xxx = to_std_vector<double>(par);
		return (*(_minimizeChi2->method()))(xxx);
	};
	double EvalSelf(){
		return (*(_minimizeChi2->method()))(&(_minimizeChi2->method())->parameters()[0]);
	};
	double fit(){
		_minimizeChi2->fit();
	};

#ifdef ADOL_ON
	bp::list Diff(bp::object par){
		std::vector<double> xxx = to_std_vector<double>(par);
		std::vector<double> dif = _minimizeChi2->method()->Diff(xxx);
		return std_vector_to_py_list(dif);
	};
	bp::list DiffSelf(){
		std::vector<double> dif = _minimizeChi2->method()->Diff((_minimizeChi2->method())->parameters());
		return std_vector_to_py_list(dif);
	};
#endif//ADOL_ON
	bp::list DiffAuto(bp::object par){
		std::vector<double> xxx = to_std_vector<double>(par);
		std::vector<double> dif = _minimizeChi2->method()->DiffAuto(xxx);
		return std_vector_to_py_list(dif);
	};

	bp::list diff_amps(size_t nAmp, size_t tbin, double mass, bp::object par){

		std::vector<double> parameters = to_std_vector<double>(par);
		std::vector<std::vector<std::complex<double> > > rett = _minimizeChi2->method()->Waveset()->diff_amps_full(tbin,mass, &parameters[0], false);
		return std_vector_to_py_list(rett[nAmp]);
	};


	double getParameter(std::string name){
		int i = _minimizeChi2->method()->getParNumber(name);
		if (i<0){
			std::cout<<"Error: Parameter '"<<name<<"' not definded"<<std::endl;
	 		return std::numeric_limits<double>::quiet_NaN();
		};
		return _minimizeChi2->method()->getParameter(i);
	};
	void setParameter(std::string name, double val){
		_minimizeChi2->method()->setParameter(name,val);
	};
	void write_plots(std::string file_name, size_t tbin){
		_minimizeChi2->method()->write_plots(file_name,tbin);
	};

	bp::list parameterNames(){
		std::vector<std::string> names = *_minimizeChi2->method()->parNames();
		return std_vector_to_py_list(names);
	};
	bp::list parameters(){
		std::vector<double> xxx = _minimizeChi2->method()->parameters();
		return std_vector_to_py_list(xxx);
	};
	bp::list fullParameters(){
		std::vector<double> xxx = _minimizeChi2->method()->fullParameters();
		return std_vector_to_py_list(xxx);
	};
	bp::list Amplitudes(double mass, int tbin, bp::object param_in){
		std::vector<double> par = to_std_vector<double>(param_in);
		std::vector<std::complex<double> > amps = _minimizeChi2->method()->amplitudes(mass, tbin, par, true);
		return std_vector_to_py_list(amps);
	};
	bp::list AmplitudesSelf(double mass, int tbin){
		std::vector<std::complex<double> > ampl =_minimizeChi2->method()->amplitudes(mass, tbin, _minimizeChi2->method()->parameters(), true);
		return std_vector_to_py_list(ampl);
	};
	bp::list waveNames(){
		std::vector<std::string> names = *_minimizeChi2->method()->Waveset()->waveNames();
		return std_vector_to_py_list(names);
	};
	bp::list borders_waves(){
		std::vector<size_t> brd = *_minimizeChi2->method()->Waveset()->borders_waves();
		return std_vector_to_py_list(brd);
	};
	bp::list getFuncParameters(int ftw){
		std::vector<size_t> params = _minimizeChi2->method()->Waveset()->getFuncParameters(ftw);
		return std_vector_to_py_list(params);
	};
	bp::list upperLims(){
		std::vector<double> lims = *_minimizeChi2->method()->Waveset()->upperLims();
		return std_vector_to_py_list(lims);
	};
	bp::list lowerLims(){
		std::vector<double> lims = *_minimizeChi2->method()->Waveset()->lowerLims();
		return std_vector_to_py_list(lims);
	};


	void setParLimits(std::string name, double upper, double lower){
		_minimizeChi2->method()->setParLimits(name,upper,lower);
	};
	void setParameterFile(std::string fileName){
		_minimizeChi2->method()->setParameterFile(fileName);
	};
	void setNout(size_t n){
		_minimizeChi2->method()->setNoutFile(n);
	};
	void writeParameters(std::string fileName){
		_minimizeChi2->method()->writeParameters(fileName);
	};
	void readParameters(std::string fileName){
		_minimizeChi2->method()->readParameters(fileName);
	};
	void initCouplings(size_t nseeds){
		_minimizeChi2->initCouplings(nseeds,-1);
	};
	void initSingleTbin(size_t nseeds, size_t tbin){
		_minimizeChi2->initCouplings(nseeds,tbin);
	};
	std::string YAML_file(){
		return _minimizeChi2->method()->Waveset()->YAML_file();
	};
	std::string get_component_name(size_t i){
		return _minimizeChi2->method()->Waveset()->get_component_name(i);
	};
	void setRandomCpl(){
		 _minimizeChi2->setRandomCpl();
	};
	void setRandomBra(){
		 _minimizeChi2->setRandomBra();
	};
	void setMaxCalls(size_t nCalls){
		 _minimizeChi2->setMaxCalls(nCalls);
	};
};

struct chi2py_wrap: chi2py, bp::wrapper<chi2py>{
	
	chi2py_wrap(const chi2py &inin) : chi2py(inin), bp::wrapper<chi2py>(){};
	chi2py_wrap(std::string inin)   : chi2py(inin), bp::wrapper<chi2py>(){};

};


int wolen(){
	return 0;
};

BOOST_PYTHON_MODULE(libchi2py){
	// Expose the std::vectors<...> to python
	
	bp::class_<std::vector<double> >("vector_double")
		.def(bp::vector_indexing_suite<std::vector<double> >())
	;
	bp::class_<std::vector<bool> >("vector_bool")
		.def(bp::vector_indexing_suite<std::vector<bool> >())
	;
	bp::class_<std::vector<std::string> >("vector_string")
		.def(bp::vector_indexing_suite<std::vector<std::string> >())
	;
	bp::class_<std::vector<int> >("vector_int")
		.def(bp::vector_indexing_suite<std::vector<int> >())
	;
	bp::class_<std::vector<size_t> >("vector_size_t")
		.def(bp::vector_indexing_suite<std::vector<size_t> >())
	;

	bp::class_<std::vector<std::complex<double> > >("vector_complex")
		.def(bp::vector_indexing_suite<std::vector<std::complex<double> > >())
	;

	bp::class_<std::vector<std::vector<double> > >("vector_array_double")
		.def(bp::vector_indexing_suite<std::vector<std::vector<double> > >())
	;


	bp::class_<std::vector<std::vector<std::vector<double> > > >("vector_array_array_double") // To be able to have a vector of covariance matrices.
		.def(bp::vector_indexing_suite<std::vector<std::vector<std::vector<double> > > >())
	;

	bp::class_<chi2py_wrap> chi2 = bp::class_<chi2py_wrap>("chi2",bp::init<std::string>());


	chi2.def("Eval",                                &chi2py::Eval                                   );
	chi2.def("EvalSelf",                            &chi2py::EvalSelf                               );

	chi2.def("Amplitudes",                          &chi2py::Amplitudes                             );
	chi2.def("AmplitudesSelf",                      &chi2py::AmplitudesSelf                         );

	chi2.def("Diff",                                &chi2py::Diff                                   );
	chi2.def("DiffSelf",                            &chi2py::DiffSelf                               );
#ifdef ADOL_ON
	chi2.def("DiffAuto",                            &chi2py::DiffAuto                               );
#endif//ADOL_ON
	chi2.def("nTot",                                &chi2py::nTot                                   );
	chi2.def("nCpl",                                &chi2py::nCpl                                   );
	chi2.def("nBrCpl",                              &chi2py::nBrCpl                                 );
	chi2.def("nPar",                                &chi2py::nPar                                   );
	chi2.def("nBra",                                &chi2py::nBra                                   );
	chi2.def("nFtw",                                &chi2py::nFtw                                   );
	chi2.def("nTbin",                               &chi2py::nTbin                                  );
	chi2.def("nWaves",                              &chi2py::nWaves                                 );

	chi2.def("diff_amps",                           &chi2py::diff_amps                              );

	chi2.def("method",                              &chi2py::method                                 );
	chi2.def("className",                           &chi2py::className                              );
	chi2.def("printStatus",                         &chi2py::printStatus                            );
	chi2.def("printParameters",                     &chi2py::printParameters                        );
	chi2.def("setParameter",                        &chi2py::setParameter                           );
	chi2.def("getParameter",                        &chi2py::getParameter                           );
	
	chi2.def("parameterNames",                      &chi2py::parameterNames                         );
	chi2.def("setRandomCpl",                        &chi2py::setRandomCpl                           );
	chi2.def("setRandomBra",                        &chi2py::setRandomBra                           );

	chi2.def("fit",                                 &chi2py::fit                                    );
	chi2.def("initCouplings",                       &chi2py::initCouplings                          );
	chi2.def("initSingleTbin",                      &chi2py::initSingleTbin                         );
	chi2.def("relPar",                              &chi2py::relPar                                 );
	chi2.def("fixPar",                              &chi2py::fixPar                                 );
	chi2.def("parameters",                          &chi2py::parameters                             );
	chi2.def("fullParameters",                      &chi2py::fullParameters                         );

	chi2.def("write_plots",                         &chi2py::write_plots                            );
	chi2.def("setParLimits",                        &chi2py::setParLimits                           );
	chi2.def("setParameterFile",                    &chi2py::setParameterFile                       );
	chi2.def("setNout",                             &chi2py::setNout                                );
	chi2.def("writeParameters",                     &chi2py::writeParameters                        );
	chi2.def("readParameters",                      &chi2py::readParameters                         );
	chi2.def("setMaxCalls",                         &chi2py::setMaxCalls                            );
	chi2.def("YAML_file",                           &chi2py::YAML_file                              );
	chi2.def("waveNames",                           &chi2py::waveNames                              );
	chi2.def("borders_waves",                       &chi2py::borders_waves                          );
	chi2.def("lowerLims",                           &chi2py::lowerLims                              );
	chi2.def("upperLims",                           &chi2py::upperLims                              );
	chi2.def("get_component_name",                  &chi2py::get_component_name                     );
	chi2.def("getFuncParameters",                   &chi2py::getFuncParameters                      );
	chi2.def("setMinimizerSpecifications",          &chi2py::setMinimizerSpecifications             );
};























