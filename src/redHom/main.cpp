#include <boost/algorithm/string/join.hpp>
#include <RedHomOptions.hpp>

#include <iostream>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/for_each.hpp>

#include <redHom/complex/scomplex/SComplex.hpp>
#include <redHom/complex/scomplex/SComplexReader.hpp>
#include <redHom/complex/scomplex/SComplexDefaultTraits.hpp>
#include <redHom/complex/scomplex/SComplexBuilderFromSimplices.hpp>

#include <redHom/complex/cubical/CubSComplex.hpp>
#include <redHom/complex/cubical/CubSComplexReader.hpp>
#include <redHom/algorithm/Algorithms.hpp>
#include <redHom/algorithm/strategy/DefaultReduceStrategy_CubSComplex.hpp>


enum ComplexMajorId {SComplexCubesMajorId = 0, CubicalMajorId = 1};

typedef boost::tuple<ComplexMajorId, int, boost::any> ComplexTuple;

template<typename Complex, ComplexMajorId MAJOR_ID, int MINOR_ID>
struct ComplexDescriptor {
  typedef Complex type;
  static const ComplexMajorId majorId = MAJOR_ID;
  static const int minorId = MINOR_ID;
};



typedef boost::mpl::vector<ComplexDescriptor<CubSComplex<2>, CubicalMajorId, 2>,
			   ComplexDescriptor<CubSComplex<3>, CubicalMajorId, 3>,
			   ComplexDescriptor<CubSComplex<4>, CubicalMajorId, 4>
			   > CubComplexDescriptors;
typedef boost::mpl::joint_view<CubComplexDescriptors,
			       boost::mpl::vector<
				 ComplexDescriptor<SComplex<SComplexDefaultTraits>, SComplexCubesMajorId, 0>
				 > > ComplexDescriptors;


class GenericCubSComplexReader {

  struct ComplexConstructor {
    ComplexTuple& complex;
    const int dim;
    RepSet<ElementaryCube>& repSet;

    ComplexConstructor(int _dim, ComplexTuple& _complex, 
		  RepSet<ElementaryCube>& _repSet): 

      dim(_dim), complex(_complex), repSet(_repSet) {
      BOOST_ASSERT(boost::get<2>(complex).empty());
    }

    template<typename Complex, ComplexMajorId MajorId, int MinorId>
    void operator()(ComplexDescriptor<Complex, MajorId, MinorId>) {
      if (MajorId == CubicalMajorId && dim == MinorId) {
	BOOST_ASSERT(boost::get<2>(complex).empty());      
	complex = ComplexTuple(MajorId, MinorId, 
			       boost::any(boost::shared_ptr<Complex>(new Complex(repSet))));
      }
    }

  };

public:
  static ComplexTuple read(const std::string& fileName) {
    ComplexTuple complex;
    ifstream file;
    file.open(fileName.c_str());
    if (!file) {
      throw std::logic_error("File not found: " + fileName);
    }


    CRef<RepSet<ElementaryCube> > repSet(new RepSet<ElementaryCube>());
    RepCubSetBuilder<RepSet<ElementaryCube> > repCubSetBuilder(repSet);

    readCubicalSet(file, repCubSetBuilder);

    ComplexConstructor constructor(repCubSetBuilder.getDim(), complex, repSet());
    boost::mpl::for_each<CubComplexDescriptors>(constructor);

    return complex;
  }

};

class SComplexReaderFromCubes {

public:
  static ComplexTuple read(const std::string& fileName) {
    ifstream file;
    file.open(fileName.c_str());
    if (!file) {
      throw std::logic_error("File not found: " + fileName);
    }

    SComplexReader<SComplexDefaultTraits> reader;
    
    return boost::make_tuple(SComplexCubesMajorId, 0, reader(fileName, 3, 1));
  }

};

template<typename ComplexReader, int MAJOR_ID>
struct ComplexReaderDescriptor {
  typedef ComplexReader type;
  static const int majorId = MAJOR_ID;
};

typedef boost::mpl::vector<ComplexReaderDescriptor<GenericCubSComplexReader, CubicalMajorId>,
			   ComplexReaderDescriptor<SComplexReaderFromCubes, SComplexCubesMajorId>
			   > ComplexReaderDescriptors;


struct ComplexReader {
  ComplexTuple& complex;
  const ComplexMajorId majorId;
  const std::string fileName;

  ComplexReader(ComplexMajorId _majorId, ComplexTuple& _complex, const std::string _fileName):
    majorId(_majorId), complex(_complex), fileName(_fileName) {
    BOOST_ASSERT(boost::get<2>(complex).empty());
  }

  template<typename Reader, int MajorId>
  void operator()(ComplexReaderDescriptor<Reader, MajorId>) {
    if (MajorId == majorId) {
      BOOST_ASSERT(boost::get<2>(complex).empty());      
      complex = Reader::read(fileName);
    }
  }

};

template<typename Derive>
class AlgorithmExecutor {
  ComplexTuple& complex;
  std::string name;

protected:
  AlgorithmExecutor(ComplexTuple& _complex, const std::string _name):
    complex(_complex), name(_name) {  
  }

public:

  template<typename Complex, ComplexMajorId MajorId, int MinorId>
  void operator()(ComplexDescriptor<Complex, MajorId, MinorId>) {
    if (boost::get<0>(complex) == MajorId && boost::get<1>(complex) == MinorId) {
      Stopwatch stopwatch;
      boost::shared_ptr<Complex>& complexPtr = boost::any_cast<boost::shared_ptr<Complex>&>(boost::get<2>(complex));
      
      std::cout << "Starting: " << name << std::endl;
      ((Derive*)this)->run(*complexPtr);
      std::cout << "Finished: " << name << " in: " << stopwatch << std::endl;
    }
  }

};

struct ShaveExecutor: AlgorithmExecutor<ShaveExecutor> {

  ShaveExecutor(ComplexTuple& _complex): AlgorithmExecutor<ShaveExecutor>(_complex, "Shave") {}

  template<typename Complex>
  void run(Complex& complex) {
    (ShaveAlgorithmFactory::createDefault(complex))();
  }

};

struct CoreductionExecutor: AlgorithmExecutor<CoreductionExecutor> {

  CoreductionExecutor(ComplexTuple& _complex): AlgorithmExecutor<CoreductionExecutor>(_complex, "Coreduction") {}

  template<typename Complex>
  void run(Complex& complex) {
    (*CoreductionAlgorithmFactory::createDefault(complex))();
  }

};

struct RFCExecutor: AlgorithmExecutor<RFCExecutor> {

  RFCExecutor(ComplexTuple& _complex): AlgorithmExecutor<RFCExecutor>(_complex, "RFC") {}

  template<typename Complex>
  void run(Complex& complex) {
    typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
    typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

    CRef<ReducibleFreeChainComplexType> RFCComplexCR = 
      (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(complex))();

    CRef<HomologySignature<int> > homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
  }

};

namespace {

  void validateOptions(const RedHomOptions& options) throw(std::logic_error) {
    switch (options.getFileType()) {
    case RedHomOptions::cubical : break;
    case RedHomOptions::simplicial : break;
    default: throw std::logic_error("Unknown file type");
    }

    if (options.getInputFilesCount() != 1) {
      throw std::logic_error("Support exactly one input file");
    }
  }

  ComplexMajorId determineComplexMajorId(const RedHomOptions& options) throw(std::logic_error) {    
    if (options.getFileType() == RedHomOptions::cubical && ! options.isGeneralSComplex()) {
      std::cout << "Using CubSComplex" << std::endl;
      return CubicalMajorId;
    } else if (options.getFileType() == RedHomOptions::cubical && options.isGeneralSComplex()) {
      std::cout << "Using SComplex for cubes" << std::endl;
      return SComplexCubesMajorId;
    }

    throw logic_error("Unsupported file reader type");
  }

}

int main(int argc, char** argv) {

  try {

    Stopwatch totalStopwatch;

    RedHomOptions options(argc, argv);

    if (options.isHelp()) {
      std::cout << options.getHelp() << std::endl;
      return 0;
    }

    validateOptions(options);

    ComplexTuple complex;
    ComplexMajorId complexMajorId = determineComplexMajorId(options);
    ComplexReader complexReader(complexMajorId, complex, options.getInputFile(0));

    Stopwatch readStopwatch;
    boost::mpl::for_each<ComplexReaderDescriptors>(complexReader);
    std::cout << "Complex read in: " << readStopwatch << std::endl;


    if (options.isRfcAtStart()) {
      boost::mpl::for_each<ComplexDescriptors>(RFCExecutor(complex));
    }

    if (options.isShave()) {
      boost::mpl::for_each<ComplexDescriptors>(ShaveExecutor(complex));
    }

    if (options.isCoreduction()) {
      boost::mpl::for_each<ComplexDescriptors>(CoreductionExecutor(complex));
    }

    if (options.isRfcAtEnd()) {
      boost::mpl::for_each<ComplexDescriptors>(RFCExecutor(complex));
    }
    
    std::cout << "Total execution: " << totalStopwatch << std::endl;
  } catch (std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "Unknown exception" << std::endl;
    return 1;
  }

  return 0;
}
