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


enum ComplexMajorId {SComplexMajorId = 0, CubicalMajorId = 1};

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
				 ComplexDescriptor<SComplex<SComplexDefaultTraits>, SComplexMajorId, 0>
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

template<typename ComplexReader, int MAJOR_ID>
struct ComplexReaderDescriptor {
  typedef ComplexReader type;
  static const int majorId = MAJOR_ID;
};

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

typedef boost::mpl::vector<ComplexReaderDescriptor<GenericCubSComplexReader, CubicalMajorId>
			   > ComplexReaderDescriptors;


struct CoreductionExecutor {
  ComplexTuple& complex;

  CoreductionExecutor(ComplexTuple& _complex):
    complex(_complex) {  
  }

  template<typename Complex, ComplexMajorId MajorId, int MinorId>
  void operator()(ComplexDescriptor<Complex, MajorId, MinorId>) {
    if (boost::get<0>(complex) == MajorId && boost::get<1>(complex) == MinorId) {
      boost::shared_ptr<Complex>& complexPtr = boost::any_cast<boost::shared_ptr<Complex>&>(boost::get<2>(complex));
      
      std::cout << "Starting correduction" << std::endl;
      (*CoreductionAlgorithmFactory::createDefault(*complexPtr))();
      std::cout << "Correduction algorithm finished" << std::endl;

//       BOOST_MESSAGE("Building RFC complex");
//       CRef<ReducibleFreeChainComplexType> RFCComplexCR=
// 	(ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex<TraitsT>, ReducibleFreeChainComplexType>(complex))();

//       //std::cerr << RFCComplexCR() << std::endl;

//       BOOST_MESSAGE("Computing homology signature");
//       CRef<HomologySignature<int> > homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

//       std::ostringstream signature;
//       signature << homSignCR();

//       std::string sig = signature.str();
//       std::replace(sig.begin(), sig.end(), '\n', '#');

    }
  }

};

struct RFCExecutor {
  ComplexTuple& complex;

  RFCExecutor(ComplexTuple& _complex):
    complex(_complex) {  
  }

  template<typename Complex, ComplexMajorId MajorId, int MinorId>
  void operator()(ComplexDescriptor<Complex, MajorId, MinorId>) {
    typedef ElementaryCell ElementaryCellType;
    typedef int ScalarType;
    typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
    typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
    typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

    if (boost::get<0>(complex) == MajorId && boost::get<1>(complex) == MinorId) {
      boost::shared_ptr<Complex>& complexPtr = boost::any_cast<boost::shared_ptr<Complex>&>(boost::get<2>(complex));
      
      std::cout << ("Building RFC complex") << std::endl;
      CRef<ReducibleFreeChainComplexType> RFCComplexCR = 
	(ReducibleFreeChainComplexOverZFromSComplexAlgorithm<Complex, ReducibleFreeChainComplexType>(*complexPtr))();

      std::cout << ("Computing homology signature") << std::endl;;
      CRef<HomologySignature<int> > homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);

      std::cout << homSignCR() << std::endl;
    }
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
    if (options.getFileType() == RedHomOptions::cubical) {
      return CubicalMajorId;
    } 

    throw logic_error("Unsupported file reader type");
  }

}

int main(int argc, char** argv) {

  try {

    RedHomOptions options(argc, argv);

    if (options.isHelp()) {
      std::cout << options.getHelp() << std::endl;
      return 0;
    }

    validateOptions(options);

    ComplexTuple complex;
    ComplexMajorId complexMajorId = determineComplexMajorId(options);
    ComplexReader complexReader(complexMajorId, complex, options.getInputFile(0));
    boost::mpl::for_each<ComplexReaderDescriptors>(complexReader);

    boost::mpl::for_each<ComplexDescriptors>(CoreductionExecutor(complex));
    boost::mpl::for_each<ComplexDescriptors>(RFCExecutor(complex));
    
  } catch (std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "Unknown exception" << std::endl;
    return 1;
  }

  return 0;
}
