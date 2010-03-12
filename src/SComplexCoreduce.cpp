#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <set>
#include <map>
#include <cassert>
// #include <array>
#include <list>
using namespace std;

#include "simple_set.h"

#include "Simplex.hpp"
#include "SimplexSComplex.hpp"

#include <iostream>
#include <queue>
#include <deque>

using namespace std;
#include <capd/auxil/Stopwatch.h>
#include <capd/auxil/CRef.h>

#include <capd/homologicalAlgebra/embeddingDim.h>

#include <capd/vectalg/MatrixSlice.h>
#include <capd/matrixAlgorithms/intMatrixAlgorithms.hpp>

#include <capd/homologicalAlgebra/homologicalAlgebra.hpp>
#include <capd/homologicalAlgebra/homAlgFunctors.hpp>
#include <capd/homologicalAlgebra/cubSetFunctors.hpp>
#include <capd/homologicalAlgebra/ReducibleFreeChainComplex.hpp>


ofstreamcout fcout; // ?

#include "SComplexAlgs.hpp"

typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

#include <SComplex.hpp>
#include <SComplexDefaultTraits.hpp>
#include <SComplexAlgs.hpp>
#include <SComplexBuilderFromSimplices.hpp>

#include "SimplexSubdivision.hpp"
#include <CrHomS.hpp>


void CrHomS_torus(int argc,char* argv[])
{
    Stopwatch swTot;
	 vector<set<int> > tris = makeTest();

	 Stopwatch swBuild;
	typedef SComplex<SComplexDefaultTraits> Complex;
	SComplexBuilderFromSimplices<long, SComplexDefaultTraits> builder(1234567);

	boost::shared_ptr<Complex> complex = builder(tris, 3, 1);
	cout << " --- built in " << swBuild << std::endl;
	testReduce(*complex);
    cout << " --- generated simplicial complex --- \n cardinality: " << complex->cardinality() << endl;
}

int main(int argc,char* argv[])
{
    try
    {
		CrHomS_torus(argc,argv);
    }
    catch (std::exception& e)
    {
        std::cout << "Caught exception: " << e.what() << endl;
    }
    catch (std::string& s)
    {
        std::cout << "Caught exception: " << s.c_str() << endl;
    }
    catch (const char* c)
    {
        std::cout << "Caught exception: " << c << endl;
    }
    catch (...)
    {
        std::cout << "Caught an unknown exception: " << endl;
    }
}
