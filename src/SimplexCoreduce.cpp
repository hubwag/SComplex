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

#include "SimplexSubdivision.hpp"
#include <CrHomS.hpp>

template<typename SComplex>
void CrHomS_torus(int argc,char* argv[])
{
    Stopwatch swTot;
    CRef<SComplex> SComplexCR(new SComplex());
	 vector<set<int> > tris = makeTest();

    for (size_t i = 0; i < tris.size(); i++)
    {
    	SComplexCR().addSimplex(tris[i]);
    }

	 testReduce(SComplexCR());
    cout << " --- generated simplicial complex --- \n cardinality: " << SComplexCR().cardinality() << endl;
}

int main(int argc,char* argv[])
{
    try
    {
		CrHomS_torus<SimplexSComplex>(argc,argv);
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
