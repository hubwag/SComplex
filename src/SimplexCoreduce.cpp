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

/*
template<typename simplex_t, typename out_iter_t>
void get_neighbours(simplex_t &s, out_iter_t out)
{
    for (auto it = s.coborder.begin(), end = s.coborder.end(); it != end; ++it)
    {
        for (auto it2 = (*it)->border.begin(), end2 = (*it)->border.end(); it2 != end2; ++it2)
            *out++ = *it2;
    }
}
*/

/****  CRHomS.cpp            ****/
/**** (c) Marian Mrozek 2009 ****/

#include <iostream>
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

//#include "CubSComplex.hpp"
#include "SComplexAlgs.hpp"

typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

template<typename SComplex>
void CrHomS(int argc,char* argv[])
{
    Stopwatch swTot;

    CRef<SComplex> SComplexCR(new SComplex());

    int n = 100000;
    int mod = n; // = 4
    unsigned d = 4;

    for (int i = 0; i < n; i++)
    {
        set<int> s;
        while (s.size() < d)
        {
            s.insert(rand()%mod);
        }

        SComplexCR().add_simplex(s);
    }

/*
    add(SComplexCR(), 0, 1, 2)->debug_output__();
    add(SComplexCR(), 1, 3);
    add(SComplexCR(), 2, 3);
    add(SComplexCR(), 3, 4);
    add(SComplexCR(), 2, 4);
    add(SComplexCR(), 3, 5);

*/

/*
	add(SComplexCR(), 0, 1);
	add(SComplexCR(), 0, 3);
	add(SComplexCR(), 0, 4);
	add(SComplexCR(), 1, 3);
	add(SComplexCR(), 3, 4);
*/

    cout << " --- generated random simplicial complex --- \n cardinality: " << SComplexCR().cardinality() << endl;

//  SComplexAlgs<CubSComplex>::test(SComplexCR());

    Stopwatch swComp,swRed;
    (ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    Stopwatch swCoRed;
    (CoreductionAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

/*
      CRef<ReducibleFreeChainComplexType> RFCComplexCR=
    		(ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType,ElementaryCellType>(SComplexCR()))();
      cout << " --- RFCC constructed  " << endl;

      CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
      cout << " --- Computation completed in " << swComp  << std::endl;
      cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;

      */

      cout << " --- Total computation time is: " << swTot  << std::endl;

}
ofstreamcout fcout;

int main(int argc,char* argv[])
{
    std::string outFname="out.txt";
    fcout.turnOn();
    fcout.open(outFname.c_str());

    try
    {
        CrHomS<SimplexSComplex>(argc,argv);
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
