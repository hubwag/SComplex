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

//#include "CubSComplex.hpp"

ofstreamcout fcout; // ?

#include "SComplexAlgs.hpp"

// typedef ElementaryCell ElementaryCellType;
typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

#include "SimplexSubdivision.hpp"

template<typename SComplex>
void CrHomS_fromTris(const vector<set<int> > &tris, const string &description)
{
    Stopwatch swTot;

    CRef<SComplex> SComplexCR(new SComplex());

    for (size_t i = 0; i < tris.size(); i++)
    {
        SComplexCR().addSimplex(tris[i]);
    }

    cout << " --- generated " << description << " simplicial complex --- \n cardinality: " << SComplexCR().cardinality() << endl;

    if (SComplexCR().cardinality() < 1000)
    {
        Stopwatch swAlgebra;
        CRef<ReducibleFreeChainComplexType> RFCComplexCR=
            (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
        cout << " --- RFCC constructed  " << endl;


        CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
        cout << " --- Computation completed in " << swAlgebra << std::endl;
        cout << " --- Computed homology is: \n\n" << homSignCR() << std::endl;

    }

    //  SComplexAlgs<CubSComplex>::test(SComplexCR());

    Stopwatch swComp,swRed;
    (ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    Stopwatch swCoRed;
    (CoreductionAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

    CRef<ReducibleFreeChainComplexType> RFCComplexCR=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
    cout << " --- RFCC constructed  " << endl;

    /*cout << "debug " << endl;

    for (SimplexSComplex::iterator it = SComplexCR().all_begin(), end = SComplexCR().all_end(); it != end; ++it)
    {
        Simplex &s = (Simplex&)*it;
        cout << distance(s.border_begin(1), s.border_end(1)) << endl;
        cout << distance(s.coborder_begin(1), s.coborder_end(1)) << endl;
    }
    */

    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
    cout << " --- Computation completed in " << swComp  << std::endl;
    cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;

    cout << " --- Total computation time is: " << swTot  << std::endl;

    cout << "\n\n\n";
}

template<typename SComplex>
void CrHomS(int argc,char* argv[])
{
    Stopwatch swTot;

    CRef<SComplex> SComplexCR(new SComplex());

    int n = 1000;
    int mod = n; // = 4
    unsigned d = 4;

    for (int i = 0; i < n; i++)
    {
        set<int> s;
        while (s.size() < d)
        {
            s.insert(rand()%mod);
        }

        SComplexCR().addSimplex(s);
    }

    cout << " --- generated random simplicial complex --- \n cardinality: " << SComplexCR().cardinality() << endl;

//  SComplexAlgs<CubSComplex>::test(SComplexCR());

    Stopwatch swComp,swRed;
    // (ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    // cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    // Stopwatch swCoRed;
    // (CoreductionAlgorithmFactory::createDefault(SComplexCR()))();
    // cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;


    CRef<ReducibleFreeChainComplexType> RFCComplexCR=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
    cout << " --- RFCC constructed  " << endl;

    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
    cout << " --- Computation completed in " << swComp  << std::endl;
    cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;
    cout << " --- Total computation time is: " << swTot  << std::endl;

    cout << "\n\n\n";
}

#include "simplexIO.hpp"

void showObj(const string &s)
{
	Stopwatch swTot;
	ifstream ifs(s.c_str());
	CRef<SimplexSComplex> SComplexCR(new SimplexSComplex());

	parseObj(ifs, SComplexCR()); //

	cout << SComplexCR().cardinality() << endl;

	Stopwatch swComp,swRed;
	(ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    Stopwatch swCoRed;
    (CoreductionAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;


    CRef<ReducibleFreeChainComplexType> RFCComplexCR=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
    cout << " --- RFCC constructed  " << endl;

    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
    cout << " --- Computation completed in " << swComp  << std::endl;
    cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;
    cout << " --- Total computation time is: " << swTot  << std::endl;

    cout << "\n\n\n";
}

int main()
{
	string pref = "c:/Users/hub/Downloads/";
	string files[] = {"buddha.obj", "bunny.obj", "dragon.obj", "toruses.obj", "s2.obj"};

	for (int i  = 0; i < sizeof(files)/sizeof(*files); i++)
	{
		cout << files[i] << endl;
		showObj(pref + files[i]);
		cout << endl;
	}
}

int __main(int argc,char* argv[])
{
    std::string outFname="out.txt";
    fcout.turnOn();
    fcout.open(outFname.c_str());

    vector<set<int> > klein = makeSpaceFromWelds(makeKleinWelds());
    vector<set<int> > torus = makeSpaceFromWelds(makeTorusWelds());
    vector<set<int> > pspace = makeSpaceFromWelds(makeProjectiveSpaceWelds());

    vector<set<int> > random;
    int n = 1000;
    int mod = n; // = 4
    unsigned d = 4;

    for (int i = 0; i < n; i++)
    {
        set<int> s;
        while (s.size() < d)
        {
            s.insert(rand()%mod);
        }

        random.push_back(s);
    }

    vector<set<int> > small(1);
    small.back().insert(0);
    small.back().insert(1);
    small.back().insert(2);


    vector<set<int> > ex;
    ex.push_back(make_int_set(1,2,3));
    ex.push_back(make_int_set(0,1));
    ex.push_back(make_int_set(0,3));
    ex.push_back(make_int_set(0,4));
    ex.push_back(make_int_set(3,4));
    ex.push_back(make_int_set(0,5));

    try
    {
        CrHomS_fromTris<SimplexSComplex>(ex, "example from the article");

        CrHomS_fromTris<SimplexSComplex>(small, "one triangle");
        CrHomS_fromTris<SimplexSComplex>(klein, "klein bottle");
        CrHomS_fromTris<SimplexSComplex>(pspace, "real projective space");

        for (int i = 0; i < 5; i++)
        {
            klein = subdivide6(klein);
            CrHomS_fromTris<SimplexSComplex>(klein, "klein bottle (after barycntric subdivisions)");
        }

        CrHomS_fromTris<SimplexSComplex>(torus, "torus");

        for (int i = 0; i < 5; i++)
        {
            torus = subdivide3(torus);
            CrHomS_fromTris<SimplexSComplex>(torus, "torus (after 3-subdivisions)");
        }

        CrHomS_fromTris<SimplexSComplex>(random, "random complex");

        // CrHomS<SimplexSComplex>(argc,argv);
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
