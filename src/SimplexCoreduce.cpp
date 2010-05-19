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
#include <queue>
#include <deque>

using namespace std;

#include "simple_set.h"

#include "Simplex.hpp"
#include "SimplexSComplex.hpp"



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

typedef int ScalarType;
typedef FreeModule<int,capd::vectalg::Matrix<int,0,0> > FreeModuleType;
typedef FreeChainComplex<FreeModuleType> FreeChainComplexType;
typedef ReducibleFreeChainComplex<FreeModuleType,int> ReducibleFreeChainComplexType;

#include "SComplexAlgs.hpp"


#include "SimplexSubdivision.hpp"
// #include "CrHomS.hpp"

template<typename SComplex>
void CrHomS_fromTris(const vector<set<int> > &tris, const string &description)
{
    Stopwatch swTot;
    CRef<SComplex> SComplexCR(new SComplex());
	// vector<set<int> > tris = makeTest();

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

    /*Stopwatch swComp,swRed;
    (ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    Stopwatch swCoRed;
    (*CoreductionAlgorithmFactory::createDefault(SComplexCR()))();

    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

    CRef<ReducibleFreeChainComplexType> RFCComplexCR=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
    cout << " --- RFCC constructed  " << endl;

    */

    /*cout << "debug " << endl;

    for (SimplexSComplex::iterator it = SComplexCR().all_begin(), end = SComplexCR().all_end(); it != end; ++it)
    {
        Simplex &s = (Simplex&)*it;
        cout << distance(s.border_begin(1), s.border_end(1)) << endl;
        cout << distance(s.coborder_begin(1), s.coborder_end(1)) << endl;
    }
    */

    /*

    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
    cout << " --- Computation completed in " << swComp  << std::endl;
    cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;

    cout << " --- Total computation time is: " << swTot  << std::endl;

    cout << "\n\n\n"; */


    boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy<SimplexSComplex> > >
    cored = CoreductionAlgorithmFactory::createDefault(SComplexCR());

    (*cored)();

//    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

    testReduce(*cored->getStrategy()->outputComplex);
    delete cored->getStrategy()->outputComplex;
}


template<typename SComplex>
void testReduce(SComplex& complex) {

	 Stopwatch swComp,swRed;
	 //(ShaveAlgorithmFactory::createDefault(SComplexCR()))();
     //cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

	 CRef<ReducibleFreeChainComplexType> RFCComplexCR=
		(ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex, ReducibleFreeChainComplexType>(complex))();
	 cout << " --- RFCC constructed  " << endl;

	 CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR(RFCComplexCR);
	 cout << " --- Computation completed in " << swComp  << std::endl;
	 cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;
}

#include "SComplex.hpp"
#include "SComplexDefaultTraits.hpp"

#include "simplexIO.hpp"

#include "OldCored.hpp"

void showObj(const string &s, const string &method = "KMS", int subdivs = 0)
{
	Stopwatch swTot;
	ifstream ifs(s.c_str());
	CRef<SimplexSComplex> SComplexCR(new SimplexSComplex());

	if (s.find(".obj") != string::npos)
	{
		cerr << "parsing obj simplex file" << endl;
		parseObj(ifs, SComplexCR(), subdivs);
	}
	else {
		cout << "parsing dat or plain txt simplex file" << endl;
		parseDat(ifs, SComplexCR(), subdivs);
	}

	cout << "parsed file, cardinality: " << SComplexCR().cardinality() << endl;
	cout << "it took: " << swTot << endl;

	Stopwatch swComp;

	if (method == "KMS")
	{
		cout << "RUNNING KMS: \n";
		testReduce(SComplexCR());

		cout << "calculations completed in: " << swComp << endl;
		return;
	} else if (method == "CORED") {
		cout << "RUNNING STANDARD REDUCTIONS THEN KMS";
		boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy<SimplexSComplex> > >
		old = CoreductionAlgorithmFactory::createDefault(SComplexCR());
		(*old)();

		testReduce(SComplexCR());

		cout << "calculations completed in: " << swComp << endl;
		return;
	}

    boost::shared_ptr<CoreductionAlgorithm<AKQReduceStrategy<SimplexSComplex> > >
    cored = CoreductionAlgorithmFactory::createAKQ(SComplexCR());

    (*cored)();

    cout << "AKQ completed in: " << swComp << endl;

	cout << "HOMOLOGIE PO MORSIE: \n";
    testReduce(*cored->getStrategy()->outputComplex);

    cout << "calculations completed in: " << swComp << endl;

    delete cored->getStrategy()->outputComplex;
}

int main(int n, char **v)
{
	ios_base::sync_with_stdio(false);

	string method = "AKQ";
	// string file = "c:/users/hub/Downloads/buddha.obj";
	// string file = "c:/users/hub/Downloads/library/knot.dat";
	string file;

	int subdivs = 0;
	if (n <= 1)
	{
		// cerr << "using default file - probably nonexistent...\n";
		cerr << "usage: filename [AKQ|CORED|MKS] [number_of_subdivisions]" << endl;
		cerr << "no arguments - quitting\n";
		exit(1);
	}

	if (n > 1)
 	{
		file = v[1];
	}

	if (n >= 3)
	{
		cerr << "using the supplied method!\n";
		method = v[2];
	}

	if (n >= 4)
	{
		cerr << "using the supplied number of subdivisions (simple, not barycentric for now)!\n";
		subdivs = atoi(v[3]);
	}

	Stopwatch all;

	showObj(file.c_str(), method, subdivs);
	cout << "total time: " << all  << endl;

	return 0;
}
using namespace std;

/*
int main(int argc,char* argv[])
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

    CrHomS_fromTris<SimplexSComplex>(klein, "klein bottle (after barycntric subdivisions)");

    try
    {
        CrHomS_fromTris<SimplexSComplex>(ex, "example from the article");

        CrHomS_fromTris<SimplexSComplex>(small, "one triangle");
        CrHomS_fromTris<SimplexSComplex>(klein, "klein bottle");
        CrHomS_fromTris<SimplexSComplex>(pspace, "real projective space");

        for (int i = 0; i < 4; i++)
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

*/
