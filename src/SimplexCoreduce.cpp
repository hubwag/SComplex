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

	cout << "\n\ntesting on morse-smale complex\n\n";
	 Stopwatch swComp,swRed;
	 //(ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    //cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;


	 CRef<ReducibleFreeChainComplexType> RFCComplexCR=
		(ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SComplex, ReducibleFreeChainComplexType>(complex))();
	 cout << " --- RFCC constructed  " << endl;

	 CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
	 cout << " --- Computation completed in " << swComp  << std::endl;
	 cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;
}

#include "SComplex.hpp"
#include "SComplexDefaultTraits.hpp"

template<typename SComplex>
void CrHomS()
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
    (ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    // Stopwatch swCoRed;
    // (CoreductionAlgorithmFactory::createDefault(SComplexCR()))();
    // cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

/*
    CRef<ReducibleFreeChainComplexType> RFCComplexCR=
        (ReducibleFreeChainComplexOverZFromSComplexAlgorithm<SimplexSComplex, ReducibleFreeChainComplexType>(SComplexCR()))();
    cout << " --- RFCC constructed  " << endl;

    CRef<HomologySignature> homSignCR=HomAlgFunctors<FreeModuleType>::homSignViaAR_Random(RFCComplexCR);
    cout << " --- Computation completed in " << swComp  << std::endl;
    cout << " --- Computed homology is: \n\n" << homSignCR()  << std::endl;
    cout << " --- Total computation time is: " << swTot  << std::endl;

    cout << "\n\n\n";*/

    boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy<SimplexSComplex> > >
    cored = CoreductionAlgorithmFactory::createDefault(SComplexCR());

    (*cored)();


    boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy< ::SComplex<SComplexDefaultTraits> > > >
    cored2 = CoreductionAlgorithmFactory::createDefault(*cored->getStrategy()->outputComplex);

    (*cored2)();

    cout << *cored2->getStrategy()->outputComplex.size() << endl;

//    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

    testReduce(*cored2->getStrategy()->outputComplex);
    delete cored->getStrategy()->outputComplex;
    delete cored2->getStrategy()->outputComplex;
}

#include "simplexIO.hpp"

#include "OldCored.hpp"

void showObj(const string &s)
{
	Stopwatch swTot;
	ifstream ifs(s.c_str());
	CRef<SimplexSComplex> SComplexCR(new SimplexSComplex());

	if (s.find(".obj") != string::npos)
	{
		cerr << "parsing obj simplex file" << endl;
		parseObj(ifs, SComplexCR());
	}
	else {
		cout << "parsing dat or plain txt simplex file" << endl;
		parseDat(ifs, SComplexCR());
	}

	// SComplexCR() = subdivide3(SComplexCR());

	cout << "parsed file, cardinality: " << SComplexCR().cardinality() << endl;
	cout << "it took: " << swTot << endl;

	Stopwatch swComp,swRed;
	//(ShaveAlgorithmFactory::createDefault(SComplexCR()))();
    // cout << " --- Shave reduced the size to " << SComplexCR().cardinality() << " in " << swRed <<  endl;

    Stopwatch swCoRed;

     //boost::shared_ptr< CoreductionAlgorithm<OldReduceStrategy<SComplex> > > old_red (new CoreductionAlgorithm<OldReduceStrategy<SComplex> >(new OldReduceStrategy<SComplex>(SComplexCR())));

boost::shared_ptr<CoreductionAlgorithm<OldReduceStrategy<SimplexSComplex> > >
    old = OldCoreductionAlgorithmFactory::createDefault(SComplexCR());

	(*old)();

	cout << "~DOBRE HOMOLOGIE PO ZWYKLYCH KOREDUKCJACH: \n";
    testReduce(SComplexCR());


    boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy<SimplexSComplex> > >
    cored = CoreductionAlgorithmFactory::createDefault(SComplexCR());

    (*cored)();

	cout << "DOBRE HOMOLOGIE PO MORSIE: \n";
    testReduce(*cored->getStrategy()->outputComplex);

	boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy< ::SComplex<SComplexDefaultTraits> > > >
    cored2 = CoreductionAlgorithmFactory::createDefault(*cored->getStrategy()->outputComplex);

    (*cored2)();

    boost::shared_ptr<CoreductionAlgorithm<DefaultReduceStrategy< ::SComplex<SComplexDefaultTraits> > > >
    cored3 = CoreductionAlgorithmFactory::createDefault(*cored2->getStrategy()->outputComplex);

    (*cored3)();

    // cout << *cored3->getStrategy()->outputComplex.cardinality() << endl;

//    cout << " --- Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;




    // cout << " --- AKQ Coreduction reduced the size to " << SComplexCR().cardinality() << " in " << swCoRed <<  endl;

    testReduce(*cored3->getStrategy()->outputComplex);
    delete cored->getStrategy()->outputComplex;
    delete cored2->getStrategy()->outputComplex;
    delete cored3->getStrategy()->outputComplex;
}

int main(int n, char **v)
{
	if (n <= 1)
	{
		cerr << "nie ma argumentow!!";
		exit(1);
	}


/*	string files[] = {"bing.dat",
"bjorner.dat",
"c-ns.dat",
"c-ns2.dat",
"c-ns3.dat",
"dunce.dat",
"dunce_hat.dat",
"gruenbaum.dat",
"knot.dat",
"lockeberg.dat",
"mani-walkup-C.dat",
"mani-walkup-D.dat",
"nc_sphere.dat",
"nonextend.dat",
"nonpl_sphere.dat",
"poincare.dat",
"projective.dat",
"rudin.dat",
"simon.dat",
"simon2.dat",
"solid_2_torus.dat",
"ziegler.dat"};
*/

	// string pref = "c:/Users/hub/Downloads/library/";
	// string pref = "c:/Users/hub/Downloads/";
	// string files[] = {"knot.dat"};

	// string files[] = {"poincare.dat"};

	// string files[] = {"first.txt", "second.txt", "third.txt", "fourth.txt", "fifth.txt"};

	string files[] = {"buddha.obj", "bunny.obj", "dragon.obj", "toruses.obj", "s2.obj"};

	//for (int i  = 0; i < sizeof(files)/sizeof(*files); i++)
	{
		//cout << files[i] << endl;
		showObj(v[1]);
		cout << endl;
	}
}
using namespace std;


int __main(int argc,char* argv[])
{
   vector<set<int> > torus = makeSpaceFromWelds(makeProjectiveSpaceWelds());
   // vector<set<int> > torus;

   // for (int i = 0; i < 4; i++)
   //	torus = subdivide6(torus);
   CrHomS_fromTris<SimplexSComplex>(torus, "torus");
   // showObj("c:/users/hub/downloads/bunny.obj");
}


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
