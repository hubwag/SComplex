#ifndef _CUBSCOMPLEX_READER_HPP
#define _CUBSCOMPLEX_READER_HPP

#include "CubSComplex.hpp"
#include "../../RedHomCAPD.h"



template<int DIM>
class CubSComplexReader {

  class BmpCubSetBuilder{
    typedef RepSet<ElementaryCube> CubicalSet;
  public:
    BmpCubSetBuilder():
      fullCubes(true),
      embDim(0)
    {}

    void setFullCubes(bool b){
      // Informuje o typie kostek. Jesli b==true,
      // dostarczane sa wylacznie kostki wymiaru rownego wymiarowi wlozenia
      fullCubes=b;
    }
    void setDim(int A_embDim){
      // Dostarcza wymiar wlozenia
      embDim=A_embDim;
    }
    void addCell(int coords[]){
      // Dostarcza n wspolrzednych kostki, gdzie n to wymiar wlozenia
      // W przypadku gdy setFullCubes dostarczylo false wspolrzedne sa parzyste
      // dla przedzialow zdegenerowanych, a w przeciwnym razie sa nieparzyste
      // W przypadku gdy setFullCubes dostarczylo true dostarczane sa wspolrzedne
      // lewego dolnego rogu pelnej kostki
      if(fullCubes){
        std::vector<int> data;
        data.reserve(embDim);
        for(int i=0;i<embDim;++i){
          data.push_back(2*coords[i]+1);
        }
        cubicalSet.insert(ElementaryCube(data));
      }else{
        std::vector<int> data;
        bool* parity=new bool[embDim];
        data.reserve(embDim);
        for(int i=0;i<embDim;++i){
          data.push_back(coords[i]/2);
          parity[i]=(coords[i] % 2);
        }
        cubicalSet.insert(ElementaryCube(&data[0],parity,embDim));
      }
    }

    CubicalSet* getRepSet() {
      return new CubicalSet(cubicalSet);
    }

  private:
    bool fullCubes;
    int embDim;
    CubicalSet cubicalSet;
  };

public:

  boost::shared_ptr<CubSComplex<DIM> > operator()(std::string fileName) {
    ifstream file;
    file.open(fileName.c_str());
    if (!file) {
      std::cerr << "File not found: " << fileName << std::endl;
      return boost::shared_ptr<CubSComplex<DIM> >();
    }

    BmpCubSetBuilder csb;
    readCubicalSet(file,csb);

    return boost::shared_ptr<CubSComplex<DIM> >(new CubSComplex<DIM>(*csb.getRepSet()));
  }

  
};


#endif //  _CUBSCOMPLEX_READER_HPP

