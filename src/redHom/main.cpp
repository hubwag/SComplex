#include <boost/algorithm/string/join.hpp>
#include <RedHomOptions.hpp>

#include <iostream>


int main(int argc, char** argv) {

  try {

    RedHomOptions options(argc, argv);

    if (options.isHelp()) {
      std::cout << options.getHelp() << std::endl;
    }

    switch (options.getFileType()) {
    case RedHomOptions::cubical : std::cout << "File type set to cubical" << std::endl; break;
    case RedHomOptions::simplicial : std::cout << "File type set to simplicial" << std::endl; break;
    default: throw std::logic_error("Unknown file type");
    }
  
    std::cout << "Input files: " << boost::algorithm::join(options.getInputFiles(), " ") << std::endl;
    
    if (options.isGeneralSComplex()) {
      std::cout << " Using general SComplex implementation" << std::endl;
    }
    
  } catch (std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "Unknown exception" << std::endl;
    return 1;
  }

  return 0;
}
