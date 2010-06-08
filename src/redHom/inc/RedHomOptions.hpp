#ifndef _RED_HOM_OPTIONS_HPP
#define _RED_HOM_OPTIONS_HPP

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

#include <sstream>
#include <string>


class RedHomOptions {
public:
  enum FileType { cubical=0, simplicial };


public:



  RedHomOptions(int ac, char** av, int inputFilesCount=3): desc("Allowed options"), 
							   inputFiles(inputFilesCount) {

    using namespace boost::program_options;
    
    desc.add_options()
      ("help,h", "produce help message")
      ;

    options_description fileDesc;
    fileDesc.add_options()
      ("fileType,f", value<FileType>(&fileType),
       "input file type: c=cubical, s=simplicial");
    
    positional_options_description positionalFileDesc;

    for (int i = 0; i < inputFilesCount; ++i) {
      std::string key = "input-file" + boost::lexical_cast<std::string>(i);
      positionalFileDesc.add(key.c_str(), 1);
      fileDesc.add_options()(key.c_str(), value<std::string>(&inputFiles[i]), "input file");
    }

    desc.add(fileDesc);

    options_description complexDesc;
    complexDesc.add_options()
      ("general-scomplex,g", bool_switch(&generalSComplex), "use general SComplex");

    desc.add(complexDesc);

    store(command_line_parser(ac, av).options(desc).positional(positionalFileDesc).run(), vm);
    notify(vm);    

  }
 
  bool isHelp() const {
    return vm.count("help");
  }

  std::string getHelp() const {
    std::ostringstream helpStream;
    helpStream << desc;

    return helpStream.str();
  }

  FileType getFileType() const {
    return fileType;
  }

  const std::string& getInputFile(int i) const {
    return inputFiles[i];
  }

  const std::vector<std::string> getInputFiles() const {
    return inputFiles;
  }

  bool isGeneralSComplex() const {
    return generalSComplex;
  }

private:

  boost::program_options::variables_map vm;
  boost::program_options::options_description desc;
  std::vector<std::string> inputFiles;
  bool generalSComplex;

  FileType fileType;
};


void validate(boost::any& v, 
              const std::vector<std::string>& values,
              RedHomOptions::FileType*, int)
{
    using namespace boost::program_options;

    // Make sure no previous assignment to 'a' was made.
    validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = validators::get_single_string(values);

    // Do regex match and convert the interesting part to 
    // int.
    if (s == "s" || s == "simplicial") {
      v = boost::any(RedHomOptions::simplicial);
    } else if (s == "c" || s == "cubical") {
      v = boost::any(RedHomOptions::cubical);
    } else {
        throw validation_error(validation_error::invalid_option_value);
    }        
}


#endif //  _RED_HOM_OPTIONS_HPP
