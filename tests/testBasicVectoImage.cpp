#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/base/Common.h>
#include <DGtal/base/ConstAlias.h>
#include <DGtal/base/Alias.h>
#include <DGtal/helpers/StdDefs.h>

#include "BasicVectoImageExporter.h"

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////
void toto()
{
  
}

int main( int argc, char** argv )
{
  using namespace DGtal;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("eps,e", "test eps toy image generatipn.")
    ;
  
  bool parseOK = true;
  po::variables_map vm;
  try {
    po::store( po::parse_command_line(argc, argv, general_opt), vm );  
  } catch ( const std::exception& ex ) {
    parseOK = false;
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
  }
    
  po::notify(vm);    
  if( ! parseOK || vm.count("help") )
    {
      trace.info()<< "Test export image with Class BasicVectoImage." <<std::endl << "Basic usage: " << std::endl
		  << "\ttestBasicVectoImage [options] -e "<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  using namespace std;
  using namespace DGtal;
  BasicVectoImageExporter exporter ("testExport.eps", 20, 20, true, 2.0);
  BasicVectoImageExporter exporterSVG ("testExport.svg", 20, 20, true, 2.0);
  // test auto detected extension:
  trace.info() << "Export type (testExport.eps) : " << exporter.getExportType() << std::endl;
  trace.info() << "Export type (testExport.svg) : " << exporterSVG.getExportType() << std::endl;
  

  
  std::vector<Z2i::RealPoint> contour = {Z2i::RealPoint(2,2), Z2i::RealPoint(18, 2),
                                  Z2i::RealPoint(18,18), Z2i::RealPoint(2,18)};
  // to be considered as hole we need to drow it on opposite direction
  std::vector<Z2i::RealPoint> hole = {Z2i::RealPoint(5,15), Z2i::RealPoint(15,15),
                                      Z2i::RealPoint(15,  5),Z2i::RealPoint(5,5)};
  
  auto contoursHole = std::vector<std::vector<Z2i::RealPoint>>();
  contoursHole.push_back(hole);
  exporter.addContourPoints(contour);
  exporterSVG.addContourPoints(contour);
  
  exporter.addRegionWithHoles(contour, contoursHole, DGtal::Color::Blue);
  exporterSVG.addRegionWithHoles(contour, contoursHole, DGtal::Color::Blue);
  exporter.closeFigure();
  exporterSVG.closeFigure();

  return 0;
}


