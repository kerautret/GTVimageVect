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
  BasicVectoImageExporter exporter ("testExport.eps", 50, 50, true, 2.0);
  BasicVectoImageExporter exporterSVG ("testExport.svg", 50, 50, true, 2.0);
  // test auto detected extension:
  trace.info() << "Export type (testExport.eps) : " << exporter.getExportType() << std::endl;
  trace.info() << "Export type (testExport.svg) : " << exporterSVG.getExportType() << std::endl;
  

  
  std::vector<Z2i::RealPoint> contour = {Z2i::RealPoint(0,2), Z2i::RealPoint(20, 2),
                                  Z2i::RealPoint(18,18), Z2i::RealPoint(2,18)};
  // to be considered as hole we need to drow it on opposite direction
  std::vector<Z2i::RealPoint> hole = {Z2i::RealPoint(5,15), Z2i::RealPoint(15,15),
                                      Z2i::RealPoint(15,  5),Z2i::RealPoint(5,5)};
  
  std::vector<Z2i::RealPoint> region = {Z2i::RealPoint(25,3), Z2i::RealPoint(30, 3),
    Z2i::RealPoint(40,25)};
  std::vector<Z2i::RealPoint> region2 = {Z2i::RealPoint(25,13), Z2i::RealPoint(30, 13),
    Z2i::RealPoint(40,35)};
  std::vector<Z2i::RealPoint> region3 = {Z2i::RealPoint(35,20), Z2i::RealPoint(43, 15),
    Z2i::RealPoint(48,20),Z2i::RealPoint(45,30)};

  std::vector<std::vector<Z2i::RealPoint>> regions;
  regions.push_back(region2);
  regions.push_back(region3);

  auto contoursHole = std::vector<std::vector<Z2i::RealPoint>>();
  contoursHole.push_back(hole);
  
  // test contour points
  exporter.addContourPoints(contour);
  exporterSVG.addContourPoints(contour);
  
  // test regions with holes
  exporter.addRegionWithHoles(contour, contoursHole, DGtal::Color::Blue);
  exporterSVG.addRegionWithHoles(contour, contoursHole, DGtal::Color::Blue);

  // test single contours:
  exporter.addContour(contoursHole[0], DGtal::Color::Yellow);
  exporterSVG.addContour(contoursHole[0], DGtal::Color::Yellow);

  // test add region
  exporter.addRegion(region, DGtal::Color::Green, 0.0);
  exporterSVG.addRegion(region, DGtal::Color::Green, 0.0);

  // test add regions
  exporter.addRegions(regions, DGtal::Color::Purple);
  exporterSVG.addRegions(regions, DGtal::Color::Purple);

  

  exporter.closeFigure();
  exporterSVG.closeFigure();

  return 0;
}


