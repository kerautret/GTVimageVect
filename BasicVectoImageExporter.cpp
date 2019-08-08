#include "BasicVectoImageExporter.h"
#include <sstream> 


BasicVectoImageExporter::BasicVectoImageExporter(const std::string &imageName,
                                                 unsigned int width, unsigned int height,
                                                 bool displayMesh, double scale): myScale(scale)
{
  myWidth=width;
  myHeight=height;
  std::string ext = imageName.substr( imageName.find_last_of(".") + 1 );
  myExportType = (ext == "svg" || ext == "SVG")? ExportType::SvgExport :
                 (ext == "eps" || ext == "EPS")? ExportType::EpsExport: UnknowExport;
  myOutputStream.open(imageName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
  fillEPSHeader();
  myDisplayMesh = displayMesh;
  LINE_COLOR  = " 0.8 0.1 0.1 ";
  POINT_COLOR = " 0.1 0.1 0.8 ";
}



std::string BasicVectoImageExporter::getExportType(){
  switch (myExportType) {
    case EpsExport:
      return "eps";
      break;
    case SvgExport:
      return "svg";
      break;
    default:
      break;
  }
  return "unknow";
}


void BasicVectoImageExporter::fillHeader()
{
  switch (myExportType) {
    case EpsExport:
      fillEPSHeader();
      break;
      case SvgExport:
      fillSVGHeader();
    default:
      break;
  }
  DGtal::trace.warning() << "Header will be empty (unknow export type), change your extension file" << std::endl;
}


void BasicVectoImageExporter::fillEPSHeader()
{
  myOutputStream << "%!PS-Adobe-2.0 EPSF-2.0"
                 << "%%Title:  testBoard.eps \n"
                 << "%%Creator: BasicVectoImageExporter \n"
                 << "%%BoundingBox: 0 -1 " << (int) (myWidth*myScale)<< " " << (int)((myHeight-1)*myScale) << "\n" 
                 << "%Magnification: 30.0000\n"
                 << "%%EOF \n" 
                 << (int) myScale << " " << (int) myScale <<" scale" <<std::endl;
}




void BasicVectoImageExporter::fillSVGHeader()
  {
    
    myOutputStream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
                   << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
                   << "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
                   <<" <svg width=\"100mm\" height=\"100mm\" \n viewBox=\"0 0 " << myWidth << " " << myHeight << "\""
                   << " xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" >"
                   << "<desc>output.svg, created with DGtalTools</desc>" << std::endl;
  }

