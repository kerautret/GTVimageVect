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
  fillHeader();
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
      break;
    default:
      DGtal::trace.warning() << "Header will be empty (unknow export type), change your extension file" << std::endl;
      break;
  }
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
    
    myOutputStream << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl;
    myOutputStream << "<!-- Created BasicImageExporter --> " << std::endl;
    myOutputStream << " " << std::endl;
    myOutputStream << "<svg" << std::endl;
    myOutputStream << "xmlns:dc=\"http://purl.org/dc/elements/1.1/\"" << std::endl;
    myOutputStream << "xmlns:cc=\"http://creativecommons.org/ns#\"" << std::endl;
    myOutputStream << "xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"" << std::endl;
    myOutputStream << "xmlns:svg=\"http://www.w3.org/2000/svg\"" << std::endl;
    myOutputStream << "xmlns=\"http://www.w3.org/2000/svg\""<< std::endl;
    myOutputStream << "xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\""<< std::endl;
    myOutputStream << "xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"" << std::endl;
    myOutputStream << "width=\""<<myWidth<<"mm\""<< std::endl;
    myOutputStream << "height=\""<<myHeight<<"mm\""<< std::endl;
    myOutputStream << "viewBox=\"0 0"<<myWidth<< " " << myHeight<< "\"" << std::endl;
    myOutputStream << "version=\"1.1\"" << std::endl;
    myOutputStream << "id=\"svg3\"" << std::endl;
    myOutputStream << "sodipodi:docname=\"" << myImageName<<"\""<< std::endl;
    myOutputStream << "inkscape:version=\"0.92.2 5c3e80d, 2017-08-06\">" << std::endl;
    myOutputStream << "<metadata " << std::endl;
    myOutputStream << "id=\"metadata9\"> " << std::endl;
    myOutputStream << "<rdf:RDF> " << std::endl;
    myOutputStream << "<cc:Work " << std::endl;
    myOutputStream << "rdf:about=\"\"> " << std::endl;
    myOutputStream << "<dc:format>image/svg+xml</dc:format> " << std::endl;
    myOutputStream << "<dc:type " << std::endl;
    myOutputStream << "rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" /> " << std::endl;
    myOutputStream << "<dc:title></dc:title> " << std::endl;
    myOutputStream << "</cc:Work> " << std::endl;
    myOutputStream << "</rdf:RDF> " << std::endl;
    myOutputStream << "</metadata> " << std::endl;
    myOutputStream << "<defs " << std::endl;
    myOutputStream << "id=\"defs7\" /> " << std::endl;
    myOutputStream << "<sodipodi:namedview" << std::endl;
    myOutputStream << "pagecolor=\"#ffffff\" " << std::endl;
    myOutputStream << "bordercolor=\"#666666\"" << std::endl;
    myOutputStream << "borderopacity=\"1\"" << std::endl;
    myOutputStream << "objecttolerance=\"10\"" << std::endl;
    myOutputStream << "gridtolerance=\"10\" " << std::endl;
    myOutputStream << "guidetolerance=\"10\"" << std::endl;
    myOutputStream << "inkscape:pageopacity=\"0\"" << std::endl;
    myOutputStream << "inkscape:pageshadow=\"2\" " << std::endl;
    myOutputStream << "inkscape:window-width=\"1250\" " << std::endl;
    myOutputStream << "inkscape:window-height=\"735\"" << std::endl;
    myOutputStream << "id=\"namedview5\" " << std::endl;
    myOutputStream << "showgrid=\"false\" " << std::endl;
    myOutputStream << "inkscape:zoom=\"0.2102413\"" << std::endl;
    myOutputStream << "inkscape:cx=\"396.85039\" " << std::endl;
    myOutputStream << "inkscape:cy=\"561.25984\" " << std::endl;
    myOutputStream << "inkscape:window-x=\"0\" " << std::endl;
    myOutputStream << "inkscape:window-y=\"0\"" << std::endl;
    myOutputStream << "inkscape:window-maximized=\"0\"" << std::endl;
    myOutputStream << "inkscape:current-layer=\"svg3\" />" << std::endl;

}

void BasicVectoImageExporter::closeFigure(){
  if(myExportType == SvgExport)
  {
    myOutputStream << "</svg>"<< std::endl;
  }
  
}
