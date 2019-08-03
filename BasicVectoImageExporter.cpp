#include "BasicVectoImageExporter.h"

#include <sstream> 







BasicVectoImageExporter::BasicVectoImageExporter(const std::string &imageName,
                                                 unsigned int width, unsigned int height,
                                                 bool displayMesh, double scale): myScale(scale)
{
  myWidth=width;
  myHeight=height;
  myOutputStream.open(imageName, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
  fillEPSHeader();
  myDisplayMesh = displayMesh;
  LINE_COLOR  = " 0.8 0.1 0.1 ";
  POINT_COLOR = " 0.1 0.1 0.8 ";
}





void BasicVectoImageExporter::drawLine(const Point2D &pt1, const Point2D &pt2,  const  DGtal::Color &color, double lineWidth)
{
  float r,g,b;
  r = color.red()/255.0;
  g = color.green()/255.0;  
  b = color.blue()/255.0;
  myOutputStream    << r << " " << g << " " << b <<  " setrgbcolor" << std::endl;
  myOutputStream    << lineWidth << " setlinewidth" << std::endl; 
  myOutputStream    << pt1[0] << " " << pt1[1] << " moveto" << std::endl;
  myOutputStream    << pt2[0] << " " << pt2[1] << " lineto" << std::endl;
  myOutputStream    << "stroke" << std::endl;
      
}


void BasicVectoImageExporter::addContour(const std::vector<BasicVectoImageExporter::Point2D> &contour,
                                         const  DGtal::Color &color, double lineWidth)
{
  myOutputStream << "newpath" << std::endl;
  addPathContent(contour);
  myOutputStream << "closepath" << std::endl;
  float r,g,b;
  r = color.red()/255.0;
  g = color.green()/255.0;  
  b = color.blue()/255.0;
  myOutputStream  << r << " " << g << " " << b <<  " setrgbcolor" << std::endl;
  
  myOutputStream << lineWidth <<"  setlinewidth stroke" << std::endl;
  
  
}




void BasicVectoImageExporter::addRegion(const std::vector<BasicVectoImageExporter::Point2D> &contour,
                                        const  DGtal::Color &color, double linewidth)
{
  myOutputStream << "newpath" << std::endl;
  addPathContent(contour);
  myOutputStream << "closepath" << std::endl;
  if(myDisplayMesh)
    {
      myOutputStream << "gsave" << std::endl;
    }
  float r,g,b;
  r = color.red()/255.0;
  g = color.green()/255.0;  
  b = color.blue()/255.0;
  myOutputStream  << r << " " << g << " " << b <<  " setrgbcolor" << std::endl;
  myOutputStream << "fill" << std::endl;
  if(myDisplayMesh)
    {
      myOutputStream << "grestore" << std::endl;
      myOutputStream  << linewidth << " setlinewidth 0.7 0.2 0.2 setrgbcolor" << std::endl;
      myOutputStream << "stroke" << std::endl;
    }
  
}




void BasicVectoImageExporter::addRegionWithHoles(const std::vector<BasicVectoImageExporter::Point2D> &contour,
                        const std::vector<BasicVectoImageExporter::Contour2D> &listHoles,
                        const DGtal::Color &color )
{
  myOutputStream << "newpath" << std::endl;
  addPathContent(contour);
  for(auto const &hole: listHoles)
    {
      addPathContent(hole);
    }
  myOutputStream << "closepath" << std::endl;
  float r,g,b;
  r = color.red()/255.0;
  g = color.green()/255.0;  
  b = color.blue()/255.0;
  myOutputStream  << r << " " << g << " " << b <<  " setrgbcolor" << std::endl;
  myOutputStream << "fill" << std::endl;
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


