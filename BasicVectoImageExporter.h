#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/geometry/curves/FrechetShortcut.h>
#include <DGtal/shapes/Mesh.h>


#ifndef BASICVECTOIMAGEEXPORTER_H
#define BASICVECTOIMAGEEXPORTER_H


class BasicVectoImageExporter{
  std::string LINE_COLOR;// = " 0.8 0.1 0.1 ";
  std::string POINT_COLOR;// = " 0.1 0.1 0.8 ";
  
  typedef enum {EpsExport, SvgExport, UnknowExport} ExportType;
  
  
public:

  typedef DGtal::Z2i::RealPoint Point2D;
  typedef std::vector<Point2D> Contour2D;
  
  void fillHeader();
  void fillSVGHeader();
  void fillEPSHeader();
  std::string getExportType();
  


  template<typename TContour2D>
  void addPathContent(const TContour2D &contour){
    if( contour.size() == 0 )
    {
      return;
    }
    myOutputStream << contour[0][0] << " " << contour[0][1] << " moveto" << std::endl;
  
    for(unsigned int i = 1; i<contour.size(); i++)
    {
      myOutputStream << contour[i][0] << " " << contour[i][1] << " lineto" << std::endl;
    }
    myOutputStream << contour[0][0] << " " << contour[0][1] << " lineto" << std::endl;
  };

  

  template<typename TContour>
  void addPathContentBezierP0P1P2P3(const TContour &contour)
    {
      if( contour.size() == 0 )
      {
        return;
      }
      myOutputStream << contour[0][0] << " " << contour[0][1] << " moveto" << std::endl;

      for( int i = 1; i<(int)(contour.size())-2; i=i+3)
      {

        myOutputStream << contour[i][0] << " " << contour[i][1] << " ";
        myOutputStream << contour[i+1][0] << " " << contour[i+1][1] << " ";
        myOutputStream << contour[i+2][0] << " " << contour[i+2][1] << " curveto" << std::endl;
      }
    };

 

  
  template<typename TContour>
  void addRegion(const TContour &contour,
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

  template<typename TContour>
  void addRegions(const std::vector<TContour> &contours, const  DGtal::Color &color)
    {
      myOutputStream << "newpath" << std::endl;
      for(auto const &cnt: contours)
      {
        addPathContent(cnt);
      }
      myOutputStream << " closepath " << std::endl;
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
        myOutputStream  << LINE_COLOR <<  "setrgbcolor" << std::endl;
        myOutputStream  <<"0.1 setlinewidth" << std::endl;
        myOutputStream << "stroke" << std::endl;
        myOutputStream  << POINT_COLOR <<  "setrgbcolor" << std::endl;
        for(auto const &cnt: contours)
        {
          addContourPoints(cnt);
        } 
      }


    };

  template<typename TContour>
  void addRegionWithHoles(const TContour &contour,
                          const std::vector<TContour> &listHoles,
                          const  DGtal::Color &color )
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
  
  


  template<typename TPoint2D>
  void drawLine(const TPoint2D &pt1, const TPoint2D &pt2,
                const  DGtal::Color &color, double lineWidth=2.0)
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

  
  void addContour(const std::vector<BasicVectoImageExporter::Point2D> &contour,
                  const  DGtal::Color &color, double lineWidth=1.0)
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


  template<typename TContour>
  void addPathContentBezier(const TContour &contour)
    {
      if( contour.size() <= 1 )
      {
        return;
      }
      // format from dominantPointPolygonalisation_Bezier() (in VectorisationHelper) 
      myOutputStream << contour[2][0] << " " << contour[2][1] << " moveto" << std::endl;  
      for( int i =0; i<(int)(contour.size()); i=i+4)
      {
        myOutputStream << contour[(i+1)%contour.size()][0] << " " << contour[(i+1)%contour.size()][1] << " ";
        myOutputStream << contour[(i+4)%contour.size()][0] << " " << contour[(i+4)%contour.size()][1] << " ";
        myOutputStream << contour[(i+6)%contour.size()][0] << " " << contour[(i+6)%contour.size()][1] << " curveto" << std::endl;
      }

    }

  
  template<typename TContour>
  void addRegionsBezier(const std::vector<TContour> &contours, const  DGtal::Color &color, bool basicOrder=false){
    myOutputStream << "newpath" << std::endl;
    for(auto const &cnt: contours)
    {
      if(basicOrder){
        addPathContentBezierP0P1P2P3(cnt);
      }
      else
      {
        addPathContentBezier(cnt);
      }
    }
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
      myOutputStream << LINE_COLOR <<  "setrgbcolor" << std::endl;
      myOutputStream <<"0.1 setlinewidth" << std::endl;
      myOutputStream << "stroke" << std::endl;
      myOutputStream << POINT_COLOR <<  "setrgbcolor" << std::endl;
      for(auto const &cnt: contours)
      {
        addContourPoints(cnt);
      } 
    }

  }
  



  template<typename TContour>
  void addContourPoints(const TContour &contour,  const  DGtal::Color &color=DGtal::Color::Red, double radius=2.0)
    {
      float r,g,b;
      r = color.red()/255.0;
      g = color.green()/255.0;  
      b = color.blue()/255.0;
      myOutputStream  << r << " " << g << " " << b <<  " setrgbcolor" << std::endl;
      if( contour.size() == 0 )
      {
        return;
      }
  
      for(const auto &p: contour)
      {
        myOutputStream << p[0] << " " << p[1] << " moveto" << std::endl;
        myOutputStream << std::fixed << p[0] << " " << p[1] << " "<< radius<< " 0  360 arc" << std::endl;
        myOutputStream << "fill" << std::endl;
      }
  
  
    };

  

  BasicVectoImageExporter(const std::string &imageName, unsigned int width, unsigned int height,
                          bool displayMesh = false, double  scale=1.0);
  
  ~BasicVectoImageExporter(){myOutputStream.close();};

  

  
  
protected:
  unsigned int myWidth = 200;
  unsigned int myHeight = 200;
  double myScale=1.0 ;

  double myShiftX =  0.0;
  double myShiftY =  0.0;

  std::vector<Contour2D> myPlainContours;
  std::vector<Contour2D> myHoleContours;
  std::string myImageName;
  std::ofstream myOutputStream;
  bool myDisplayMesh;
  
  // Associate for each plain contours a set of index representing the contour holes.
  std::map<unsigned int, std::vector<unsigned int> > mapHoles; 
  std::map<unsigned int, DGtal::Color> colorMap;
  ExportType myExportType=UnknowExport;
  
};


#endif // BASICVECTOIMAGEEXPORTER_H
