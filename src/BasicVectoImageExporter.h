#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/geometry/curves/FrechetShortcut.h>
#include <DGtal/shapes/Mesh.h>


#ifndef BASICVECTOIMAGEEXPORTER_H
#define BASICVECTOIMAGEEXPORTER_H


class BasicVectoImageExporter {
    std::string LINE_COLOR;// = " 0.8 0.1 0.1 ";
    std::string POINT_COLOR;// = " 0.1 0.1 0.8 ";
    const float emptyCntWidth = 0.1; // used in % (larger for better inskape rendering)
    const float meshCntWidth = 0.1;

    typedef enum {
        EpsExport, SvgExport, UnknowExport
    } ExportType;

public:
    typedef DGtal::Z2i::RealPoint Point2D;
    typedef std::vector <Point2D> Contour2D;

    void closeFigure();

    void fillHeader();

    void fillSVGHeader();

    void fillEPSHeader();

    std::string getHexCode(const DGtal::Color &c);

    std::string getExportType();

    template<typename CoordType>
    CoordType reverseYCoord(CoordType y) {
        return myHeight - y;
    }

    int mod(int value, int m) {
        int res = value % m;
        
        return res < 0 ? res + m : res;
    }

    template<typename TContour2D>
    void addPathContent(const TContour2D &contour) {
        if (contour.size() == 0) {
            return;
        }
        if (myExportType == EpsExport) {
            myOutputStream << contour[0][0] << " " << contour[0][1] << " moveto" << std::endl;
            for (unsigned int i = 1; i < contour.size(); i++) {
                myOutputStream << contour[i][0] << " " << contour[i][1] << " lineto" << std::endl;
            }
            myOutputStream << contour[0][0] << " " << contour[0][1] << " lineto" << std::endl;
        } else if (myExportType == SvgExport) {
            if (contour.size() >= 3) {
                // Array of points for the path
                std::vector<std::pair<double, double>> points;
                // Compute all tangents
                int size = contour.size();
                int scale = 4;
                for (int i = 0; i < size; i++)
                {
                    // Get points
                    double point[2] = {contour[mod(i, size)][0], contour[mod(i, size)][1]};
                    double pointBefore[2] = {contour[mod(i-1, size)][0], contour[mod(i-1, size)][1]};
                    double pointAfter[2] = {contour[mod(i+1, size)][0], contour[mod(i+1, size)][1]};
                    // Compute tangent vector
                    double vector[2] = {(pointAfter[0] - pointBefore[0]), (pointAfter[1] - pointBefore[1])};
                    vector[0] /= scale;
                    vector[1] /= scale;
                    // Compute control points
                    double controlBefore[2] = {(point[0] - vector[0]), (point[1] - vector[1])};
                    double controlAfter[2] = {(point[0] + vector[0]), (point[1] + vector[1])};
                    // Add points to vector
                    points.push_back(std::make_pair(controlBefore[0], reverseYCoord(controlBefore[1])));
                    points.push_back(std::make_pair(point[0], reverseYCoord(point[1])));
                    points.push_back(std::make_pair(controlAfter[0], reverseYCoord(controlAfter[1])));
                }
                

                // Move to first point
                myOutputStream << "M " << points[1].first << " " << points[1].second;
                // Create path
                size = points.size();
                for (int i = 2; i < size; i += 3)
                {
                    // Bezier curve order 3
                    myOutputStream << " C " << points[mod(i, size)].first << " " << points[mod(i, size)].second
                                << ", " << points[mod(i+1, size)].first << " " << points[mod(i+1, size)].second
                                << ", " << points[mod(i+2, size)].first << " " << points[mod(i+2, size)].second;
                }
                // Close path
                myOutputStream << " Z";
            } else {
                myOutputStream << "M ";
                myOutputStream << contour[0][0] << "," << reverseYCoord(contour[0][1]) << " ";
                for (unsigned int i = 1; i < contour.size(); i++) {
                    myOutputStream << contour[i][0] << "," << reverseYCoord(contour[i][1]) << " ";
                }
                myOutputStream << contour[0][0] << "," << reverseYCoord(contour[0][1]) << " z ";
            }
        }
    };


    // Not yet used (@todo later add SVG case)
    template<typename TContour>
    void addPathContentBezierP0P1P2P3(const TContour &contour) {
        if (myExportType == EpsExport) {
            if (contour.size() == 0) {
                return;
            }
            myOutputStream << contour[0][0] << " " << contour[0][1] << " moveto" << std::endl;

            for (int i = 1; i < (int) (contour.size()) - 2; i = i + 3) {

                myOutputStream << contour[i][0] << " " << contour[i][1] << " ";
                myOutputStream << contour[i + 1][0] << " " << contour[i + 1][1] << " ";
                myOutputStream << contour[i + 2][0] << " " << contour[i + 2][1] << " curveto" << std::endl;
            }
        } else if (myExportType == SvgExport) {
            if (contour.size < 1) return;

            myOutputStream << "M " << contour[0][0] << " " << contour[0][1];

            for (int i = 1; i < contour.size(); i += 2) {
                myOutputStream << " Q " << contour[i % contour.size()][0] << " " << contour[i % contour.size()][1]
                               << ", " << contour[(i + 1) % contour.size()][0] << " "
                               << contour[(i + 1) % contour.size()][1];
            }
        }
    };


    // Used (@todo tp add SVG case)
    template<typename TContour>
    void addRegion(const TContour &contour,
                   const DGtal::Color &color, double linewidth) {
        if (myExportType == EpsExport) {
            myOutputStream << "newpath" << std::endl;
            addPathContent(contour);
            myOutputStream << "closepath" << std::endl;
            if (myDisplayMesh) {
                myOutputStream << "gsave" << std::endl;
            }
            float r, g, b;
            r = color.red() / 255.0;
            g = color.green() / 255.0;
            b = color.blue() / 255.0;
            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;
            myOutputStream << "fill" << std::endl;
            if (myDisplayMesh) {
                myOutputStream << "grestore" << std::endl;
                myOutputStream << linewidth << " setlinewidth 0.7 0.2 0.2 setrgbcolor" << std::endl;
                myOutputStream << "stroke" << std::endl;
            } else {
                myOutputStream << "grestore" << std::endl;
                myOutputStream << emptyCntWidth << " setlinewidth " << r << " " << g << " " << b << " setrgbcolor"
                               << std::endl;
                myOutputStream << "stroke" << std::endl;
            }
        } else if (myExportType == SvgExport) {
            if (!myDisplayMesh) {
                myOutputStream << "<path \n style=\"fill:#" << getHexCode(color);
                myOutputStream << "; fill-opacity:1,fill-rule:evenodd;stroke:#" << getHexCode(color) << ";stroke-width:"
                               << emptyCntWidth << "px;";
                myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
                myOutputStream << "d=\"";
            } else {
                myOutputStream << "<path \n style=\"fill:#" << getHexCode(color);
                myOutputStream << "; fill-opacity:1,fill-rule:evenodd;stroke:red;stroke-width:" << meshCntWidth
                               << "px;";
                myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
                myOutputStream << "d=\"";

            }

            addPathContent(contour);
            myOutputStream << "\"\n";
            myOutputStream << "id=\"path" << myCurrentIdPath << "\" \n";
            myOutputStream << "inkscape:connector-curvature=\"0\" ",
                    myOutputStream << "sodipodi:nodetypes=\"cccccccccc\" />\n";
            myCurrentIdPath++;
        }
    }

    // Used (@todo tp add SVG case)
    template<typename TContour>
    void addRegions(const std::vector <TContour> &contours, const DGtal::Color &color) {
        if (myExportType == EpsExport) {
            myOutputStream << emptyCntWidth << " setlinewidth" << std::endl;
            float r, g, b;
            r = color.red() / 255.0;
            g = color.green() / 255.0;
            b = color.blue() / 255.0;
            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;

            myOutputStream << "newpath" << std::endl;
            for (auto const &cnt: contours) {
                addPathContent(cnt);
            }
            myOutputStream << " closepath " << std::endl;
            // myOutputStream << "gsave" << std::endl;


            myOutputStream << "gsave" << std::endl;
            myOutputStream << "fill" << std::endl;
            myOutputStream << "grestore" << std::endl;
            myOutputStream << "stroke" << std::endl;


            if (myDisplayMesh) {
                myOutputStream << "grestore" << std::endl;
                myOutputStream << LINE_COLOR << " setrgbcolor" << std::endl;
                myOutputStream << meshCntWidth << " setlinewidth" << std::endl;
                myOutputStream << "stroke" << std::endl;
                myOutputStream << POINT_COLOR << " setrgbcolor" << std::endl;
                for (auto const &cnt: contours) {
                    addContourPoints(cnt);
                }
            }


        } else if (myExportType == SvgExport) {
            if (!myDisplayMesh) {
                myOutputStream << "<path \n style=\"fill:#" << getHexCode(color);
                myOutputStream << "; fill-opacity:1,fill-rule:evenodd;stroke:#" << getHexCode(color) << ";stroke-width:"
                               << emptyCntWidth << ";";
                myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
                myOutputStream << "d=\"";
            } else {
                myOutputStream << "<path \n style=\"fill:#" << getHexCode(color);
                myOutputStream << "; fill-opacity:1,fill-rule:evenodd;stroke:red;stroke-width:" << meshCntWidth << ";";
                myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
                myOutputStream << "d=\"";

            }
            for (auto const &cnt: contours) {
                addPathContent(cnt);
            }
            myOutputStream << "\"\n";
            myOutputStream << "id=\"path" << myCurrentIdPath << "\" \n";
            myOutputStream << "inkscape:connector-curvature=\"0\" ",
                    myOutputStream << "sodipodi:nodetypes=\"cccccccccc\" />\n";
            myCurrentIdPath++;
            if (myDisplayMesh) {
                for (auto const &cnt: contours) {
                    addContourPoints(cnt);
                }
            }
        }

    };

    template<typename TContour>
    void addRegionWithHoles(const TContour &contour,
                            const std::vector <TContour> &listHoles,
                            const DGtal::Color &color) {
        float r, g, b;
        r = color.red() / 255.0;
        g = color.green() / 255.0;
        b = color.blue() / 255.0;
        if (myExportType == EpsExport) {
            myOutputStream << "newpath" << std::endl;
            addPathContent(contour);
            for (auto const &hole: listHoles) {
                addPathContent(hole);
            }
            myOutputStream << "closepath" << std::endl;

            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;
            myOutputStream << "fill" << std::endl;
        } else if (myExportType == SvgExport) {
            myOutputStream << "<path \n style=\"fill:#" << getHexCode(color);
            myOutputStream << "; fill-opacity:1,fill-rule:evenodd;stroke:#" << getHexCode(color) << ";stroke-width:"
                           << emptyCntWidth << "px;";
            myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
            myOutputStream << "d=\"";
            addPathContent(contour);
            for (auto const &hole: listHoles) {
                addPathContent(hole);
            }
            myOutputStream << "\"\n";
            myOutputStream << "id=\"path" << myCurrentIdPath << "\" \n";
            myOutputStream << "inkscape:connector-curvature=\"0\" ",
                    myOutputStream << "sodipodi:nodetypes=\"cccccccccc\" />\n";
            myCurrentIdPath++;
        }

    }


    // Used (@todo tp add SVG case)
    template<typename TPoint2D>
    void drawLine(const TPoint2D &pt1, const TPoint2D &pt2,
                  const DGtal::Color &color, double lineWidth = 2.0) {
        float r, g, b;
        r = color.red() / 255.0;
        g = color.green() / 255.0;
        b = color.blue() / 255.0;
        if (myExportType == EpsExport) {
            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;
            myOutputStream << lineWidth << " setlinewidth" << std::endl;
            myOutputStream << pt1[0] << " " << pt1[1] << " moveto" << std::endl;
            myOutputStream << pt2[0] << " " << pt2[1] << " lineto" << std::endl;
            myOutputStream << "stroke" << std::endl;
        } else if (myExportType == SvgExport) {
            myOutputStream << "draw line not implemented in SVG" << std::endl;
        }
    }

    template<typename TPoint2D>
    void addContour(const std::vector <TPoint2D> &contour,
                    const DGtal::Color &color, double lineWidth = 1.0) {
        if (myExportType == EpsExport) {
            myOutputStream << "newpath" << std::endl;
            addPathContent(contour);
            myOutputStream << "closepath" << std::endl;
            float r, g, b;
            r = color.red() / 255.0;
            g = color.green() / 255.0;
            b = color.blue() / 255.0;
            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;

            myOutputStream << lineWidth << "  setlinewidth stroke" << std::endl;
        } else if (myExportType == SvgExport) {

            myOutputStream << "<path \n style=\"stroke:#" << getHexCode(color);
            myOutputStream << "; fill:none; stroke-opacity:1;stroke-width:" << lineWidth << ";";
            myOutputStream << "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
            myOutputStream << "d=\"";

            addPathContent(contour);

            myOutputStream << "\"\n";
            myOutputStream << "id=\"path" << myCurrentIdPath << "\" \n";
            myOutputStream << "inkscape:connector-curvature=\"0\" ",
                    myOutputStream << "sodipodi:nodetypes=\"cccccccccc\" />\n";
            myCurrentIdPath++;
        }
    }

    // Not yet used (@todo later add SVG case)
    template<typename TContour>
    void addPathContentBezier(const TContour &contour) {
        if (myExportType == EpsExport) {
            if (contour.size() <= 1) {
                return;
            }
            // format from dominantPointPolygonalisation_Bezier() (in VectorisationHelper)
            myOutputStream << contour[2][0] << " " << contour[2][1] << " moveto" << std::endl;
            for (int i = 0; i < (int) (contour.size()); i = i + 4) {
                myOutputStream << contour[(i + 1) % contour.size()][0] << " " << contour[(i + 1) % contour.size()][1]
                               << " ";
                myOutputStream << contour[(i + 4) % contour.size()][0] << " " << contour[(i + 4) % contour.size()][1]
                               << " ";
                myOutputStream << contour[(i + 6) % contour.size()][0] << " " << contour[(i + 6) % contour.size()][1]
                               << " curveto" << std::endl;
            }
        } else if (myExportType == SvgExport) {
            if (contour.size() < 1) return;

            myOutputStream << "M " << contour[0][0] << " " << contour[0][1];
            for (int i = 1; i < contour.size(); i += 2) {
                myOutputStream << " Q " << contour[i % contour.size()][0] << " " << contour[i % contour.size()][1]
                               << ", " << contour[(i + 1) % contour.size()][0] << " "
                               << contour[(i + 1) % contour.size()][1];
            }
        }
    }

    // Not yet used (@todo later add SVG case)
    template<typename TContour>
    void addRegionsBezier(const std::vector <TContour> &contours, const DGtal::Color &color, bool basicOrder = false) {
        myOutputStream << "newpath" << std::endl;
        for (auto const &cnt: contours) {
            if (basicOrder) {
                addPathContentBezierP0P1P2P3(cnt);
            } else {
                addPathContentBezier(cnt);
            }
        }
        myOutputStream << "closepath" << std::endl;
        if (myDisplayMesh) {
            myOutputStream << "gsave" << std::endl;
        }

        float r, g, b;
        r = color.red() / 255.0;
        g = color.green() / 255.0;
        b = color.blue() / 255.0;
        myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        myOutputStream << "fill" << std::endl;
        if (myDisplayMesh) {
            myOutputStream << "grestore" << std::endl;
            myOutputStream << LINE_COLOR << "setrgbcolor" << std::endl;
            myOutputStream << "0.1 setlinewidth" << std::endl;
            myOutputStream << "stroke" << std::endl;
            myOutputStream << POINT_COLOR << "setrgbcolor" << std::endl;
            for (auto const &cnt: contours) {
                addContourPoints(cnt);
            }
        }

    }


    template<typename TContour>
    void addContourPoints(const TContour &contour, const DGtal::Color &color = DGtal::Color::Red, double radius = 2.0) {
        float r, g, b;
        r = color.red() / 255.0;
        g = color.green() / 255.0;
        b = color.blue() / 255.0;
        if (myExportType == EpsExport) {
            myOutputStream << r << " " << g << " " << b << " setrgbcolor" << std::endl;
            if (contour.size() == 0) {
                return;
            }

            for (const auto &p: contour) {
                myOutputStream << p[0] << " " << p[1] << " moveto" << std::endl;
                myOutputStream << std::fixed << p[0] << " " << p[1] << " " << radius << " 0  360 arc" << std::endl;
                myOutputStream << "fill" << std::endl;
            }
        } else if (myExportType == SvgExport) {
            for (const auto &p: contour) {
                myOutputStream << "<circle cx=\"" << p[0] << "\" cy=\"" << reverseYCoord(p[1]) << "\" r=\"" << radius
                               << "\" fill = \"#" << getHexCode(color) << "\"/>";
            }
        }

    };


    BasicVectoImageExporter(const std::string &imageName, unsigned int width, unsigned int height,
                            bool displayMesh = false, double scale = 1.0);

    ~BasicVectoImageExporter() { myOutputStream.close(); };


protected:
    unsigned int myWidth = 200;
    unsigned int myHeight = 200;
    double myScale = 1.0;

    double myShiftX = 0.0;
    double myShiftY = 0.0;
    int myCurrentIdPath = 1;
    std::vector <Contour2D> myPlainContours;
    std::vector <Contour2D> myHoleContours;
    std::string myImageName;
    std::ofstream myOutputStream;
    bool myDisplayMesh;

    // Associate for each plain contours a set of index representing the contour holes.
    std::map<unsigned int, std::vector<unsigned int> > mapHoles;
    std::map<unsigned int, DGtal::Color> colorMap;
    ExportType myExportType = UnknowExport;

};


#endif // BASICVECTOIMAGEEXPORTER_H
