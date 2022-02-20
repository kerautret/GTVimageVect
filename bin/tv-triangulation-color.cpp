#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/base/Common.h>
#include <DGtal/base/ConstAlias.h>
#include <DGtal/base/Alias.h>
#include <DGtal/kernel/BasicPointPredicates.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/SimpleThresholdForegroundPredicate.h>
#include <DGtal/shapes/TriangulatedSurface.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/math/linalg/SimpleMatrix.h>
#include <DGtal/geometry/helpers/ContourHelper.h>
#include <DGtal/base/Clock.h>

#include "CairoViewer.h"
#include "BasicVectoImageExporter.h"
#include "ImageConnecter.h"
#include "ImageTVRegularization.h"
#include "LeastSquares.h"
#include "BezierCurve.h"
#include "BreadthFirstVisitorWithParent.h"

#include "FourierPoisson.h"
#include "cairo.h"

#include "CLI11.hpp"


static const std::string version = "0.1.3";

///////////////////////////////////////////////////////////////////////////////
double randomUniform()
{
  return (double)rand() / (double)RAND_MAX;
}

namespace DGtal
{

  /// We use a specific representation of multivariate polynomial to
  /// speed up its computation.
  /// p( x, y ) := a x^2 + b xy + c y^2 + d x + e y
  template <typename TPoint>
  struct QuadraticPolynomial
  {
    typedef TPoint Point;
    typedef QuadraticPolynomial<Point> Self;
    typedef typename Point::Coordinate Scalar;

    std::array<Scalar, 6> _c;
    Scalar _x0, _y0;

    QuadraticPolynomial()
        : _c{0, 0, 0, 0, 0, 0},
          _x0(0), _y0(0)
    {
    }
    QuadraticPolynomial(Scalar a, Scalar b, Scalar c,
                        Scalar d, Scalar e, Scalar f,
                        Scalar x0, Scalar y0)
        : _c{a, b, c, d, e, f},
          _x0(x0), _y0(y0)
    {
    }

    ~QuadraticPolynomial() {}

    // Global coordinates
    Scalar operator()(Point p) const
    {
      return eval(p[0] - _x0, p[1] - _y0);
    }
    // Local coordinates
    Scalar eval(Scalar x, Scalar y) const
    {
      return x * (_c[3] * x + _c[4] * y + _c[1]) + y * (_c[5] * y + _c[2]) + _c[0];
    }

    void fit(std::vector<Point> X,
             std::vector<Scalar> V)
    {
      _c = LeastSquares<Point>::linearFit(X, V);
      _x0 = X[0][0];
      _y0 = X[0][1];
    }
  };

  /// We use a specific representation to represent a linear fit that
  /// can be straightened around its mid-value.
  /// p( x, y ) := a x^2 + b xy + c y^2 + d x + e y
  template <typename TPoint>
  struct StraightableLinearFit
  {
    typedef TPoint Point;
    typedef StraightableLinearFit<Point> Self;
    typedef typename Point::Coordinate Scalar;

    std::array<Scalar, 6> _c;
    Scalar _x0, _y0;
    Scalar _mv, _midv, _Mv;

    StraightableLinearFit()
        : _c{0, 0, 0, 0, 0, 0},
          _x0(0), _y0(0)
    {
    }
    StraightableLinearFit(Scalar a, Scalar b, Scalar c,
                          Scalar d, Scalar e, Scalar f,
                          Scalar x0, Scalar y0)
        : _c{a, b, c, d, e, f},
          _x0(x0), _y0(y0)
    {
    }

    ~StraightableLinearFit() {}

    // Global coordinates
    Scalar operator()(Point p) const
    {
      return eval(p[0] - _x0, p[1] - _y0);
    }
    // Local coordinates
    Scalar eval(Scalar x, Scalar y) const
    {
      return x * _c[1] + y * _c[2] + _c[0];
    }
    // Global coordinates
    Scalar operator()(Point p, Scalar amp, Scalar stiff) const
    {
      return eval(p[0] - _x0, p[1] - _y0, amp, stiff);
    }
    // Local coordinates
    Scalar eval(Scalar x, Scalar y,
                Scalar amp, Scalar stiff) const
    {
      Scalar v = x * _c[1] + y * _c[2] + _c[0];
      if (v <= _midv)
        return _mv;
      else if ((v - _mv) <= stiff * (_midv - _mv))
        return _mv;
      // return ( v - _mv )*(1.0-amp) + _mv;
      else if ((v - _mv) <= (1.0 - stiff) * (_midv - _mv))
        return _midv;
      else
        return _Mv;
      // from ( stiff*(_midv - _mv)+_mv, 0.5*( _Mv - _mv) *(1.0-amp) + _mv )
      // to   ( _midv, _midv )
      // return ( v - _mv )*(1.0-amp) + _mv;
    }

    void fit(std::vector<Point> X,
             std::vector<Scalar> V)
    {
      // Compute min and max value.
      int m = 0;
      int M = 0;
      if (V[1] < V[m])
        m = 1;
      if (V[2] < V[m])
        m = 2;
      if (V[1] > V[M])
        M = 1;
      if (V[2] > V[M])
        M = 2;
      auto c = LeastSquares<Point>::linearFit(X, V);
      _c[0] = c[0];
      _c[1] = c[1];
      _c[2] = c[2];
      _x0 = X[0][0];
      _y0 = X[0][1];
      _mv = V[m];
      _Mv = V[M];
      _midv = 0.5 * (_Mv + _mv);
    }
  };

  /// This class represents a triangulation with a Total Variation
  /// (TV) energy computed per triangle. The main idea of the class is
  /// to flip edges in order to decrease the total TV energy.
  struct TVTriangulation
  {

    typedef Z2i::Integer Integer;
    typedef Z2i::Point Point;
    typedef Z2i::RealPoint RealPoint;
    typedef Z2i::RealVector RealVector;
    typedef Z2i::Domain Domain;
    typedef TriangulatedSurface<RealPoint> Triangulation;
    typedef Triangulation::VertexIndex VertexIndex;
    typedef Triangulation::Vertex Vertex;
    typedef Triangulation::Arc Arc;
    typedef Triangulation::Face Face;
    typedef Triangulation::VertexRange VertexRange;
    typedef Triangulation::FaceRange FaceRange;
    typedef double Scalar;
    typedef PointVector<3, Scalar> Value;
    typedef std::pair<DGtal::Color,
                      std::vector<std::vector<std::pair<RealPoint, bool>>>>
        ColorContours;

    // typedef std::array< Value, 2 >     VectorValue;
    struct VectorValue
    {
      Value x;
      Value y;
    };
    typedef std::vector<Scalar> ScalarForm;
    typedef std::vector<Value> ValueForm;
    typedef std::vector<VectorValue> VectorValueForm;
    typedef QuadraticPolynomial<RealPoint> QPolynomial;
    typedef std::array<QPolynomial, 3> QPolynomial3;
    typedef StraightableLinearFit<RealPoint> SLFit;
    typedef std::array<SLFit, 3> SLFit3;
    /// A symmetric definite positive matrix [[a,b],[b,c]]
    struct Metric
    {
      Scalar _a, _b, _c;
      Metric() {}
      Metric(RealPoint p, RealPoint q, RealPoint r)
      {
        RealVector sides[3] = {q - p, r - q, p - r};
        Scalar lengths[3] = {sides[0].norm(),
                             sides[1].norm(),
                             sides[2].norm()};
        int l = 0;
        if (lengths[1] > lengths[l])
          l = 1;
        if (lengths[2] > lengths[l])
          l = 2;
        const RealVector cst = sides[l] / lengths[l]; // [ cos t, sin t ]
        const RealVector ort = {-cst[1], cst[0]};
        const Scalar l1 = lengths[l];
        const Scalar l2 = fabs(sides[(l + 1) % 3].dot(ort));
        _a = cst[0] * cst[0] / (l1) + cst[1] * cst[1] / (l2);
        _b = cst[0] * cst[1] / (l1)-cst[0] * cst[1] / (l2);
        _c = cst[1] * cst[1] / (l1) + cst[0] * cst[0] / (l2);
        // _a = cst[0] * cst[0] / ( l1 * l1 ) + cst[1] * cst[1] / ( l2 * l2 );
        // _b = cst[0] * cst[1] / ( l1 * l1 ) - cst[0] * cst[1] / ( l2 * l2 );
        // _c = cst[1] * cst[1] / ( l1 * l1 ) + cst[0] * cst[0] / ( l2 * l2 );
      }
      Scalar operator()(const RealVector &v, const RealVector &w) const
      {
        return v[0] * (_a * w[0] + _b * w[1]) + v[1] * (_b * w[0] + _c * w[1]);
      }
    };

    /// Contains information per face for 2nd order reconstruction.
    struct FaceInformation
    {
      Value v;          ///< the value at the barycenter of the face.
      VectorValue grad; ///< the constant gradient on the face.
      int order;        ///< Order according to TV discontinuity (0:highest)
      Scalar rel_tv;    ///< relative tv energy wrt average
      Scalar cumul_tv;  ///< cumulative percent of total TV energy
      Scalar disc;      ///< between 0 and 1, means that you consider the face as a discontinuity, more means that you do not consider this face as discontinuous.
      Scalar cont;      ///< below 0 is discontinuous, otherwise between 0 (almost discontinuous) and 1 (continuous).
    };

    /// The domain triangulation
    Triangulation T;
    /// The image values at each vertex
    ValueForm _I;
    /// Initial vertices
    Integer _nbV;
    /// Tells if it is a color image (just for optimization).
    bool _color;
    /// Power for gradient computation
    Scalar _power;
    /// List of arcs whose energy may have been changed by a surrounding flip.
    std::vector<Arc> _Queue;
    /// List of arcs that may be flipped with the same energy.
    std::vector<Arc> _Q_equal;

    /// The TV-regularized values
    ValueForm _u;
    /// The TV-regularized vectors
    VectorValueForm _p;
    /// The vector storing the tv energy of each triangle.
    ScalarForm _tv_per_triangle;
    /// The total variation energy of T.
    Scalar _tv_energy;
    /// true if arc is flippable.
    std::vector<bool> _flippable;

    /// Contains the polynomials per vertex that fits the values.
    std::vector<QPolynomial3> _u_approx;
    /// Contains the polynomials per faces that fits the values.
    std::vector<SLFit3> _u_approx_f;
    /// Width of the image
    Integer _width;
    /// Domain of the input image
    Point _lo;
    Point _up;
    /// Contains the metric of each triangle.
    std::vector<Metric> _metrics;
    /// Contains information per face for 2nd order reconstruction.
    std::vector<FaceInformation> _finfos;
    /// Tells if arcs needs to be updated in regularization
    std::vector<bool> _arc_updatable;
    /// Tells if faces needs to be updated in regularization
    std::vector<bool> _face_updatable;
    /// Contains arc similarities.
    std::vector<Scalar> _arc_dissimilarities;

    // Data needed for vectorization. Each contour is a succession of
    // arc, where head points outside.

    /// Contains the relative positions of the contour at each arc.
    /// A value between 0 and 1: 0 is tail and 1 is head of arc
    std::vector<Scalar> _t;
    /// Contains the coordinates of the barycenter of each face.
    std::vector<RealPoint> _b;
    /// Contains the area for each vertex.
    std::vector<Scalar> _A;

    /// @return the polynomial at vertex \a v
    const QPolynomial3 &uApprox(const Vertex v) const
    {
      ASSERT(v < _u_approx.size());
      return _u_approx[v];
    }

    /// @return the face information associated to face \a f.
    const FaceInformation &finfo(const Face f) const
    {
      ASSERT(f < _finfos.size());
      return _finfos[f];
    }
    /// @return the contour point at arc \a a.
    Scalar arcDissimilarity(const Arc a) const
    {
      return _arc_dissimilarities[a];
    }

    /// @return the regularized value at vertex v.
    const Value &u(const VertexIndex v) const
    {
      return _u[v];
    }

    /// @return the regularized barycenter at face f.
    const RealPoint &barycenter(const Face f) const
    {
      return _b[f];
    }

    /// @return the contour point at arc \a a.
    RealPoint contourPoint(const Arc a) const
    {
      RealPoint B = T.position(T.head(a));
      RealPoint A = T.position(T.tail(a));
      return (1.0 - _t[a]) * A + _t[a] * B;
    }
    /// @return the contour point at arc \a a.
    RealPoint edgeContourPoint(const Arc a) const
    {
      return 0.5 * (contourPoint(a) + contourPoint(T.opposite(a)));
    }

    /// invalidate a vertex by using specific value (to process image border).
    void invalidate(const VertexIndex v)
    {
      _u[v][0] = std::numeric_limits<double>::min();
      _u[v][1] = std::numeric_limits<double>::min();
      _u[v][2] = std::numeric_limits<double>::min();
    }
    bool isinvalid(const VertexIndex v)
    {
      return _u[v][0] == std::numeric_limits<double>::min() &&
             _u[v][1] == std::numeric_limits<double>::min() &&
             _u[v][2] == std::numeric_limits<double>::min();
    }

    /// The norm used for the scalars induced by vector-value space (RGB)
    std::function<Scalar(const Value &v)> _normX;
    /// The norm used for the 2d vectors induced by vector-value space (RGB)
    std::function<Scalar(const VectorValue &v)> _normY;

    static Scalar square(Scalar x) { return x * x; }

    // ------------- Discrete operators ---------------------

    /// @return the vector of the norms of each vector in p (one per triangle).
    ScalarForm norm(const VectorValueForm &p) const
    {
      ScalarForm S(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); f++)
        S[f] = _normY(p[f]);
      return S;
    }
    // Definition of a global gradient operator that assigns vectors to triangles
    VectorValueForm grad(const ValueForm &u) const
    { // it suffices to traverse all (valid) triangles.
      VectorValueForm G(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        G[f] = grad(f, u);
      }
      return G;
    }

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad(Face f, const ValueForm &u) const
    {
      VertexRange V = T.verticesAroundFace(f);
      return grad(V[0], V[1], V[2], u);
    }

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad(VertexIndex i, VertexIndex j, VertexIndex k,
                     const ValueForm &u) const
    {
      // [ yj-yk yk-yi yi-yk ] * [ ui ]
      // [ xk-xj xi-xk xj-xi ]   [ uj ]
      //                         [ uk ]
      const RealPoint &pi = T.position(i);
      const RealPoint &pj = T.position(j);
      const RealPoint &pk = T.position(k);
      const Value &ui = u[i];
      const Value &uj = u[j];
      const Value &uk = u[k];
      VectorValue G;
      const int d = _color ? 3 : 1;
      for (int m = 0; m < d; ++m)
      {
        G.x[m] = ui[m] * (pj[1] - pk[1]) + uj[m] * (pk[1] - pi[1]) + uk[m] * (pi[1] - pj[1]);
        G.y[m] = ui[m] * (pk[0] - pj[0]) + uj[m] * (pi[0] - pk[0]) + uk[m] * (pj[0] - pi[0]);
      }
      G.x *= 0.5;
      G.y *= 0.5;
      return G;
    }

    // Definition of a (global) divergence operator that assigns
    // scalars to vertices from a vector field.
    ValueForm div(const VectorValueForm &G) const
    {
      const int d = _color ? 3 : 1;
      ValueForm S(T.nbVertices());
      for (VertexIndex v = 0; v < T.nbVertices(); ++v)
        S[v] = Value{0, 0, 0};
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        auto V = T.verticesAroundFace(f);
        const RealPoint &pi = T.position(V[0]);
        const RealPoint &pj = T.position(V[1]);
        const RealPoint &pk = T.position(V[2]);
        const VectorValue &Gf = G[f];
        for (int m = 0; m < d; ++m)
        {
          S[V[0]][m] -= (pj[1] - pk[1]) * Gf.x[m] + (pk[0] - pj[0]) * Gf.y[m];
          S[V[1]][m] -= (pk[1] - pi[1]) * Gf.x[m] + (pi[0] - pk[0]) * Gf.y[m];
          S[V[2]][m] -= (pi[1] - pj[1]) * Gf.x[m] + (pj[0] - pi[0]) * Gf.y[m];
        }
      }
      for (VertexIndex v = 0; v < T.nbVertices(); ++v)
        S[v] *= 0.5;
      return S;
    }

    /// @return the scalar form lambda.u
    ValueForm multiplication(Scalar lambda, const ValueForm &u) const
    {
      ValueForm S(T.nbVertices());
      for (VertexIndex v = 0; v < T.nbVertices(); ++v)
        S[v] = lambda * u[v];
      return S;
    }

    /// @return the scalar form u - v
    ValueForm subtraction(const ValueForm &u, const ValueForm &v) const
    {
      ValueForm S(T.nbVertices());
      for (VertexIndex i = 0; i < T.nbVertices(); ++i)
        S[i] = u[i] - v[i];
      return S;
    }
    /// u -= v
    void subtract(ValueForm &u, const ValueForm &v) const
    {
      for (VertexIndex i = 0; i < T.nbVertices(); ++i)
        u[i] -= v[i];
    }

    /// @return the scalar form a.u + b.v
    ValueForm combination(const Scalar a, const ValueForm &u,
                          const Scalar b, const ValueForm &v) const
    {
      ValueForm S(T.nbVertices());
      // trace.info() << "[combination] variation is ["
      // 		   << ( b * *( std::min_element( v.begin(), v.end() ) ) )
      // 		   << " " << ( b * *( std::max_element( v.begin(), v.end() ) ) )
      // 		   << std::endl;
      for (VertexIndex i = 0; i < T.nbVertices(); ++i)
        S[i] = a * u[i] + b * v[i];
      return S;
    }

    /// The function that evaluates the energy at each triangle.
    /// It is now just the norm of the gradient.
    Scalar computeEnergyTV(VertexIndex v1, VertexIndex v2, VertexIndex v3) const
    {
      return _normY(grad(v1, v2, v3, _u));
    }

    /// @return the tv energy stored at this face.
    Scalar computeEnergyTV(const Face f)
    {
      VertexRange V = T.verticesAroundFace(f);
      return (_tv_per_triangle[f] = computeEnergyTV(V[0], V[1], V[2]));
    }

    /// @return the tv energy stored at this face.
    Scalar energyTV(const Face f) const
    {
      return _tv_per_triangle[f];
    }

    /// @return the tv energy stored at this face.

    Scalar &energyTV(const Face f)
    {
      return _tv_per_triangle[f];
    }

    /// Compute (and store in _tv_per_triangle) the TV-energy per triangle.
    Scalar computeEnergyTV()
    {
      Scalar E = 0;
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        E += computeEnergyTV(f);
      }
      _tv_energy = E;
      // trace.info() << "TV(u) = " << E << std::endl;
      return E;
    }

    /// Gets the current TV energy of the triangulation.
    Scalar getEnergyTV()
    {
      return _tv_energy;
    }

    /// @return the aspect ratio of a face (the greater, the most elongated it is.
    Scalar aspectRatio(const Face f) const
    {
      VertexRange P = T.verticesAroundFace(f);
      const RealPoint &a = T.position(P[0]);
      const RealPoint &b = T.position(P[1]);
      const RealPoint &c = T.position(P[2]);
      RealVector ab = b - a;
      RealVector bc = c - b;
      RealVector ca = a - c;
      Scalar dab = ab.norm();
      Scalar dbc = bc.norm();
      Scalar dca = ca.norm();
      RealVector uab = ab / dab;
      RealVector ubc = bc / dbc;
      RealVector uca = ca / dca;
      Scalar ha = (ab - ab.dot(ubc) * ubc).norm();
      Scalar hb = (bc - bc.dot(uca) * uca).norm();
      Scalar hc = (ca - ca.dot(uab) * uab).norm();
      return std::max(dab / hc, std::max(dbc / ha, dca / hb));
    }

    /// @return the diameter of a face (the greater, the most elongated it is.
    Scalar diameter(const Face f) const
    {
      VertexRange P = T.verticesAroundFace(f);
      const RealPoint &a = T.position(P[0]);
      const RealPoint &b = T.position(P[1]);
      const RealPoint &c = T.position(P[2]);
      RealVector ab = b - a;
      RealVector bc = c - b;
      RealVector ca = a - c;
      Scalar dab = ab.norm();
      Scalar dbc = bc.norm();
      Scalar dca = ca.norm();
      return std::max(dab, std::max(dbc, dca));
    }

    Value meanValue() const
    {
      Value m;
      for (VertexIndex v = 0; v < _u.size(); ++v)
        m += _u[v];
      m /= _u.size();
      return m;
    }

    Domain getFaceDomain(Face f) const
    {
      VertexRange P = T.verticesAroundFace(f);
      const RealPoint &a = T.position(P[0]);
      const RealPoint &b = T.position(P[1]);
      const RealPoint &c = T.position(P[2]);
      RealPoint lo = a.inf(b).inf(c);
      RealPoint hi = a.sup(b).sup(c);
      return Domain(Point(floor(lo[0]), floor(lo[1])),
                    Point(ceil(hi[0]), ceil(hi[1])));
    }

    template <typename ScalarImage>
    void fillZoomedImage(ScalarImage &output,
                         const image::Zoom &Z, Dimension d) const
    {
      output = ScalarImage(Z.zoomedDomain());
      for (auto &&zp : Z.zoomedDomain())
      {
        const Point p = Z.unzoom(zp);
        const VertexIndex v = Z.index(p);
        output.setValue(zp, _u[v][d]);
      }
    }

    template <typename ScalarImage>
    void fillZoomedGradient(ScalarImage &output_x, ScalarImage &output_y,
                            const image::Zoom &Z, Dimension d) const
    {
      output_x = ScalarImage(Z.zoomedDomain());
      output_y = ScalarImage(Z.zoomedDomain());
      ScalarImage nb = ScalarImage(Z.zoomedDomain());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        fillTriangleZoomedGradientSharp(nb, output_x, output_y, f, Z, d);
      }
      // Diminue la dynamique de l'image.
      // auto it_gx = output_x.begin();
      // auto it_gy = output_y.begin();
      // for ( auto it_nb = nb.cbegin(), itE_nb = nb.cend(); it_nb != itE_nb; ++it_nb )
      // 	{
      // 	  Scalar n = *it_nb;
      // 	  if ( n > 0.0 ) {
      // 	    *it_gx++ /= n;
      // 	    *it_gy++ /= n;
      // 	  }
      // 	}
    }
    template <typename ScalarImage>
    void fillTriangleZoomedGradient(ScalarImage &nb,
                                    ScalarImage &output_x, ScalarImage &output_y,
                                    Face f, const image::Zoom &Z, Dimension d) const
    {
      VectorValue G = grad(f, _u);
      Domain D = getFaceDomain(f);
      Domain ZD = Domain(Z.zoom(D.lowerBound()), Z.zoom(D.upperBound()));
      double weight = 0.25 * (double)Z._zoom; // * Z._zoom;
      for (auto &&p : ZD)
      {
        RealPoint pp = Z.project(p);
        if (isApproximatelyInTriangle(f, pp))
        {
          output_x.setValue(p, -G.x[d] / weight); //+ output_x( p ) );
          output_y.setValue(p, -G.y[d] / weight); //+ output_y( p ) );
          nb.setValue(p, 1.0 + nb(p));
        }
      }
    }

    template <typename ScalarImage>
    void
    fillTriangleZoomedGradientSharp(ScalarImage &nb,
                                    ScalarImage &output_x, ScalarImage &output_y,
                                    Face f, const image::Zoom &Z, Dimension d) const
    { // JACO
      VectorValue G = grad(f, _u);
      RealVector gn(-G.x[d], -G.y[d]); // normal = grad
      bool flat = gn.norm() == 0.0;
      if (flat)
      {
        //fillTriangleZoomedGradient( nb, output_x, output_y, f, Z, d );
        return;
      }
      Domain D = getFaceDomain(f);
      Domain ZD = Domain(Z.zoom(D.lowerBound()), Z.zoom(D.upperBound()));
      VertexRange P = T.verticesAroundFace(f);
      RealPoint V[3] = {T.position(P[0]),
                        T.position(P[1]),
                        T.position(P[2])};
      const double Gnorm = gn.norm();
      gn /= Gnorm;
      double nmax = V[0].dot(gn);
      double nmin = V[0].dot(gn);
      for (unsigned int i = 1; i < 3; ++i)
      {
        nmax = std::max(V[i].dot(gn), nmax);
        nmin = std::min(V[i].dot(gn), nmin);
      }
      double distance = nmax - nmin;
      double middle = 0.5 * (nmin + nmax);
      double lambda = Gnorm > 20.0 ? 1.0 : Gnorm / 20.0;
      double sigma = 1.0 * (1.0 - lambda) + lambda / (2.0 * (double)Z._zoom);
      double c1 = -1.0 / (2.0 * sigma * sigma);
      double c2 = 1.0 / (sigma * sqrt(2.0 * M_PI));
      const double ph = 0.5 / (double)Z._zoom;
      std::vector<Point> points;
      std::vector<double> weights;
      for (auto &&p : ZD)
      {
        RealPoint pp = Z.project(p);
        if (isApproximatelyInTriangle(f, pp))
        {
          double xx = gn.dot(pp) - middle;
          double w = c2 * exp(c1 * xx * xx);
          // Une vraie intégration ne change quasi rien.
          // double x[ 4 ] = { gn.dot( pp + RealVector( ph, ph ) ) - middle,
          // 		      gn.dot( pp + RealVector( ph,-ph ) ) - middle,
          // 		      gn.dot( pp + RealVector(-ph,-ph ) ) - middle,
          // 		      gn.dot( pp + RealVector(-ph, ph ) ) - middle };
          // double minx = std::min( std::min( x[0], x[1] ), std::min( x[2], x[3] ) );
          // double maxx = std::max( std::max( x[0], x[1] ), std::max( x[2], x[3] ) );
          // const int n = 8; // change pas grand chose
          // double    h = (maxx-minx) / (double) n;
          // double    w = 0.0;
          // // integration by middle point
          // for ( double xx = minx + 0.5*h; xx < maxx; xx += h )
          //   w += c2 * exp( c1 * xx * xx );
          points.push_back(p);
          weights.push_back(w);
        }
      }
      double total_weight = 0;
      for (auto &&w : weights)
        total_weight += w;
      const double target = 1.0 * (double)Z._zoom; // area of a triangle
      total_weight /= target;
      for (unsigned int i = 0; i < points.size(); ++i)
      {
        const Point &p = points[i];
        const double w = weights[i] / total_weight;
        output_x.setValue(p, output_x(p) - G.x[d] * w);
        output_y.setValue(p, output_y(p) - G.y[d] * w);
        nb.setValue(p, 1.0 + nb(p));
      }
    }

    // template <typename Image>
    // void outputPoissonImage( Image& image,
    //     		     const image::Zoom& Z )
    // {
    //   typedef ImageContainerBySTLVector< Domain, Scalar > ScalarImage;
    //   const Dimension max_d = _color ? 3 : 1;
    //   image = Image( Z.zoomedDomain() );
    //   std::vector< ScalarImage*> Id( max_d );
    //   for ( Dimension d = 0; d < max_d; ++d )
    //     {
    //       Id[ d ] = new ScalarImage( Z.zoomedDomain() );
    //       ScalarImage Gxd( Z.zoomedDomain() );
    //       ScalarImage Gyd( Z.zoomedDomain() );
    //       fillZoomedImage( *( Id[ d ] ), Z, d );
    //       fillZoomedGradient( Gxd, Gyd, Z, d );
    //       if ( d == 0 ) {
    //         image::functions::exportGradPGM( "gradx.pgm", Gxd );
    //         image::functions::exportGradPGM( "grady.pgm", Gyd );
    //       }
    //       image::functions::poisson( *( Id[ d ] ), Gxd, Gyd );
    //     }
    //   if ( _color )
    //     for ( auto&& p : image.domain() ) {
    //       Color col( (unsigned char) round( (*Id[ 0 ])( p ) ),
    //     	     (unsigned char) round( (*Id[ 1 ])( p ) ),
    //     	     (unsigned char) round( (*Id[ 2 ])( p ) ) );
    //       image.setValue( p, col );
    //     }
    //   else
    //     for ( auto&& p : image.domain() ) {
    //       image.setValue( p, (unsigned char) round( (*Id[ 0 ])( p ) ) );
    //     }
    //   for ( Dimension d = 0; d < max_d; ++d )
    //     delete Id[ d ];
    // }

    // -------------- Construction services -------------------------

    // Constructor from domain, triangulation, and value form.
    TVTriangulation(const Domain &aDomain,
                    const ValueForm &anI,
                    const Triangulation &aT,
                    bool color,
                    Scalar p = 0.5, Scalar sim = 0.0)
    {
      _color = color;
      _power = p;
      _lo = aDomain.lowerBound();
      _up = aDomain.upperBound();
      _width = _up[0] - _lo[0] + 1;
      _I = anI;
      T = aT;
      // Defining norms.
      if (color)
      {
        _normX = [p](const Value &v) -> Scalar
        {
          return pow(square(v[0]) + square(v[1]) + square(v[2]), p);
        };
        // _normY = [p] ( const VectorValue& v ) -> Scalar
        //   {
        //     return pow( square( v.x[ 0 ] ) + square( v.y[ 0 ] ), p )
        //     + pow( square( v.x[ 1 ] ) + square( v.y[ 1 ] ), p )
        //     + pow( square( v.x[ 2 ] ) + square( v.y[ 2 ] ), p );
        //   };
        // Standard ColorTV is
        _normY = [p](const VectorValue &v) -> Scalar
        {
          return pow(square(v.x[0]) + square(v.y[0]) + square(v.x[1]) + square(v.y[1]) + square(v.x[2]) + square(v.y[2]), p);
        };
      }
      else
      {
        _normX = [p](const Value &v) -> Scalar
        {
          return pow(square(v[0]), p);
        };
        _normY = [p](const VectorValue &v) -> Scalar
        {
          return pow(square(v.x[0]) + square(v.y[0]), p);
        };
      }
      _nbV = T.nbVertices();
      // Building forms.
      _u = _I;                // u = image at initialization
      _p.resize(T.nbFaces()); // p = 0     at initialization
      // TV-energy is computed and stored per face to speed-up computations.
      _tv_per_triangle.resize(T.nbFaces());
      computeEnergyTV();

      // Fix some arcs;
      _flippable.resize(T.nbArcs());
      int nbFlippable = 0;
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        _flippable[a] = true;
        nbFlippable += _flippable[a] ? 1 : 0;
      }
      trace.info() << "Nb arcs flippable = " << nbFlippable
                   << "/" << T.nbArcs() << std::endl;
    }

    // Constructor from color image.
    template <typename Image>
    TVTriangulation(const Image &I, bool color,
                    Scalar p = 0.5, Scalar sim = 0.0,
                    int connectivity_strategy = 1, bool debug = false)
    {
      _color = color;
      _power = p;
      _lo = I.domain().lowerBound();
      _up = I.domain().upperBound();
      _width = _up[0] - _lo[0] + 1;
      // Defining norms.
      if (color)
      {
        _normX = [p](const Value &v) -> Scalar
        {
          return pow(square(v[0]) + square(v[1]) + square(v[2]), p);
        };
        // _normY = [p] ( const VectorValue& v ) -> Scalar
        //   {
        //     return pow( square( v.x[ 0 ] ) + square( v.y[ 0 ] ), p )
        //     + pow( square( v.x[ 1 ] ) + square( v.y[ 1 ] ), p )
        //     + pow( square( v.x[ 2 ] ) + square( v.y[ 2 ] ), p );
        //   };
        // Standard ColorTV is
        _normY = [p](const VectorValue &v) -> Scalar
        {
          return pow(square(v.x[0]) + square(v.y[0]) + square(v.x[1]) + square(v.y[1]) + square(v.x[2]) + square(v.y[2]), p);
        };
      }
      else
      {
        _normX = [p](const Value &v) -> Scalar
        {
          return pow(square(v[0]), p);
        };
        _normY = [p](const VectorValue &v) -> Scalar
        {
          return pow(square(v.x[0]) + square(v.y[0]), p);
        };
      }
      // Creates image form _I
      image::Utils im_utils;
      _I = im_utils.getValueForm(I, color);

      // Building connections.
      typedef ImageContainerBySTLVector<Domain, Value> ValueImage;
      typedef ImageConnecter<ValueImage> Connecter;
      ValueImage tmpI(I.domain());
      auto it = tmpI.begin();
      for (auto v : _I)
        *it++ = v;

      Connecter connecter(debug);
      typename Connecter::Comparator comp = [this](Value v1, Value v2)
      { return _normX(v1 - v2); };
      trace.info() << "Compute connections ... ";
      auto ccn_strategy = connectivity_strategy == 0
                              ? Connecter::SizeOfConnectedComponents
                              : Connecter::OrderOfConnectedComponents;
      connecter.init(ccn_strategy, tmpI, comp, sim);
      trace.info() << "ended." << std::endl;

      // Building triangulation
      const Point taille = I.extent();
      // Creates vertices
      for (auto p : I.domain())
        T.addVertex(p);
      // Creates triangles
      for (Integer y = 0; y < taille[1] - 1; ++y)
      {
        for (Integer x = 0; x < taille[0] - 1; ++x)
        {
          const VertexIndex v00 = y * taille[0] + x;
          const VertexIndex v10 = v00 + 1;
          const VertexIndex v01 = v00 + taille[0];
          const VertexIndex v11 = v01 + 1;
          auto how = connecter.howConnected(Point(x, y));
          bool diag00_11 = (how.diagonal == Connecter::Diagonal00_11);
          if (diag00_11)
          {
            T.addTriangle(v00, v01, v11);
            T.addTriangle(v00, v11, v10);
          }
          else
          {
            T.addTriangle(v00, v01, v10);
            T.addTriangle(v10, v01, v11);
          }
        }
      }
      bool ok = T.build();
      trace.info() << "Build triangulation: "
                   << (ok ? "OK" : "ERROR") << std::endl;
      _nbV = T.nbVertices();
      // Building forms.
      _u = _I;                // u = image at initialization
      _p.resize(T.nbFaces()); // p = 0     at initialization
      // TV-energy is computed and stored per face to speed-up computations.
      _tv_per_triangle.resize(T.nbFaces());
      computeEnergyTV();

      // Fix some arcs;
      _flippable.resize(T.nbArcs());
      int nbFlippable = 0;
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        RealPoint p = T.position(T.head(a));
        RealPoint q = T.position(T.tail(a));
        RealPoint l = p.inf(q);
        RealPoint u = p.sup(q);
        auto how = connecter.howConnected(Point((int)round(l[0]), (int)round(l[1])));
        if ((l - u).dot(l - u) == 2)
          _flippable[a] = (how.diagonal == Connecter::Default);
        else if (u[0] != l[0])
          _flippable[a] = !how.horizontal;
        else
          _flippable[a] = !how.vertical;
        nbFlippable += _flippable[a] ? 1 : 0;
      }
      trace.info() << "Nb arcs flippable = " << nbFlippable
                   << "/" << T.nbArcs() << std::endl;
    }

    // -------------------------- Laplacian services --------------------
  public:
    /// Computes and outputs the norm of the image Laplacian.
    template <typename Image>
    void outputLaplacian(Image &L, const Scalar factor) const
    {
      auto lap_bel = div(grad(_u));
      Vertex vtx = 0;
      for (unsigned int &val : L)
      {
        val = (unsigned int)std::min(255, std::max(0, (int)round(factor * _normX(lap_bel[vtx]))));
        vtx++;
      }
    }

    template <typename Image>
    bool outputU(Image &J) const
    {
      VertexIndex v = 0;
      for (unsigned int &val : J)
      {
        val = _color
                  ? (((int)_u[v][0]) << 16) + (((int)_u[v][1]) << 8) + ((int)_u[v][2])
                  : (((int)_u[v][0]) << 16) + (((int)_u[v][0]) << 8) + ((int)_u[v][0]);
        v += 1;
      }
      return v == T.nbVertices();
    }

    // @return the determinant of \a pq and \a qr.
    static Scalar det(const RealPoint &pq, const RealPoint &qr)
    {
      return pq[0] * qr[1] - pq[1] * qr[0];
    }

    static Scalar doesTurnLeft(const RealPoint &p, const RealPoint &q, const RealPoint &r)
    {
      const RealPoint pq = q - p;
      const RealPoint qr = r - q;
      return det(pq, qr);
    }
    static Scalar doesTurnLeft(const RealPoint &pq, const RealPoint &qr)
    {
      return det(pq, qr);
    }

    // Check strict convexity of quadrilateron.
    bool isConvex(const VertexRange &V) const
    {
      RealPoint P[] = {T.position(V[1]) - T.position(V[0]),
                       T.position(V[2]) - T.position(V[1]),
                       T.position(V[3]) - T.position(V[2]),
                       T.position(V[0]) - T.position(V[3])};
      bool cvx = (doesTurnLeft(P[0], P[1]) < 0) && (doesTurnLeft(P[1], P[2]) < 0) && (doesTurnLeft(P[2], P[3]) < 0) && (doesTurnLeft(P[3], P[0]) < 0);
      return cvx;
    }
    /**
       NB: process only arcs (s,t) where ( t > s ).
       
       @return 1 is energy is lowered (a has been flipped), 0 if arc
       is flippable but it does not change the energy, negative
       otherwise (-1: boundary, -2 s > t, -3 non convex, -4 increases
       the energy.
    */
    int updateArc(const Arc a)
    {
      // Checks that edge can be flipped.
      if (!_flippable[a])
        return -4;
      VertexRange P = T.verticesOfFacesAroundArc(a);
      if (P.size() != 4)
        return -1;
      if (P[0] < P[2])
        return -2;
      if (!isConvex(P))
        return -3;
      // Computes energies
      const Face f012 = T.faceAroundArc(a);
      const Face f023 = T.faceAroundArc(T.opposite(a));
      const Scalar E012 = energyTV(f012); //P[ 0 ], P[ 1 ], P[ 2 ] );
      const Scalar E023 = energyTV(f023); //P[ 0 ], P[ 2 ], P[ 3 ] );
      const Scalar E013 = computeEnergyTV(P[0], P[1], P[3]);
      const Scalar E123 = computeEnergyTV(P[1], P[2], P[3]);
      const Scalar Ecurr = E012 + E023;
      const Scalar Eflip = E013 + E123;
      if (Eflip < Ecurr)
      {
        // Save arcs that may be affected.
        queueSurroundingArcs(a);
        T.flip(a);
        _tv_per_triangle[f012] = E123; // f012 -> f123
        _tv_per_triangle[f023] = E013; // f023 -> f013
        _tv_energy += Eflip - Ecurr;
        return 1;
      }
      else if (Eflip == Ecurr)
      {
        return (Eflip > 0.0) ? 0 : -6;
      }
      else
        return -7;
    }

    void queueSurroundingArcs(const Arc a)
    {
      Arc around[4];
      around[0] = T.next(a);
      around[1] = T.next(around[0]);
      around[2] = T.next(T.opposite(a));
      around[3] = T.next(around[2]);
      for (int i = 0; i < 4; ++i)
      {
        _Queue.push_back(around[i]);
        _Queue.push_back(T.opposite(around[i]));
      }
    }

    /// Quantify the regularized image _u.
    void quantify(const int level)
    {
      const Scalar factor = 255.0 / (level - 1);
      for (VertexIndex i = 0; i < T.nbVertices(); ++i)
        for (int m = 0; m < 3; ++m)
        {
          _u[i][m] = round((_u[i][m]) / factor) * factor;
          _u[i][m] = std::min(255.0, std::max(0.0, _u[i][m]));
        }
      computeEnergyTV();
    }

    /// Does one pass of TV regularization (u, p and I must have the
    /// meaning of the previous iteration).
    Scalar tvPass(Scalar lambda, Scalar dt, Scalar tol, int N = 10)
    {
      const Scalar E = getEnergyTV();
      const Scalar perimeter = E / (255.0 * (_color ? sqrt(3.0) : 1.0));
      trace.info() << "TV( u ) = " << getEnergyTV()
                   << " Perimeter = " << perimeter << std::endl;
      //trace.info() << "lambda.f" << std::endl;
      ScalarForm diam(T.nbFaces());
      ScalarForm diam_v(T.nbVertices());
      for (Vertex vtx = 0; vtx < T.nbVertices(); vtx++)
      {
        auto F = T.facesAroundVertex(vtx);
        Scalar d = 0.0;
        for (auto f : F)
          d = std::max(d, diameter(f));
        diam_v[vtx] = d;
      }
      for (Face f = 0; f < T.nbFaces(); f++)
      {
        auto V = T.verticesAroundFace(f);
        diam[f] = std::max(diam_v[V[0]],
                           std::max(diam_v[V[1]], diam_v[V[2]]));
      }
      _p = VectorValueForm(_p.size());
      // grad( combination( lambda, _u, -lambda, _I ) );
      VectorValueForm p(_p.size()); // this is p^{n+1}
      ValueForm lf = multiplication(lambda, _I);
      Scalar diff_p = 0.0;
      int n = 0; // iteration number
      Scalar last_E = -1.0;
      do
      {
        // trace.info() << "div( p ) - lambda.f" << std::endl;
        ValueForm dp_lf = subtraction(div(_p), lf);
        // trace.info() << "G := grad( div( p ) - lambda.f)" << std::endl;
        VectorValueForm gdp_lf = grad(dp_lf);
        // trace.info() << "N := | G |" << std::endl;
        ScalarForm ngdp_lf = norm(gdp_lf);
        // trace.info() << "p^n+1 := ( p + dt * G ) / ( 1 + dt | G | )" << std::endl;
        diff_p = 0.0;
        for (Face f = 0; f < T.nbFaces(); f++)
        {
          const Scalar ldt = 2.0 * dt / (diam[f] * diam[f]);
          const Scalar alpha = 1.0 / (1.0 + ldt * ngdp_lf[f]);
          if (alpha <= 0.0)
            trace.warning() << "Face " << f << " alpha=" << alpha << std::endl;
          p[f].x = alpha * (_p[f].x + ldt * gdp_lf[f].x);
          p[f].y = alpha * (_p[f].y + ldt * gdp_lf[f].y);
          VectorValue delta = {p[f].x - _p[f].x,
                               p[f].y - _p[f].y};
          diff_p = std::max(diff_p, _normY(delta));
        }
        std::swap(p, _p);
        _u = combination(1.0, _I, -1.0 / lambda, div(_p));
        if (!_color)
        {
          for (VertexIndex i = 0; i < _u.size(); ++i)
            _u[i][2] = _u[i][1] = _u[i][0];
        }
        Scalar new_E = computeEnergyTV();
        trace.info() << "TV( u ) = " << new_E << " ";
        trace.info() << "Iter n=" << (n++) << " diff_p=" << diff_p
                     << " tol=" << tol << std::endl;
        if ((last_E >= 0.0) && (new_E > last_E))
          break;
        else
          last_E = new_E;
      } while ((diff_p > tol) && (n < N));
      return diff_p;
    }

    // equal_strategy:
    // 0: do nothing
    // 1: subdivide all
    // 2: flip all
    // 3: flip all only if flipped = 0
    // @return either (nbflip,0), (nbflip, nbsub_eq) or (nbflip, nbflip_eq).
    std::pair<Integer, Integer>
    onePass(Scalar &total_energy, int equal_strategy = 0)
    {
      Integer nbflipped = 0;
      Integer nbequal = 0;
      total_energy = 0;
      std::vector<Arc> Q_process;
      std::swap(_Queue, Q_process);
      // Taking care of first pass
      if (Q_process.size() == 0)
        for (Arc a = 0; a < T.nbArcs(); ++a)
          if (_flippable[a])
            Q_process.push_back(a);
      // Processing arcs
      for (Arc a : Q_process)
      {
        int update = updateArc(a);
        if (update > 0)
          nbflipped++;
        else if (update == 0)
          _Q_equal.push_back(a);
      }
      total_energy = getEnergyTV();
      const Scalar perimeter = total_energy / (255.0 * (_color ? sqrt(3.0) : 1.0));
      trace.info() << "TV( u ) = " << total_energy
                   << " Perimeter = " << perimeter
                   << " nbflipped=" << nbflipped
                   << "/" << Q_process.size();
      if (equal_strategy == 1)
      {
        nbequal = subdivide(_Q_equal);
        trace.info() << " nbsubequal=" << nbequal
                     << "/" << _Q_equal.size();
        _Q_equal.clear();
      }
      else if (equal_strategy == 2)
      {
        nbequal = flipEqual(_Q_equal);
        trace.info() << " nbflipequal=" << nbequal
                     << "/" << _Q_equal.size();
        _Q_equal.clear();
      }
      else if ((equal_strategy == 3) && (nbflipped == 0))
      {
        nbequal = flipEqual(_Q_equal);
        trace.info() << " nbflipequal=" << nbequal
                     << "/" << _Q_equal.size();
        _Q_equal.clear();
      }
      else if (((equal_strategy >= 4) || (equal_strategy <= 5)) && (nbflipped == 0))
      {
        nbequal = flipEqualWithProb(_Q_equal, 0.5);
        trace.info() << " nbflipequalP=" << nbequal
                     << "/" << _Q_equal.size();
        _Q_equal.clear();
      }
      trace.info() << std::endl;
      return std::make_pair(nbflipped, nbequal);
    }

    template <typename Range>
    Integer flipEqual(const Range &range)
    {
      Integer nbflip = 0;
      for (Arc a : range)
      {
        int update = updateArc(a);
        if (update == 0)
        {
          // Save arcs that may be affected.
          queueSurroundingArcs(a);
          T.flip(a);
          nbflip++;
        }
      }
      return nbflip;
    }

    template <typename Range>
    Integer flipEqualWithProb(const Range &range, double p)
    {
      Integer nbflip = 0;
      for (Arc a : range)
      {
        int update = updateArc(a);
        if (update == 0)
        {
          // Put arc back to potentially process it again.
          _Queue.push_back(a);
          _Queue.push_back(T.opposite(a));
          if (randomUniform() < p)
          {
            // Save arcs that may be affected.
            queueSurroundingArcs(a);
            T.flip(a);
            nbflip++;
          }
        }
      }
      return nbflip;
    }

    template <typename Range>
    Integer subdivide(const Range &range)
    {
      Scalar energy = 0.0;
      Integer nbsubdivided = 0;
      for (Arc a : range)
      {
        int update = updateArc(a);
        if (update == 0)
        {
          VertexRange P = T.verticesOfFacesAroundArc(a);
          // Allow one level of subdivision.
          if (std::max(std::max(P[0], P[1]),
                       std::max(P[2], P[3])) >= _nbV)
            continue;
          // Save arcs that may be affected.
          queueSurroundingArcs(a);
          // Remember faces.
          const Face f012 = T.faceAroundArc(a);
          const Face f023 = T.faceAroundArc(T.opposite(a));
          Scalar Ebefore = energyTV(f012) + energyTV(f023);
          Scalar Eafter = 0.0;
          RealPoint B = (T.position(P[0]) + T.position(P[1]) + T.position(P[2]) + T.position(P[3])) * 0.25;
          VertexIndex v = T.split(a, B);
          FaceRange F = T.facesAroundVertex(v);
          Face new_f = _tv_per_triangle.size();
          _tv_per_triangle.resize(new_f + 2);
          _p.resize(new_f + 2);
          for (Face f : F)
            Eafter += computeEnergyTV(f);
          _tv_energy += Eafter - Ebefore;

          Value V = (_u[P[0]] + _u[P[1]] + _u[P[2]] + _u[P[3]]) * 0.25;
          Value VI = (_I[P[0]] + _I[P[1]] + _I[P[2]] + _I[P[3]]) * 0.25;
          _u.push_back(V);
          _I.push_back(VI);
          ++nbsubdivided;
        }
      }
      return nbsubdivided;
    }

    std::vector<Arc> queueAndSortArcs()
    {
      std::vector<Arc> tv_arcs;
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        const Face f1 = T.faceAroundArc(a);
        const Arc opp_a = T.opposite(a);
        const Face f2 = T.faceAroundArc(opp_a);
        if ((T.head(a) > T.head(opp_a)) && (f1 != T.INVALID_FACE) && (f2 != T.INVALID_FACE))
          if (T.isMergeable(a))
            tv_arcs.push_back(a);
      }
      trace.info() << "[queueAndSortArcs] sorting #arcs=" << tv_arcs.size()
                   << std::endl;
      computeEnergyTV();
      std::sort(tv_arcs.begin(), tv_arcs.end(),
                [&](Arc a1, Arc a2) -> bool
                {
                  const Face f_a1 = T.faceAroundArc(a1);
                  const Face f_a1_opp = T.faceAroundArc(T.opposite(a1));
                  const Face f_a2 = T.faceAroundArc(a2);
                  const Face f_a2_opp = T.faceAroundArc(T.opposite(a2));
                  return (energyTV(f_a1) + energyTV(f_a1_opp)) < (energyTV(f_a2) + energyTV(f_a2_opp));
                });
      return tv_arcs;
    }

    void markSurroundingArcs(const Arc a, std::set<Arc> &M)
    {
      const Vertex v1 = T.head(a);
      const Vertex v2 = T.tail(a);
      auto A1 = T.outArcs(v1);
      auto A2 = T.outArcs(v2);
      for (Arc a1 : A1)
      {
        M.insert(a1);
        M.insert(T.opposite(a1));
        M.insert(T.next(a1));
        M.insert(T.opposite(T.next(a1)));
      }
      for (Arc a2 : A2)
      {
        M.insert(a2);
        M.insert(T.opposite(a2));
        M.insert(T.next(a2));
        M.insert(T.opposite(T.next(a2)));
      }
    }

    VertexRange localNeighborhood(Arc a) const
    {
      const Arc a_next = T.next(a);
      const Arc opp_a_next = T.next(T.opposite(a));
      const Vertex v3 = T.head(a_next);
      const Vertex v4 = T.head(opp_a_next);
      VertexRange V;
      Arc current_a = a_next;
      while (T.head(current_a) != v4)
      {
        V.push_back(T.head(current_a));
        current_a = T.next(T.opposite(current_a));
      }
      current_a = opp_a_next;
      while (T.head(current_a) != v3)
      {
        V.push_back(T.head(current_a));
        current_a = T.next(T.opposite(current_a));
      }
      return V;
    }

    bool isLocallyMergeableAt(Arc a, RealPoint p) const
    {
      auto V = localNeighborhood(a);
      bool ok = true;
      for (int i = 0; ok && (i < V.size()); ++i)
      {
        bool cw = det(T.position(V[i]) - p,
                      T.position(V[(i + 1) % V.size()]) - p) >= 0;
        if (!cw)
          ok = false;
      }
      return ok;
    }

    void simplify(Scalar compression)
    {
      // Objective.
      VertexIndex nb_zip = (VertexIndex)ceil(T.nbVertices() * compression);
      // We need first to sort arcs according to their energyTV, the lower, the first to merge.
      trace.info() << "T=" << T
                   << " valid=" << (T.heds().isValid(false) ? "OK" : "ERROR")
                   << " validT=" << (T.heds().isValidTriangulation() ? "OK" : "ERROR")
                   << std::endl;
      Scalar limit = 0.1;
      while (T.nbVertices() > nb_zip)
      {
        std::set<Arc> M;
        std::vector<Arc> tv_arcs = queueAndSortArcs();
        int current = 0;
        int nb_merged = 0;
        while (current <= (int)ceil(limit * tv_arcs.size()))
        {
          if (T.nbVertices() <= nb_zip)
            break;
          if (current >= tv_arcs.size())
            break;
          trace.info() << "current=" << current << " #V=" << T.nbVertices()
                       << " / " << nb_zip << std::endl;
          const Arc a = tv_arcs[current++];
          if (M.find(a) != M.end())
            continue;
          if (a >= T.nbArcs())
            continue;
          // if ( T.head( a ) >= T.nbVertices() ) continue;
          // if ( T.tail( a ) >= T.nbVertices() ) continue;
          // const Face   f1 = T.faceAroundArc( a );
          // if ( f1 >= T.nbFaces() || f1 == T.INVALID_FACE ) continue;
          const Arc opp_a = T.opposite(a);
          if (opp_a >= T.nbArcs())
            continue;
          // VertexRange P = T.verticesAroundArc( a );
          // if ( ! isConvex( P ) )               continue;

          // const Face   f2 = T.faceAroundArc( opp_a );
          // if ( f2 >= T.nbFaces() || f2 == T.INVALID_FACE ) continue;
          // if ( T.head( a ) < T.head( opp_a ) ) continue;
          RealPoint p = 0.5 * (T.position(T.head(a)) + T.position(T.tail(a)));
          if (isLocallyMergeableAt(a, p) && T.isMergeable(a))
          {
            const Vertex rv = T.head(a);
            markSurroundingArcs(a, M);
            T.merge(a, p);
            _u[rv] = _u.back();
            _u.pop_back();
            _I[rv] = _I.back();
            _I.pop_back();
            nb_merged++;
            trace.info() << "T=" << T
                         // 		 << " valid=" << ( T.heds().isValid( false ) ? "OK" : "ERROR" )
                         // 		 << " validT=" << ( T.heds().isValidTriangulation() ? "OK" : "ERROR" )
                         << std::endl;
          }
        } // while ( current <= (int) ceil( 0.1 * tv_arcs.size() ) ) {
        if (nb_merged == 0)
          limit *= 2.0;
        if (current >= tv_arcs.size())
          break;
      } // while ( T.nbVertices() > nb_zip ) {
    }

    // -------------------------- Regularization services --------------------
  public:
    /// Initializes the regularization process for vectorization.
    void regularizeContours(Scalar max_dt = 0.001, int max_iter = 10)
    {
      computeArcDissimilarities();
      initContours();
      for (int n = 0; n < max_iter; ++n)
      {
        std::vector<Scalar> former_t = _t;
        updateBarycenters();
        updateContours();
        Scalar dt = enforceContoursArea(former_t);
        if (dt < max_dt)
          break;
        std::cout << n << " dt = " << dt << " max_dt = " << max_dt
                  << std::endl;
      }
    }

    /// Initializes the regularization process for vectorization.
    void initContours()
    {
      _arc_updatable.resize(T.nbArcs());
      _face_updatable.resize(T.nbFaces());
      // Computes the barycenter for each valid face.
      _b.resize(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        VertexRange V = T.verticesAroundFace(f);
        _b[f] = (T.position(V[0]) + T.position(V[1]) + T.position(V[2])) / 3.0;
        auto arcs = arcsAroundFace(f);
        // _face_updatable[ f ] = _flippable[ arcs[ 0 ] ]
        //   || _flippable[ arcs[ 1 ] ]
        //   || _flippable[ arcs[ 2 ] ];
        _face_updatable[f] =
            (_arc_dissimilarities[arcs[0]] != 0.0) || (_arc_dissimilarities[arcs[1]] != 0.0) || (_arc_dissimilarities[arcs[2]] != 0.0);
        // if ( _face_updatable[ f ] )
        // _faces_to_update.push_back( f );
      }
      // Computes the contour intersections at each arc.
      _t.resize(T.nbArcs());
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        Arc opp_a = T.opposite(a);
        if (T.isArcBoundary(a) || T.isArcBoundary(opp_a))
        {
          _t[a] = 0.5;
          continue;
        }
        Face f_a = T.faceAroundArc(a);
        Face f_opp_a = T.faceAroundArc(opp_a);
        if (!_face_updatable[f_a] && !_face_updatable[f_opp_a])
        {
          _arc_updatable[a] = false;
          _arc_updatable[opp_a] = false;
          _t[a] = 0.5;
        }
        else
        {
          RealPoint B = T.position(T.head(a));
          RealPoint A = T.position(T.head(opp_a));
          auto I = intersect(A, B, _b[f_a], _b[f_opp_a]);
          _arc_updatable[a] = true;
          _arc_updatable[opp_a] = true;
          _t[a] = std::min(0.999, std::max(0.001, I.first));
        }
        // std::cout << "t[" << a << "]=" << _t[ a ] << std::endl;
      }
      // Computes the areas associated with each vertex
      _A.resize(T.nbVertices());
      for (Vertex v = 0; v < T.nbVertices(); ++v)
      {
        // std::cout << v << " area=" << areaAtVertex( v ) << std::endl;
        // _A[ v ] = areaAtVertex( v );
        _A[v] = (Scalar)T.degree(v) / 6.0; //areaAtVertex( v );
      }
    }

    bool isArcContour(const Arc a) const
    {
      return _u[T.head(a)] != _u[T.tail(a)];
    }
    Scalar computeArcDissimilarity(const Arc a) const
    {
      Value diff = _u[T.head(a)] - _u[T.tail(a)];
      return diff.dot(diff);
      // return ( _u[ T.head( a ) ] - _u[ T.tail( a ) ] ).norm();
    }

    /// @param[in] f any valid face.
    /// @return the three arcs belonging to the face f.
    std::array<Arc, 3> arcsAroundFace(const Face f) const
    {
      VertexRange V = T.verticesAroundFace(f);
      std::array<Arc, 3> arcs;
      arcs[0] = T.arc(V[0], V[1]);
      arcs[1] = T.next(arcs[0]);
      arcs[2] = T.next(arcs[1]);
      return arcs;
    }

    /// @param[out] dissim the (dis)similarity coefficient for each arc
    /// @param[out] arcs the three arcs belonging to face \a f.
    /// @param[in] f any valid face
    ///
    /// @return the arc whose vertices have similar values if it
    /// exists, or -1 otherwise.
    int arcDissimilaritiesAtFace(std::array<Scalar, 3> &dissim,
                                 std::array<Arc, 3> &arcs,
                                 Face f) const
    {
      arcs = arcsAroundFace(f);
      for (int i = 0; i < 3; ++i)
      {
        dissim[i] = computeArcDissimilarity(arcs[i]);
      }
      // Sort similarities...
      int m = std::min_element(dissim.begin(), dissim.end()) - dissim.begin();
      bool id = ((2 * dissim[m]) < dissim[(m + 1) % 3]) && ((2 * dissim[m]) < dissim[(m + 2) % 3]);
      if (id)
        dissim[m] = 0;
      return id ? m : -1;
    }

    /// Computes the barycenters for each triangle. The idea is to
    /// count only arcs with different values.
    void updateBarycenters()
    {
      // for ( Face f : _faces_to_update ) {
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        if (!_face_updatable[f])
          continue;
        auto arcs = arcsAroundFace(f);
        Scalar w = 0.0;
        RealPoint B = RealPoint::zero;
        // std::array<Scalar,3> s;
        // int m = arcDissimilarities( s, arcs );
        for (int i = 0; i < 3; ++i)
        {
          const Scalar s = _arc_dissimilarities[arcs[i]];
          B += s * contourPoint(arcs[i]);
          w += s;
        }
        if (w > 0.0)
          _b[f] += 0.5 * (B / w - _b[f]);
      }
    }

    /// Computes the contour points for each arc from the barycenters.
    /// @return the maximum displacement.
    void updateContours()
    {
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        if (!_arc_updatable[a])
          continue;
        Arc opp_a = T.opposite(a);
        // if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) continue;
        Face f_a = T.faceAroundArc(a);
        Face f_opp_a = T.faceAroundArc(opp_a);
        RealPoint B = T.position(T.head(a));
        RealPoint A = T.position(T.head(opp_a));
        if ((f_a != T.INVALID_FACE) && (f_opp_a != T.INVALID_FACE))
        {
          auto I = intersect(A, B, _b[f_a], _b[f_opp_a]);
          Scalar t = std::min(0.999, std::max(0.001, I.first));
          _t[a] += 0.5 * (t - _t[a]);
        }
        else
        {
          _t[a] = 0.5;
        }
      }
    }

    /// Moves contour points to keep constant areas
    /// @return the maximum displacement.
    Scalar enforceContoursArea(const std::vector<Scalar> &prev_t)
    {
      Scalar max_t = 0.0;
      // Update t to keep volume constant
      for (Vertex v = 0; v < T.nbVertices(); ++v)
      {
        auto out_arcs = T.outArcs(v);
        for (Arc a : out_arcs)
        {
          if (!_arc_updatable[a])
            continue;
          // Checks the evolution of the arc area.
          Scalar ratio = areaAtArc(a) / (_A[v] / T.degree(v));
          if (ratio <= 0.001)
          {
            trace.warning() << "Negative ratio " << ratio
                            << " for area[ " << v << " ]=" << _A[v]
                            << " at_arc=" << areaAtArc(a)
                            << " d=" << T.degree(v)
                            << " A/d=" << (_A[v] / T.degree(v))
                            << std::endl;
            ratio = 0.001;
          }
          // _t[ a ] /= ( 0.2 + 0.8 * ratio );
          _t[a] = 0.5 * _t[a] * (1.0 + 1.0 / ratio);
        }
      }
      // Averages movements on each edge.
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        if (!_arc_updatable[a])
          continue;
        const Arc opp_a = T.opposite(a);
        // if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) continue;
        if (T.head(a) < T.head(opp_a))
          continue;
        Scalar t = 0.5 * (_t[a] + (1.0 - _t[opp_a]));
        t = std::min(0.999, std::max(0.001, t));
        Scalar dt = fabs(t - prev_t[a]);
        max_t = std::max(max_t, dt);
        _t[a] = t;
        _t[opp_a] = 1.0 - t;
      }
      return max_t;
    }

    /// The area at an arc \a a is the area associated to the vertex
    /// tail of \a a within the face of \a a. It corresponds to the
    /// zone between the tail, the barycenter and the two intersection
    /// points at the edges containing the tail vertex.
    Scalar areaAtArc(Arc a) const
    {
      const Face f = T.faceAroundArc(a);
      const Arc opp_a = T.opposite(a);
      if (T.isArcBoundary(a) || T.isArcBoundary(opp_a))
        return 1.0 / 6.0;

      const Arc a2 = T.next(opp_a);
      const Face f2 = T.faceAroundArc(a2);
      const RealPoint B = T.position(T.head(a));
      const RealPoint A = T.position(T.head(opp_a));
      const Scalar t = _t[a];
      return -0.5 * (det(t * (B - A), _b[f] - A) + det(_b[f2] - A, t * (B - A)));
      // const Face   f = T.faceAroundArc( a );
      // const Arc   an = T.next( a );
      // const Arc  ann = T.next( an );
      // const RealPoint  B = T.position( T.head( a ) );
      // const RealPoint  C = T.position( T.head( an ) );
      // const RealPoint  A = T.position( T.head( ann ) );
      // const Arc   a2 = T.opposite( ann );
      // const Scalar t = _t[ a ];
      // const Scalar u = _t[ a2 ];
      // return - 0.5 * ( det( t * ( B - A ), _b[ f ] - A )
      // 		       +  det( _b[ f ] - A , u * ( C - A ) ) );
    }

    Scalar areaAtVertex(const Vertex v) const
    {
      auto out_arcs = T.outArcs(v);
      Scalar area = 0;
      for (Arc a : out_arcs)
        area += areaAtArc(a);
      return area;
    }

    /// Given two straight lines (AB) and (CD), returns their
    /// intersection point as two scalars (t,u) such that:
    /// I = A + t AB = C + u CD
    static std::pair<Scalar, Scalar> intersect(const RealPoint &A, const RealPoint &B,
                                               const RealPoint &C, const RealPoint &D)
    {
      RealVector AB = B - A;
      RealVector DC = C - D;
      RealVector AC = C - A;
      Scalar d = det(AB, DC);
      if (d == 0)
        return std::make_pair(0.5, 0.5);
      Scalar t = (DC[1] * AC[0] - DC[0] * AC[1]) / d;
      Scalar u = (AB[0] * AC[1] - AB[1] * AC[0]) / d;
      return std::make_pair(t, u);
    }

    /// Computes information per face
    void compute2ndOrderInformation(Scalar discontinuities)
    {
      _finfos.resize(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        auto V = T.verticesAroundFace(f);
        _finfos[f].grad = grad(f, _u);
        _finfos[f].v = (_u[V[0]] + _u[V[1]] + _u[V[2]]) / 3.0;
      }
      // We need first to sort faces according to their energyTV.
      std::vector<Face> tv_faces(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
        tv_faces[f] = f;
      std::sort(tv_faces.begin(), tv_faces.end(),
                [&](Face f1, Face f2) -> bool
                // { return ( tvT.energyTV( f1 ) ) > ( tvT.energyTV( f2 ) ); }
                // { return ( tvT.aspectRatio( f1 ) )
                //     > ( tvT.aspectRatio( f2 ) ); }
                { return (energyTV(f1) * diameter(f1)) > (energyTV(f2) * diameter(f2)); }
                // { return ( tvT.energyTV( f1 ) * tvT.aspectRatio( f1 ) )
                //     > ( tvT.energyTV( f2 ) * tvT.aspectRatio( f2 ) ); }
      );
      Scalar Etv = getEnergyTV();
      Scalar Ctv = 0.0;
      Scalar Otv = Etv * discontinuities;
      for (int i = 0; i < tv_faces.size(); ++i)
      {
        Face f = tv_faces[i];
        Ctv += energyTV(f);
        _finfos[f].order = i;
        _finfos[f].rel_tv = energyTV(f) / (Etv / tv_faces.size());
        _finfos[f].cumul_tv = Ctv / Etv;
        _finfos[f].disc = Ctv / Otv;
        _finfos[f].cont = (Ctv - Otv) / Etv;
      }
    }

    // Computes approximation of u at each vertex.
    void computeMetrics()
    {
      const int d = _color ? 3 : 1;
      _metrics.resize(T.nbFaces());
      _p.resize(T.nbFaces());
      _b.resize(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        auto V = T.verticesAroundFace(f);
        _metrics[f] = Metric(T.position(V[0]),
                             T.position(V[1]),
                             T.position(V[2]));
        _p[f] = grad(f, _u);
        _b[f] = (T.position(V[0]) + T.position(V[1]) + T.position(V[2])) / 3.0;
      }
      trace.info() << "Finishing computing metrics" << std::endl;
      // computeUApproximations();
      computeUApproximationsAtFaces();
      trace.info() << "Finishing computing U approx" << std::endl;
    }

    Value evalMetrics(const RealPoint &p) const
    {
      const int d = _color ? 3 : 1;
      //      trace.info() << "[evalMetrics] " << p
      //           << " #V=" << T.nbVertices()
      //           << " F=" << T.nbFaces()
      //           << " d=" << d << std::endl;
      const Scalar r = 1.5;
      Point lo((int)(round(p[0]) - r), (int)(round(p[1]) - r));
      Point hi((int)(round(p[0]) + r), (int)(round(p[1]) + r));
      lo = lo.sup(_lo);
      hi = hi.inf(_up);
      Value R;
      Scalar wsum = 0;
      // Find surrounding faces.
      std::set<Face> faces;
      Domain local(lo, hi);
      for (auto qi : local)
      {
        Vertex v = linearize(qi);
        auto F = T.facesAroundVertex(v);
        for (auto f : F)
          faces.insert(f);
        // Scalar       w = fM( ( p - RealPoint( qi[ 0 ], qi[ 1 ] ) ).norm() );
        // for ( int m = 0; m < d; ++m ) {
        //   R[ m ] += w * _u_approx[ v ][ m ]( p );
        // }
        // wsum   += w;
      }
      // Compute value by metrized partition of unity
      for (auto f : faces)
      {
        RealVector bp = _b[f] - p;
        Scalar dist2 = _metrics[f](bp, bp);
        Scalar w = fM(sqrt(dist2));
        for (int m = 0; m < d; ++m)
        {
          Scalar v_stiff = _u_approx_f[f][m](p, 1.0, 0.95);
          Scalar v_grad = _u_approx_f[f][m](p);
          // R[ m ] += w * _u_approx_f[ f ][ m ]( p );
          R[m] += w * (w * v_stiff + (1.0 - w) * v_grad);
        }
        wsum += w;
      }
      R /= wsum;
      if (d == 1)
      {
        R[1] = R[0];
        R[2] = R[0];
      }
      return R;
    }

    // Computes approximation of u at each face.
    void computeUApproximationsAtFaces()
    {
      const int d = _color ? 3 : 1;
      _u_approx_f.resize(T.nbFaces());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        auto vertices = T.verticesAroundFace(f);
        std::vector<RealPoint> X;
        std::array<std::vector<Scalar>, 3> V;
        for (auto w : vertices)
        {
          X.push_back(T.position(w));
          for (int m = 0; m < d; ++m)
            V[m].push_back(_u[w][m]);
        }
        // Fit values in the least-square sense for each RGB channel.
        for (int m = 0; m < d; ++m)
          _u_approx_f[f][m].fit(X, V[m]);
      }
    }

    // Computes approximation of u at each vertex.
    void computeUApproximations()
    {
      const int d = _color ? 3 : 1;
      _u_approx.resize(T.nbVertices());
      for (Vertex v = 0; v < T.nbVertices(); ++v)
      {
        auto arcs = T.outArcs(v);
        std::vector<RealPoint> X;
        std::array<std::vector<Scalar>, 3> V;
        X.push_back(T.position(v));
        for (int m = 0; m < d; ++m)
          V[m].push_back(_u[v][m]);
        for (auto a : arcs)
        {
          Vertex w = T.head(a);
          X.push_back(T.position(w));
          for (int m = 0; m < d; ++m)
            V[m].push_back(_u[w][m]);
        }
        // Fit values in the least-square sense for each RGB channel.
        for (int m = 0; m < d; ++m)
          _u_approx[v][m].fit(X, V[m]);
      }
    }

    /// s close to 0 means that the coefficient should be important.
    static Scalar fM(Scalar s)
    {
      //return std::max( 1.0 - s/2, 0.0 ); //exp( -s*s*2 );
      // between -s^2 and -2s^2 is a good trade-off
      return exp(-s * s * 2);
    }
    /// s close to 0 means that the coefficient should be important.
    static Scalar fPOU(Scalar s)
    {
      return std::max(1.0 - s, 0.0); //exp( -s*s*2 );
      // return exp( -s*s*4 );
    }

    Value evalPOU(const RealPoint &p) const
    {
      const int d = _color ? 3 : 1;
      const Scalar r = 3.0;
      Point lo((int)(round(p[0]) - r), (int)(round(p[1]) - r));
      Point hi((int)(round(p[0]) + r), (int)(round(p[1]) + r));
      lo = lo.sup(_lo);
      hi = hi.inf(_up);
      // Start scanning of values for POU
      Domain local(lo, hi);
      Value R;
      Scalar wsum = 0;
      for (auto qi : local)
      {
        RealPoint q(qi[0], qi[1]);
        Scalar pq = (q - p).norm();
        if (pq <= r)
        {
          Integer v = linearize(qi);
          Scalar w = fPOU(pq);
          for (int m = 0; m < d; ++m)
            R[m] += w * _u_approx[v][m](p);
          wsum += w;
        }
      }
      R /= wsum;
      if (d == 1)
      {
        R[1] = R[0];
        R[2] = R[0];
      }
      return R;
    }

    /// @param p any point
    /// @return the corresponding index in the 1-dimensional array.
    Integer linearize(Point p) const
    {
      return p[1] * _width + p[0];
    }

    /// @return 'true' iff the given barycentric coordinates \a bc
    /// indicates a point inside or on the boundary of the triangle.
    bool isInTriangle(const Face f, const RealPoint &p) const
    {
      const auto V = T.verticesAroundFace(f);
      RealPoint b[3] = {T.position(V[0]),
                        T.position(V[1]),
                        T.position(V[2])};
      Value bc = {det(b[1] - p, b[2] - p),
                  det(b[2] - p, b[0] - p),
                  det(b[0] - p, b[1] - p)};
      bc /= (bc[0] + bc[1] + bc[2]);
      return (0 <= bc[0]) && (bc[0] <= 1) && (0 <= bc[1]) && (bc[1] <= 1) && (0 <= bc[2]) && (bc[2] <= 1);
    }

    /// @return 'true' iff the given barycentric coordinates \a bc
    /// indicates a point inside or on the boundary of the triangle.
    bool isApproximatelyInTriangle(const Face f, const RealPoint &p,
                                   Scalar tol = 0.05) const
    {
      const auto V = T.verticesAroundFace(f);
      RealPoint b[3] = {T.position(V[0]),
                        T.position(V[1]),
                        T.position(V[2])};
      Value bc = {det(b[1] - p, b[2] - p),
                  det(b[2] - p, b[0] - p),
                  det(b[0] - p, b[1] - p)};
      bc /= (bc[0] + bc[1] + bc[2]);
      return (-tol <= bc[0]) && (bc[0] <= 1.0 + tol) && (-tol <= bc[1]) && (bc[1] <= 1.0 + tol) && (-tol <= bc[2]) && (bc[2] <= 1.0 + tol);
    }

    /// Precompute arcs similarities
    void computeArcDissimilarities()
    {
      _arc_dissimilarities.resize(T.nbArcs());
      for (Face f = 0; f < T.nbFaces(); ++f)
      {
        std::array<Scalar, 3> s;
        std::array<Arc, 3> arcs;
        ;
        const int m = arcDissimilaritiesAtFace(s, arcs, f);
        _arc_dissimilarities[arcs[0]] = s[0];
        _arc_dissimilarities[arcs[1]] = s[1];
        _arc_dissimilarities[arcs[2]] = s[2];
      }
      for (Arc a = 0; a < T.nbArcs(); ++a)
      {
        Arc oa = T.opposite(a);
        _arc_dissimilarities[a] = std::max(_arc_dissimilarities[a],
                                           _arc_dissimilarities[oa]);
      }
    }
  };

  // Useful function for viewing triangulations.

  /**
     This class is intended for visualizing Affine Valued
     triangulation with CAIRO.
  */
  class CairoViewerTV : public CairoViewer<Z2i::Space>
  {
  public:
    typedef CairoViewer<Z2i::Space> Base;
    typedef CairoViewerTV Self;
    typedef TVTriangulation TVT;
    typedef TVT::Value Value;
    typedef TVT::VectorValue VectorValue;
    typedef TVT::Triangulation Triangulation;
    typedef TVT::VertexIndex VertexIndex;
    typedef TVT::Vertex Vertex;
    typedef TVT::Arc Arc;
    typedef TVT::Face Face;
    typedef TVT::VertexRange VertexRange;
    typedef TVT::Scalar Scalar;
    typedef TVT::Point Point;
    typedef TVT::RealPoint RealPoint;
    typedef Z2i::Domain Domain;

    using Base::_am;
    using Base::_color;
    using Base::_draw_domain;
    using Base::_xf;
    using Base::_yf;
    using Base::i;
    using Base::ij;
    using Base::j;
    using Base::x;
    using Base::xy;
    using Base::y;

    /// Contains information per boundary pixel for 2nd order reconstruction.
    struct PixelInformation
    {
      // VectorValue grad; ///< grad in pixel coordinates
      Value v; ///< value at pixel
      // Value      min_v; ///< min value in face around pixel
      // Value      max_v; ///< max value in face around pixel
      Scalar crisp; ///< 0 is linear gradient wrt to distance d, otherwise straightened the contour around.
      bool operator==(const PixelInformation &other) const
      {
        return (v == other.v) && (crisp == other.crisp);
      }
      bool operator!=(const PixelInformation &other) const
      {
        return (v != other.v) || (crisp != other.crisp);
      }
    };

    typedef ImageContainerBySTLVector<Domain, Arc> ArcImage;
    typedef ImageContainerBySTLVector<Domain, PixelInformation>
        PixelInformationImage;
    typedef ImageContainerBySTLVector<Domain, Value> OutputImage;

  public:
    Point _lo;
    Point _up;

    PixelInformationImage _pixinfo;       ///< image storing pixel values/crispness.
    PixelInformationImage _pixinfo_disc;  ///< image storing pix values/crispness.
    ArcImage _arcimage_disc;              ///< image storing discontinuities
    ArcImage _arcimage_simi;              ///< image storing similarities
    OutputImage _output;                  ///< image storing output pixel
    std::function<Point(RealPoint)> _dig; ///< discretization function (input image domain -> output draw image domain)

    /**
       Constructor. 
    */
    CairoViewerTV(int x0, int y0, int width, int height,
                  double xfactor = 1.0, double yfactor = 1.0,
                  int shading = 0,
                  bool color = true,
                  double disc_stiffness = 0.5,
                  double disc_amplitude = 0.75)
        : Base(x0, y0, width, height, xfactor, yfactor, shading, color,
               disc_stiffness, disc_amplitude),
          _pixinfo(Domain(Point(), Point())),
          _pixinfo_disc(Domain(Point(), Point())),
          _arcimage_disc(Domain(Point(), Point())),
          _arcimage_simi(Domain(Point(), Point())),
          _output(Domain(Point(), Point()))
    {
      _lo = Point(x0, y0);
      _up = Point(width, height);
      _pixinfo = PixelInformationImage(_draw_domain);
      _pixinfo_disc = PixelInformationImage(_draw_domain);
      _arcimage_disc = ArcImage(_draw_domain);
      _arcimage_simi = ArcImage(_draw_domain);
      _output = OutputImage(_draw_domain);
      _dig = [&](const RealPoint &p)
      {
        RealPoint q = ij(p);
        return Point(floor(q[0]), floor(q[1]));
      };
    }

    /// Destructor.
    ~CairoViewerTV() {}

    void viewTVTLinearGradientTriangle(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      drawLinearGradientTriangle(RealPoint(a[0], a[1]),
                                 RealPoint(b[0], b[1]),
                                 RealPoint(c[0], c[1]),
                                 tvT.u(V[0]),
                                 tvT.u(V[1]),
                                 tvT.u(V[2]));
    }
    void viewTVTNonLinearGradientTriangle(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      drawNonLinearGradientTriangle(RealPoint(a[0], a[1]),
                                    RealPoint(b[0], b[1]),
                                    RealPoint(c[0], c[1]),
                                    tvT.u(V[0]),
                                    tvT.u(V[1]),
                                    tvT.u(V[2]));
    }
    void viewTVTGouraudTriangle(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      drawGouraudTriangle(RealPoint(a[0], a[1]),
                          RealPoint(b[0], b[1]),
                          RealPoint(c[0], c[1]),
                          tvT.u(V[0]),
                          tvT.u(V[1]),
                          tvT.u(V[2]));
    }
    void viewTVTFlatTriangle(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      Value val = tvT.u(V[0]) + tvT.u(V[1]) + tvT.u(V[2]);
      val /= 3.0;
      drawFlatTriangle(RealPoint(a[0], a[1]),
                       RealPoint(b[0], b[1]),
                       RealPoint(c[0], c[1]), val);
    }

    void viewTVTTriangleDiscontinuity(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      Value val = {255, 0, 0};
      drawFlatTriangle(RealPoint(a[0], a[1]),
                       RealPoint(b[0], b[1]),
                       RealPoint(c[0], c[1]), val);
    }

    void viewTVTPartitionOfUnityTriangle(TVT &tvT, Face f)
    {
      VertexRange V = tvT.T.verticesAroundFace(f);
      RealPoint a = tvT.T.position(V[0]);
      RealPoint b = tvT.T.position(V[1]);
      RealPoint c = tvT.T.position(V[2]);
      auto fa = tvT.uApprox(V[0]);
      auto fb = tvT.uApprox(V[1]);
      auto fc = tvT.uApprox(V[2]);
      drawPartitionOfUnityTriangle(RealPoint(a[0], a[1]),
                                   RealPoint(b[0], b[1]),
                                   RealPoint(c[0], c[1]),
                                   fa, fb, fc);
    }

    /**
       Displays the AVT with flat or Gouraud shading.
    */
    void view(TVT &tvT)
    {
      if (_shading == 3)
      {
        trace.beginBlock("Compute u approximations");
        // tvT.computeUApproximations();
        tvT.computeMetrics();
        cairo_set_operator(_cr, CAIRO_OPERATOR_OVER);
        cairo_set_line_width(_cr, 0.0);
        cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
        cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
        drawPartitionOfUnity(_draw_domain, tvT);
        trace.endBlock();
      }
      else
      {
        cairo_set_operator(_cr, CAIRO_OPERATOR_ADD);
        cairo_set_line_width(_cr, 0.0);
        cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
        cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
        for (Face f = 0; f < tvT.T.nbFaces(); ++f)
        {
          if (_shading == 1)
            viewTVTGouraudTriangle(tvT, f);
          else if (_shading == 2)
            viewTVTLinearGradientTriangle(tvT, f);
          // else if ( _shading == 3 ) viewTVTPartitionOfUnityTriangle( tvT, f );
          else
            viewTVTFlatTriangle(tvT, f);
        }
      }
    }

    /**
       Displays the Laplacian with flat shading.
    */
    void viewLaplacian(TVT &tvT)
    {
      cairo_set_operator(_cr, CAIRO_OPERATOR_ADD);
      cairo_set_line_width(_cr, 0.0);
      cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
      cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
      auto Lp = tvT.div(tvT.grad(tvT._u));
      for (Face f = 0; f < tvT.T.nbFaces(); ++f)
      {
        VertexRange V = tvT.T.verticesAroundFace(f);
        RealPoint a = tvT.T.position(V[0]);
        RealPoint b = tvT.T.position(V[1]);
        RealPoint c = tvT.T.position(V[2]);
        Scalar val = (tvT._normX(Lp[V[0]]) + tvT._normX(Lp[V[1]]) + tvT._normX(Lp[V[2]])) / 3.0;
        val = std::min(255.0, std::max(0.0, val));
        drawFlatTriangle(RealPoint(a[0], a[1]),
                         RealPoint(b[0], b[1]),
                         RealPoint(c[0], c[1]),
                         Value(val, val, val));
      }
    }

    /**
       Displays the Laplacian with flat shading.
    */
    void viewGradientLaplacian(TVT &tvT)
    {
      cairo_set_operator(_cr, CAIRO_OPERATOR_ADD);
      cairo_set_line_width(_cr, 0.0);
      cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
      cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
      auto Lp = tvT.div(tvT.grad(tvT._u));
      for (Face f = 0; f < tvT.T.nbFaces(); ++f)
      {
        VertexRange V = tvT.T.verticesAroundFace(f);
        RealPoint a = tvT.T.position(V[0]);
        RealPoint b = tvT.T.position(V[1]);
        RealPoint c = tvT.T.position(V[2]);
        Scalar v0 = tvT._normX(Lp[V[0]]);
        Scalar v1 = tvT._normX(Lp[V[1]]);
        Scalar v2 = tvT._normX(Lp[V[2]]);
        Scalar vm = std::min(v0, std::min(v1, v2));
        // drawFlatTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
        // 		  RealPoint( b[ 0 ], b[ 1 ] ),
        // 		  RealPoint( c[ 0 ], c[ 1 ] ),
        // 		  Value( vm, vm, vm ) );
        drawLinearGradientTriangle(RealPoint(a[0], a[1]),
                                   RealPoint(b[0], b[1]),
                                   RealPoint(c[0], c[1]),
                                   Value(v0, v0, v0),
                                   Value(v1, v1, v1),
                                   Value(v2, v2, v2));
      }
      // cairo_set_operator( _cr,  CAIRO_OPERATOR_OVER );
      // cairo_set_line_width( _cr, 1.0 );
      // for ( Face f = 0; f < tvT.T.nbFaces(); ++f ) {
      // 	VertexRange V = tvT.T.verticesAroundFace( f );
      // 	RealPoint       a = tvT.T.position( V[ 0 ] );
      // 	RealPoint       b = tvT.T.position( V[ 1 ] );
      // 	RealPoint       c = tvT.T.position( V[ 2 ] );
      // 	Scalar         v0 = tvT._normX( Lp[ V[ 0 ] ] ) / sqrt( 3.0 );
      // 	Scalar         v1 = tvT._normX( Lp[ V[ 1 ] ] ) / sqrt( 3.0 );
      // 	Scalar         v2 = tvT._normX( Lp[ V[ 2 ] ] ) / sqrt( 3.0 );
      // 	Scalar        v01 = std::min( v0, v1 );
      // 	Scalar        v02 = std::min( v0, v2 );
      // 	Scalar        v12 = std::min( v1, v2 );
      // 	drawFlatLine(  RealPoint( a[ 0 ], a[ 1 ] ),
      // 		       RealPoint( b[ 0 ], b[ 1 ] ),
      // 		       Value( v01, v01, v01 ) );
      // 	drawFlatLine(  RealPoint( a[ 0 ], a[ 1 ] ),
      // 		       RealPoint( c[ 0 ], c[ 1 ] ),
      // 		       Value( v02, v02, v02 ) );
      // 	drawFlatLine(  RealPoint( c[ 0 ], c[ 1 ] ),
      // 		       RealPoint( b[ 0 ], b[ 1 ] ),
      // 		       Value( v12, v12, v12 ) );
      // }
    }

    /**
       Displays the AVT with flat or Gouraud shading, and displays a
       set of discontinuities as a percentage of the total energy.
    */
    void view(TVT &tvT, Scalar discontinuities)
    {
      if (_shading == 3)
      {
        trace.beginBlock("Compute u approximations");
        // tvT.computeUApproximations();
        tvT.computeMetrics();
        cairo_set_operator(_cr, CAIRO_OPERATOR_OVER);
        cairo_set_line_width(_cr, 0.0);
        cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
        cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
        drawPartitionOfUnity(_draw_domain, tvT);
        trace.endBlock();
      }
      else
      {
        // We need first to sort faces according to their energyTV.
        std::vector<Face> tv_faces(tvT.T.nbFaces());
        for (Face f = 0; f < tvT.T.nbFaces(); ++f)
          tv_faces[f] = f;
        std::sort(tv_faces.begin(), tv_faces.end(),
                  [&tvT](Face f1, Face f2) -> bool
                  // { return ( tvT.energyTV( f1 ) ) > ( tvT.energyTV( f2 ) ); }
                  // { return ( tvT.aspectRatio( f1 ) )
                  //     > ( tvT.aspectRatio( f2 ) ); }
                  { return (tvT.energyTV(f1) * tvT.diameter(f1)) > (tvT.energyTV(f2) * tvT.diameter(f2)); }
                  // { return ( tvT.energyTV( f1 ) * tvT.aspectRatio( f1 ) )
                  //     > ( tvT.energyTV( f2 ) * tvT.aspectRatio( f2 ) ); }
        );
        Scalar Etv = tvT.getEnergyTV();
        Scalar Ctv = 0.0;
        Scalar Otv = Etv * discontinuities;
        cairo_set_operator(_cr, CAIRO_OPERATOR_ADD);
        cairo_set_line_width(_cr, 0.0);
        cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
        cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
        for (int i = 0; i < tv_faces.size(); ++i)
        {
          Face f = tv_faces[i];
          Ctv += tvT.energyTV(f);
          if (Ctv < Otv)
          { // display discontinuity
            // viewTVTTriangleDiscontinuity( tvT, f );
            viewTVTNonLinearGradientTriangle(tvT, f);
          }
          else
          {
            if (_shading == 1)
              viewTVTGouraudTriangle(tvT, f);
            else if (_shading == 2)
              viewTVTLinearGradientTriangle(tvT, f);
            // else if ( _shading == 3 ) viewTVTPartitionOfUnityTriangle( tvT, f );
            else
              viewTVTFlatTriangle(tvT, f);
          }
        }
      }
    }

    /// @return the normalized tangent of the contour crossing arc \a a.
    RealVector tangentAtArc(TVT &tvT, const Arc a)
    {
      // The tangent is arbitrary defined as parallel to the vector
      // joining the two barycenters and points inward the face
      const Face f = tvT.T.faceAroundArc(a);
      const Arc oa = tvT.T.opposite(a);
      if (!tvT.T.isArcBoundary(oa))
      {
        const Face of = tvT.T.faceAroundArc(oa);
        return (tvT.barycenter(f) - tvT.barycenter(of)).getNormalized();
      }
      else
      {
        const RealPoint v1 = tvT.T.position(tvT.T.head(a));
        const RealPoint v2 = tvT.T.position(tvT.T.tail(a));
        return (tvT.barycenter(f) - 0.5 * (v1 + v2)).getNormalized();
      }
    }

    /// @return the crispness at arc \a a
    Scalar crispnessAtArc(TVT &tvT, const Arc a)
    {
      const Face f = tvT.T.faceAroundArc(a);
      const Arc oa = tvT.T.opposite(a);
      if (!tvT.T.isArcBoundary(oa))
      {
        const Face of = tvT.T.faceAroundArc(oa);
        return 0.5 * (crispness(tvT, of) + crispness(tvT, f));
      }
      else
      {
        return crispness(tvT, f);
      }
    }
    /// @return the crispness at vertex \a v
    Scalar crispnessAtVertex(TVT &tvT, const Vertex v)
    {
      Scalar c = 0.0;
      const auto F = tvT.T.facesAroundVertex(v);
      for (auto f : F)
        c += crispness(tvT, f);
      return c / F.size();
    }

    VectorValue combine(const VectorValue &g1, const VectorValue &g2,
                        Scalar t = 0.5)
    {
      return VectorValue{(1.0 - t) * g1.x + t * g2.x,
                         (1.0 - t) * g1.y + t * g2.y};
    }
    Value combine(const Value &g1, const Value &g2,
                  Scalar t = 0.5)
    {
      return (1.0 - t) * g1 + t * g2;
    }
    Scalar combine(const Scalar &g1, const Scalar &g2,
                   Scalar t = 0.5)
    {
      return (1.0 - t) * g1 + t * g2;
    }
    /// @return the value at an arc \a a.
    Value valueAtArc(TVT &tvT, const Arc a)
    {
      return 0.5 * (tvT.u(tvT.T.head(a)) + tvT.u(tvT.T.tail(a)));
    }

    /// @return the value at an arc \a a.
    std::pair<Value, Value> boundAtArc(TVT &tvT, const Arc a)
    {
      Vertex h = tvT.T.head(a);
      Vertex t = tvT.T.tail(a);
      Value m, M;
      M = m = tvT.u(h);
      M = M.sup(tvT.u(t));
      m = m.inf(tvT.u(t));
      return std::make_pair(m, M);
    }

    /// @return the value at an arc \a a.
    std::pair<Value, Value> boundAtFace(TVT &tvT, const Face f)
    {
      auto V = tvT.T.verticesAroundFace(f);
      Value m, M;
      M = m = tvT.u(V[0]);
      M = M.sup(tvT.u(V[1])).sup(tvT.u(V[2]));
      m = m.inf(tvT.u(V[1])).inf(tvT.u(V[2]));
      return std::make_pair(m, M);
    }
    /// @return the value at the barycenter of a face \a f.
    Value valueAtBarycenter(TVT &tvT, const Face f)
    {
      return tvT.finfo(f).v;
    }
    /// @return the gradient at the barycenter of a face \a f.
    VectorValue gradientAtBarycenter(TVT &tvT, const Face f)
    {
      return tvT.finfo(f).grad;
    }

    /// @return the gradient at an arc \a a.
    VectorValue gradientAtArc(TVT &tvT, const Arc a)
    {
      const Face f = tvT.T.faceAroundArc(a);
      VectorValue gf = tvT.finfo(f).grad;
      const Arc oa = tvT.T.opposite(a);
      if (!tvT.T.isArcBoundary(oa))
      {
        const Face of = tvT.T.faceAroundArc(oa);
        VectorValue gof = tvT.finfo(of).grad;
        return combine(gf, gof, 0.5);
      }
      else
      {
        return gf;
      }
    }

    Scalar crispness(TVT &tvT, const Face f) const
    {
      auto &info = tvT.finfo(f);
      Scalar discontinuities = (info.cumul_tv == 0)
                                   ? 1.0
                                   : info.cumul_tv / info.disc;
      return std::min(10.0,
                      std::max(0.0, 10.0 * discontinuities * info.rel_tv));
      // return ( info.disc <= 1.0 )
      // 	? 10.0
      // 	: std::min( 10.0,
      // 		    std::max( 0.0, discontinuities*info.rel_tv ) ); //(1.0 - info.cont, ) * 10.0 );
    }

    void drawVertex(TVT &tvT, const Vertex vtx)
    {
      const RealPoint rpv = tvT.T.position(vtx);
      const Value vv = tvT.u(vtx);
      const Point pv = _dig(rpv);
      // const Scalar       c = crispnessAtVertex( tvT, vtx );
      _arcimage_simi.setValue(pv, 0);
      _pixinfo.setValue(pv, PixelInformation{vv, 0.0});
    }

    void drawFaceSimilarities(TVT &tvT, const Face f)
    {
      const Arc invalid_arc = tvT.T.nbArcs();
      typedef SimpleMatrix<Scalar, 2, 2> Matrix;
      typedef typename Matrix::ColumnVector CVector;
      std::array<Scalar, 3> s;
      const auto arcs = tvT.arcsAroundFace(f);
      // const int            m = tvT.arcDissimilarities( s, arcs );
      const RealPoint b = tvT.barycenter(f);
      for (Dimension mk = 0; mk < 3; ++mk)
      {
        if (tvT.arcDissimilarity(arcs[mk]) != 0.0)
          continue;
        // Draw similarity arc
        const Vertex vtx1 = tvT.T.head(arcs[mk]);
        const Vertex vtx2 = tvT.T.tail(arcs[mk]);
        std::vector<RealPoint> bp = {tvT.T.position(vtx1),
                                     tvT.T.position(vtx2)};
        const Value v1 = tvT.u(vtx1);
        const Value v2 = tvT.u(vtx2);
        // const Scalar           c1 = crispnessAtVertex( tvT, vtx1 );
        // const Scalar           c2 = crispnessAtVertex( tvT, vtx2 );
        // Traces the digital Bezier curve.
        BezierCurve<Space> B(bp);
        std::vector<Point> dp;
        std::vector<RealPoint> rp;
        std::vector<Scalar> dt;
        B.traceDirect(_dig, dp, rp, dt);
        for (int k = 0; k < dp.size(); ++k)
        {
          if (_draw_domain.isInside(dp[k]))
          {
            PixelInformation pi = {
                combine(v1, v2, dt[k]), // value at pixel
                0.0                     // combine( c1, c2, dt[ k ] )
            };
            _pixinfo.setValue(dp[k], pi);
            _arcimage_simi.setValue(dp[k], arcs[mk]);
            //_arcimage_disc.setValue( dp[ k ], invalid_arc );
          }
        }
      }
    }

    void drawFaceDiscontinuities(TVT &tvT, const Face f)
    {
      typedef SimpleMatrix<Scalar, 2, 2> Matrix;
      typedef typename Matrix::ColumnVector CVector;
      std::array<Scalar, 3> s;
      const Arc invalid_arc = tvT.T.nbArcs();
      const auto arcs = tvT.arcsAroundFace(f);
      //      const int            m = tvT.arcDissimilarities( s, arcs );
      const RealPoint x[3] = {tvT.edgeContourPoint(arcs[0]),
                              tvT.edgeContourPoint(arcs[1]),
                              tvT.edgeContourPoint(arcs[2])};
      const RealPoint b = tvT.barycenter(f);
      const Scalar crisp_f = crispness(tvT, f);
      int m = -1;
      for (int i = 0; i < 3; ++i)
      {
        if ((tvT.arcDissimilarity(arcs[i]) == 0.0) && (tvT.arcDissimilarity(arcs[(i + 1) % 3]) > 0.0) && (tvT.arcDissimilarity(arcs[(i + 2) % 3]) > 0.0))
          m = i;
      }
      if (m >= 0)
      {
        // This is a contour face, we have to connect two sides.
        // Computes the control points of the Bezier curve
        const Dimension i1 = (m + 1) % 3;
        const Dimension j1 = (m + 2) % 3;
        const RealVector u = tangentAtArc(tvT, arcs[i1]);
        const RealVector v = tangentAtArc(tvT, arcs[j1]);
        const auto tu = tvT.intersect(x[i1], x[i1] + u,
                                      x[j1], x[j1] + v);
        RealPoint bb = 0.55 * b + 0.45 * (x[i1] + tu.first * u);
        if (!tvT.isInTriangle(f, bb))
          bb = b;
        Matrix M = {u[0], v[0], u[1], v[1]};
        CVector R = {(8.0 * bb[0] - 4.0 * (x[i1][0] + x[j1][0])) / 3.0,
                     (8.0 * bb[1] - 4.0 * (x[i1][1] + x[j1][1])) / 3.0};
        CVector S = (fabs(M.determinant()) > 1e-4)
                        ? (M.inverse() * R)
                        : CVector{0, 0};
        std::vector<RealPoint> bp = {x[i1],
                                     (S[0] != 0.0) ? (x[i1] + S[0] * u) : (0.66 * x[i1] + 0.34 * x[j1]),
                                     (S[1] != 0.0) ? (x[j1] + S[1] * v) : (0.34 * x[i1] + 0.66 * x[j1]),
                                     x[j1]};
        // Traces the digital Bezier curve.
        BezierCurve<Space> B(bp);
        std::vector<Point> dp;
        std::vector<RealPoint> rp;
        std::vector<Scalar> dt;
        B.traceDirect(_dig, dp, rp, dt);
        const Value vi = valueAtArc(tvT, arcs[i1]);
        const Value vj = valueAtArc(tvT, arcs[j1]);
        const Scalar ci = crispnessAtArc(tvT, arcs[i1]);
        const Scalar cj = crispnessAtArc(tvT, arcs[j1]);
        for (int k = dp.size() / 2; k >= 0; --k)
        {
          if (!_draw_domain.isInside(dp[k]))
            continue;
          if (!tvT.isApproximatelyInTriangle(f, rp[k], 0.05))
            continue;
          // if ( ! tvT.isInTriangle( f, rp[ k ] ) )         continue;
          _arcimage_disc.setValue(dp[k], (dt[k] < 0.5) ? arcs[i1] : arcs[j1]);
          //if ( _arcimage_simi( dp[ k ] ) != invalid_arc ) continue; //break;
          PixelInformation pi = {
              combine(vi, vj, dt[k]), // value at pixel
              combine(ci, cj, dt[k])  // crisp_f
          };
          _pixinfo_disc.setValue(dp[k], pi);
        }
        for (int k = dp.size() / 2 + 1; k < dp.size(); ++k)
        {
          if (!_draw_domain.isInside(dp[k]))
            continue;
          if (!tvT.isApproximatelyInTriangle(f, rp[k], 0.05))
            continue;
          // if ( ! tvT.isInTriangle( f, rp[ k ] ) )         continue;
          _arcimage_disc.setValue(dp[k], (dt[k] < 0.5) ? arcs[i1] : arcs[j1]);
          //if ( _arcimage_simi( dp[ k ] ) != invalid_arc ) continue; //break;
          PixelInformation pi = {
              combine(vi, vj, dt[k]), // value at pixel
              combine(ci, cj, dt[k])  // crisp_f
          };
          _pixinfo_disc.setValue(dp[k], pi);
        }
      }
      else
      {
        // This is a generic face, we connect the three sides to the barycenter.
        // We draw the barycenter in all cases.
        const Point db = _dig(b);
        Value vj = valueAtBarycenter(tvT, f);
        if (_draw_domain.isInside(db))
        {
          _arcimage_disc.setValue(db, 0);
          //if ( _arcimage_simi( db ) == invalid_arc ) {
          PixelInformation pi = {vj, crisp_f};
          _pixinfo_disc.setValue(db, pi);
          //}
        }
        for (Dimension i1 = 0; i1 < 3; ++i1)
        {
          if (tvT.arcDissimilarity(arcs[i1]) == 0.0)
            continue;
          // Computes the control points of the Bezier curve
          const RealVector u = tangentAtArc(tvT, arcs[i1]);
          const Scalar a = 2.0 * (x[i1] - b).dot(u);
          RealPoint bb = (fabs(a) < 1e-4)
                             ? (0.5 * (x[i1] + b))
                             : (x[i1] - ((x[i1] - b).dot(x[i1] - b) / a) * u);
          if (!tvT.isInTriangle(f, bb))
            bb = (0.5 * (x[i1] + b));
          std::vector<RealPoint> bp = {b, bb, x[i1]};
          // Traces the digital Bezier curve.
          BezierCurve<Space> B(bp);
          std::vector<Point> dp;
          std::vector<RealPoint> rp;
          std::vector<Scalar> dt;
          B.traceDirect(_dig, dp, rp, dt);
          Value vi = valueAtArc(tvT, arcs[i1]);
          Scalar crisp_i1 = crispnessAtArc(tvT, arcs[i1]);
          for (int k = 0; k < dp.size(); ++k)
          {
            if (!_draw_domain.isInside(dp[k]))
              continue;
            if (!tvT.isApproximatelyInTriangle(f, rp[k], 0.05))
              continue;
            _arcimage_disc.setValue(dp[k], arcs[i1]);
            //if ( _arcimage_simi( dp[ k ] ) != invalid_arc ) continue; //break;
            PixelInformation pi = {
                combine(vj, vi, dt[k]), // value at pixel
                combine(crisp_f, crisp_i1, dt[k])};
            _pixinfo_disc.setValue(dp[k], pi);
          }
        }
      }
    }

    /// This class is a point functor that returns 'true' only if
    /// their is an invalid arc at the given point, i.e. it is outside
    /// a discontinuity or similarity line of the image.
    struct OutsideLines
    {
      typedef Self::Point Point;
      const ArcImage *_arcimage;
      Arc _invalid_arc;
      OutsideLines(const ArcImage &arcimage, Arc inv)
          : _arcimage(&arcimage), _invalid_arc(inv) {}
      bool operator()(Point p) const
      {
        return _arcimage->operator()(p) == _invalid_arc;
      }
    };

    /**
       Displays the TV triangulation with discontinuities and
       second-order reconstruction.
    */
    void view2ndOrder(TVT &tvT, Scalar discontinuities,
                      bool display_lines = false)
    {
      cairo_set_operator(_cr, CAIRO_OPERATOR_OVER);
      cairo_set_line_width(_cr, 0.0);
      cairo_set_line_cap(_cr, CAIRO_LINE_CAP_BUTT);
      cairo_set_line_join(_cr, CAIRO_LINE_JOIN_BEVEL);
      // Pre-compute soem information per faces (gradient, etc)
      trace.info() << "Compute 2nd order information per face" << std::endl;
      tvT.compute2ndOrderInformation(discontinuities);
      // Resets the image that will be used to compute distance to contours.
      trace.info() << "Invalidate images" << std::endl;
      const Arc invalid_arc = tvT.T.nbArcs();
      for (auto p : _draw_domain)
      {
        _arcimage_disc.setValue(p, invalid_arc);
        _arcimage_simi.setValue(p, invalid_arc);
      }
      trace.info() << "Draw similarity contours" << std::endl;
      // Compute and trace the contour lines of the image.
      for (Face f = 0; f < tvT.T.nbFaces(); ++f)
        drawFaceSimilarities(tvT, f);
      trace.info() << "Draw discontinuity contours" << std::endl;
      for (Face f = 0; f < tvT.T.nbFaces(); ++f)
        drawFaceDiscontinuities(tvT, f);
      trace.info() << "Draw vertices" << std::endl;
      for (Vertex v = 0; v < tvT.T.nbVertices(); ++v)
        drawVertex(tvT, v);
      // Compute Voronoi map and distance transformation.
      trace.info() << "Compute Voronoi diagram of discontinuity lines" << std::endl;
      OutsideLines out_disclines(_arcimage_disc, invalid_arc);
      OutsideLines out_similines(_arcimage_simi, invalid_arc);
      BreadthFirstVisitorWithParent<Space> voronoimap_disc(_draw_domain,
                                                           [&](Point p)
                                                           { return _arcimage_disc(p) != invalid_arc && _arcimage_simi(p) == invalid_arc; });
      while (!voronoimap_disc.finished())
        voronoimap_disc.expand(out_similines);
      voronoimap_disc.computeRoots();
      trace.info() << "Compute Voronoi diagram of similarity lines" << std::endl;
      BreadthFirstVisitorWithParent<Space> voronoimap_simi(_draw_domain,
                                                           [&](Point p)
                                                           { return _arcimage_simi(p) != invalid_arc && _arcimage_disc(p) == invalid_arc; });
      while (!voronoimap_simi.finished())
        voronoimap_simi.expand(out_disclines);
      voronoimap_simi.computeRoots();
      trace.info() << "Interpolate value between lines" << std::endl;
      for (auto p : _draw_domain)
      {
        if ((_arcimage_disc(p) != invalid_arc) && (_arcimage_simi(p) != invalid_arc))
        {
          drawPixelInOutput(p, _pixinfo(p).v);
          continue;
        }
        const auto cp_disc = voronoimap_disc(p); // get closest contour points
        const auto cp_simi = voronoimap_simi(p); // get closest contour points
        drawPixelInOutput(p, evaluateValue(p, cp_disc, cp_simi));
      }
      if (!display_lines)
      {
        trace.info() << "Antialias discontinuity pixels" << std::endl;
        antialiasDiscontinuities(tvT);
        antialiasDiscontinuities2ndPass(tvT);
      }
      // To display contour lines
      if (display_lines)
      {
        trace.info() << "Debug similarity and discontinuity lines" << std::endl;
        for (auto p : _draw_domain)
        {
          if ((_arcimage_disc(p) != invalid_arc) && (_arcimage_simi(p) != invalid_arc))
            drawPixelInOutput(p, Value(180, 0, 180));
          else if (_arcimage_disc(p) != invalid_arc)
            drawPixelInOutput(p, Value(255, 0, 0));
          else if (_arcimage_simi(p) != invalid_arc)
            drawPixelInOutput(p, Value(0, 0, 255));
        }
      }
      trace.info() << "Draw image in cairo buffer" << std::endl;
      drawImage();
    }

    void antialiasDiscontinuities(TVT &tvT)
    {
      const Arc invalid_arc = tvT.T.nbArcs();
      std::array<Vector, 8> deltas =
          {Vector(1, 0), Vector(1, 1), Vector(0, 1),
           Vector(-1, 1), Vector(-1, 0), Vector(-1, -1),
           Vector(0, -1), Vector(1, -1)};
      for (auto p : _draw_domain)
      {
        if (_arcimage_disc(p) == invalid_arc)
          continue;
        if (_arcimage_simi(p) != invalid_arc)
          continue;
        Value acc = Value::zero;
        Scalar w = 0;
        for (auto v : deltas)
        {
          const Point q = p + v;
          if (_draw_domain.isInside(q))
          {
            Scalar qw = 1.0;
            if (_arcimage_disc(q) != invalid_arc)
              qw /= 4.0;
            if (_arcimage_simi(q) != invalid_arc)
              qw *= 4.0;
            acc += qw * _output(q);
            w += qw;
          }
        }
        _output.setValue(p, acc / w);
      }
    }

    void antialiasDiscontinuities2ndPass(TVT &tvT)
    {
      const Arc invalid_arc = tvT.T.nbArcs();
      std::array<Vector, 9> deltas =
          {Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1),
           Vector(-1, 1), Vector(-1, 0), Vector(-1, -1),
           Vector(0, -1), Vector(1, -1)};
      std::array<Scalar, 9> weights =
          {2.0, 1.0, 0.7, 1.0, 0.7, 1.0, 0.7, 1.0, 0.7};
      for (auto p : _draw_domain)
      {
        if (_arcimage_disc(p) == invalid_arc)
          continue;
        Value acc = Value::zero;
        Scalar w = 0;
        for (Dimension i = 0; i < deltas.size(); i++)
        {
          const Point q = p + deltas[i];
          if (_draw_domain.isInside(q))
          {
            acc += weights[i] * _output(q);
            w += weights[i];
          }
        }
        _output.setValue(p, acc / w);
      }
    }

    void drawPixelInOutput(Point p, Value v)
    {
      _output.setValue(p, v);
    }

    void drawImage()
    {
      for (auto p : _draw_domain)
      {
        drawPixel(p, _output(p));
      }
    }

    Value
    evaluateValue(Point p, Point cp_disc, Point cp_simi)
    {
      const PixelInformation &pi_disc = _pixinfo_disc(cp_disc);
      const PixelInformation &pi_simi = _pixinfo(cp_simi);
      Scalar ldisc = (cp_disc - p).norm();
      Scalar lsimi = (cp_simi - p).norm();
      ldisc *= pi_disc.crisp + 1.0;
      // lsimi        *= pi_simi.crisp + 1.0;
      if (lsimi == 0.0 && ldisc == 0.0)
      {
        lsimi = ldisc = 1.0;
      }
      if (pi_disc.crisp <= 1.0)
      {
        Scalar l = ldisc + lsimi;
        return (lsimi / l) * pi_disc.v + (ldisc / l) * pi_simi.v;
      }
      else
      {
        Scalar l = (ldisc + lsimi) / 2.0;
        Value mid = pi_disc.v + _am * (pi_simi.v - pi_disc.v);
        if (ldisc <= lsimi)
        {
          Scalar t = ldisc / l;
          return (1.0 - t) * pi_disc.v + t * mid;
        }
        else
        {
          Scalar t = lsimi / l;
          return (1.0 - t) * pi_simi.v + t * mid;
        }
      }
    }

    // Value evaluateValue( Point p, Point cp ) {
    //   PixelInformation& pi = _pixinfo[ cp ];
    //   RealVector ru( xy( RealVector( cp[ 0 ], cp[ 1 ] ) )
    // 		     - xy( RealVector( p[ 0 ], p[ 1 ] ) ) );
    //   Value       v = pi.v;
    //   VectorValue g = pi.grad;
    //   Scalar      c = pi.crisp;
    //   Scalar    lru = ru.norm();
    //   ru /= lru;
    //   // if ( c >= 1.0 ) lru = 1.0/(1.0+exp(-c*lru) );
    //   if ( c >= 1.0 ) lru = 2.0/(1.0+exp(-c*lru) ) - 1.0;
    //   ru *= lru;
    //   for ( Dimension m = 0; m < 3; ++m )
    // 	v[ m ] += 2.0 * ( ru[ 0 ]*g.x[ m ] + ru[ 1 ]*g.y[ m ] );
    //   return v.inf( pi.max_v ).sup( pi.min_v );
    // }
  };

  // shading; 0:flat, 1:gouraud, 2:linear gradient, 3: partition of unity
  // 4 : 2nd order reconstruction.
  void viewTVTriangulation(TVTriangulation &tvT, double b, double x0, double y0, double x1, double y1,
                           int shading, bool color, std::string fname, double discontinuities,
                           double stiffness, double amplitude)
  {
    CairoViewerTV cviewer(x0, y0, x1, y1,
                          b, b, shading, color, stiffness, amplitude);
    // CairoViewerTV cviewer
    //   ( (int) round( x0 ), (int) round( y0 ),
    // 	(int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ),
    // 	b, b, shading, color, stiffness, amplitude );
    if (shading < 4)
      cviewer.view(tvT, discontinuities);
    else if (shading < 6)
      cviewer.view2ndOrder(tvT, discontinuities,
                           shading == 5);
    else
      cviewer.viewGradientLaplacian(tvT);
    //      cviewer.viewLaplacian( tvT );
    cviewer.save(fname.c_str());
  }

  void viewTVTriangulationAll(TVTriangulation &tvT, double b, double x0, double y0, double x1, double y1,
                              bool color, std::string fname, int display, double discontinuities,
                              double stiffness, double amplitude)
  {
    if (display & 0x1)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 0, color, fname + ".png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x2)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 1, color, fname + "-g.png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x4)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 2, color, fname + "-lg.png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x8)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 3, color, fname + "-pu.png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x10)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 4, color, fname + "-2nd.png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x20)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 5, color, fname + "-2nd-with-lines.png",
                          discontinuities, stiffness, amplitude);
    if (display & 0x40)
      viewTVTriangulation(tvT, b, x0, y0, x1, y1, 6, color, fname + "-laplacian.png",
                          discontinuities, stiffness, amplitude);
  }

  void exportVectMesh(TVTriangulation &tvT, const std::string &name, unsigned int width,
                      unsigned int height, bool displayMesh = true, double scale = 1.0)
  {
    BasicVectoImageExporter exp(name, width, height, displayMesh, 100);
    BasicVectoImageExporter expMean("mean.eps", width, height, displayMesh, scale);

    for (TVTriangulation::Face f = 0; f < tvT.T.nbFaces(); f++)
    {
      TVTriangulation::VertexRange V = tvT.T.verticesAroundFace(f);
      std::vector<std::pair<TVTriangulation::RealPoint, bool>> tr;
      tr.push_back(std::make_pair(tvT.T.position(V[0]), false));
      tr.push_back(std::make_pair(tvT.T.position(V[1]), false));
      tr.push_back(std::make_pair(tvT.T.position(V[2]), false));
      tr.push_back(std::make_pair(tvT.T.position(V[0]), false));
      TVTriangulation::Value valMedian;
      TVTriangulation::Value valMean = tvT.u(V[0]) + tvT.u(V[1]) + tvT.u(V[2]);
      valMean /= 3.0;
      double val1 = sqrt(tvT.u(V[0])[0] * tvT.u(V[0])[0] +
                         tvT.u(V[0])[1] * tvT.u(V[0])[1] +
                         tvT.u(V[0])[2] * tvT.u(V[0])[2]);
      double val2 = sqrt(tvT.u(V[1])[0] * tvT.u(V[1])[0] +
                         tvT.u(V[1])[1] * tvT.u(V[1])[1] +
                         tvT.u(V[1])[2] * tvT.u(V[1])[2]);
      double val3 = sqrt(tvT.u(V[2])[0] * tvT.u(V[2])[0] +
                         tvT.u(V[2])[1] * tvT.u(V[2])[1] +
                         tvT.u(V[2])[2] * tvT.u(V[2])[2]);
      if ((val1 >= val2 && val1 <= val3) ||
          (val1 <= val2 && val1 >= val3))
      {
        valMedian = tvT.u(V[0]);
      }
      else if ((val2 <= val1 && val2 >= val3) ||
               (val2 >= val1 && val2 <= val3))
      {
        valMedian = tvT.u(V[1]);
      }
      else
      {
        valMedian = tvT.u(V[2]);
      }

      exp.addRegion(tr, DGtal::Color(valMedian[0], valMedian[1], valMedian[2]), 0.001);
      expMean.addRegion(tr, DGtal::Color(valMean[0], valMean[1], valMean[2]), 0.001);
    }

    exp.closeFigure();
  }

  TVTriangulation::Arc pivotNext(TVTriangulation &tvT, TVTriangulation::Arc a,
                                 const TVTriangulation::Value &valTrack)
  {
    TVTriangulation::Value currentHead = tvT.u(tvT.T.head(a));
    while (currentHead[0] == valTrack[0] && currentHead[1] == valTrack[1] && currentHead[2] == valTrack[2])
    {
      a = tvT.T.next(a);
      currentHead = tvT.u(tvT.T.head(a));
    }
    return tvT.T.opposite(a);
  }

  std::vector<std::pair<TVTriangulation::RealPoint, bool>>
  trackBorderFromFace(TVTriangulation &tvT,
                      TVTriangulation::Face startArc,
                      TVTriangulation::Value valInside,
                      std::vector<bool> &markedArcs) // TODO
  {
    std::vector<std::pair<TVTriangulation::RealPoint, bool>> res;

    // starting ext point: arc tail
    TVTriangulation::Face faceIni = tvT.T.faceAroundArc(startArc);

    TVTriangulation::Arc currentArc = startArc;
    TVTriangulation::Face currentFace = faceIni;
    markedArcs[startArc] = true;

    do
    {
      TVTriangulation::VertexRange V = tvT.T.verticesAroundFace(currentFace);
      TVTriangulation::RealPoint center = tvT.barycenter(currentFace); //(tvT.T.position(V[0])+tvT.T.position(V[1])+tvT.T.position(V[2]))/3.0;

      // Get the colors around the face
      TVTriangulation::Value a = tvT.u(V[0]);
      TVTriangulation::Value b = tvT.u(V[1]);
      TVTriangulation::Value c = tvT.u(V[2]);
      // The point is fixed (is line) if the 3 colors are differents
      res.push_back(std::make_pair(center, a != b && b != c && a != c));

      currentArc = pivotNext(tvT, currentArc, valInside);
      currentFace = tvT.T.faceAroundArc(currentArc);
      if (currentFace == TVTriangulation::Triangulation::INVALID_FACE)
      {
        break;
      }
      markedArcs[currentArc] = true;
    } while (currentFace != faceIni);
    return res;
  }

  std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>>
  trackBorders(TVTriangulation &tvT, unsigned int num)
  {
    typedef std::map<DGtal::Color, std::vector<unsigned int>> MapColorContours;
    MapColorContours mapContours;
    std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>> resAll;
    std::vector<bool> markedArcs(tvT.T.nbArcs());
    for (unsigned int i = 0; i < markedArcs.size(); i++)
    {
      markedArcs[i] = false;
    }
    bool found = true;
    while (found)
    {
      found = false;
      for (unsigned int a = 0; a < markedArcs.size(); a++)
      {
        // tracking Head color
        TVTriangulation::Value valH = tvT.u(tvT.T.head(a));
        TVTriangulation::Value valT = tvT.u(tvT.T.tail(a));

        found = !markedArcs[a] && (valH[0] != valT[0] || valH[1] != valT[1] || valH[2] != valT[2]) && !tvT.T.isArcBoundary(a);
        if (found)
        {
          resAll.push_back(trackBorderFromFace(tvT, a, valH, markedArcs));
          if (mapContours.count(DGtal::Color(valH[0], valH[1], valH[2])) == 0)
          {
            std::vector<unsigned int> indexC;
            indexC.push_back(resAll.size() - 1);
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])] = indexC;
          }
          else
          {
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])].push_back(resAll.size() - 1);
          }
        }
      }
    }
    auto itMap = mapContours.begin();
    for (unsigned int i = 0; i < num; i++)
      itMap++;
    std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>> res;
    for (unsigned int i = 0; i < (itMap->second).size(); i++)
    {
      res.push_back(resAll[(itMap->second)[i]]);
    }
    return res;
  }

  /**
   * Invalidate contour border and compute the median image color.
   * 
   **/
  DGtal::Color invalidateImageBorder(TVTriangulation &tvT)
  {
    std::vector<TVTriangulation::Value> vectColors;
    for (TVTriangulation::Arc a = 0; a < tvT.T.nbArcs(); a++)
    {
      if (tvT.T.isArcBoundary(a))
      {
        vectColors.push_back(tvT.u(tvT.T.head(a)));
      }
    }
    for (TVTriangulation::Arc a = 0; a < tvT.T.nbArcs(); a++)
    {
      if (tvT.T.isArcBoundary(a))
      {

        tvT.invalidate(tvT.T.head(a));
        tvT.invalidate(tvT.T.tail(a));
      }
    }
    std::sort(vectColors.begin(), vectColors.end(), [](const TVTriangulation::Value &a, const TVTriangulation::Value &b)
              { return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] > b[0] * b[0] + b[1] * b[1] + b[2] * b[2]; });
    TVTriangulation::Value valMed = vectColors[vectColors.size() / 2];
    return DGtal::Color(valMed[0], valMed[1], valMed[2]);
  }

  std::vector<TVTriangulation::ColorContours> trackAllBorders(TVTriangulation &tvT, unsigned int width, unsigned int height)
  {
    typedef std::map<DGtal::Color, std::vector<unsigned int>> MapColorContours;
    std::vector<TVTriangulation::ColorContours> res;
    MapColorContours mapContours;
    std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>> resAll;
    std::vector<bool> markedArcs(tvT.T.nbArcs());
    for (unsigned int i = 0; i < markedArcs.size(); i++)
    {
      markedArcs[i] = false;
    }
    DGtal::Color med = invalidateImageBorder(tvT);

    // Adding background:
    TVTriangulation::ColorContours c;
    c.first = med;
    std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>>
        bg = {{std::make_pair(TVTriangulation::RealPoint(0, height), true),
               std::make_pair(TVTriangulation::RealPoint(width, height), true),
               std::make_pair(TVTriangulation::RealPoint(width, 0), true),
               std::make_pair(TVTriangulation::RealPoint(0, 0), true)}};
    c.second = bg;
    res.push_back(c);

    bool found = true;
    while (found)
    {
      found = false;
      for (unsigned int a = 0; a < markedArcs.size(); a++)
      {
        // tracking Head color
        TVTriangulation::Value valH = tvT.u(tvT.T.head(a));
        TVTriangulation::Value valT = tvT.u(tvT.T.tail(a));

        // Si on est pas passé sur l'arc et que les 2 points sont différents et que ce n'est pas un arc frontière et que le point en tête est valide
        found = !markedArcs[a] && (valH[0] != valT[0] || valH[1] != valT[1] || valH[2] != valT[2]) && !tvT.T.isArcBoundary(a) && !tvT.isinvalid(tvT.T.head(a));
        if (found)
        {
          resAll.push_back(trackBorderFromFace(tvT, a, valH, markedArcs));
          if (mapContours.count(DGtal::Color(valH[0], valH[1], valH[2])) == 0)
          {
            std::vector<unsigned int> indexC;
            indexC.push_back(resAll.size() - 1);
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])] = indexC;
          }
          else
          {
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])].push_back(resAll.size() - 1);
          }
        }
      }
    }

    for (auto cc : mapContours)
    {
      TVTriangulation::ColorContours c;
      c.first = cc.first;
      for (unsigned int i = 0; i < (cc.second).size(); i++)
      {
        c.second.push_back(resAll[(cc.second)[i]]);
      }
      res.push_back(c);
    }
    return res;
  }

  void exportVectMeshDual(TVTriangulation &tvT, const std::string &name, unsigned int width,
                          unsigned int height, bool displayMesh, unsigned int numColor, double scale)
  {
    BasicVectoImageExporter exp(name, width, height, displayMesh, scale);
    for (TVTriangulation::VertexIndex v = 0; v < tvT.T.nbVertices(); v++)
    {
      std::vector<std::pair<TVTriangulation::RealPoint, bool>> tr;
      auto outArcs = tvT.T.outArcs(v);
      for (auto rit = outArcs.rbegin(), ritEnd = outArcs.rend();
           rit != ritEnd; ++rit)
      {
        auto a = *rit;
        tr.push_back(std::make_pair(tvT.contourPoint(a), false));
        if (tvT.T.isArcBoundary(a))
        {
          trace.warning() << "Boundary arc" << std::endl;
        }
        else
        {
          auto f = tvT.T.faceAroundArc(a);
          tr.push_back(std::make_pair(tvT.barycenter(f), false));
        }
      }
      TVTriangulation::Value val = tvT.u(v);
      exp.addRegion(tr, DGtal::Color(val[0], val[1], val[2]), 0.001);
    }
    if (displayMesh)
    {
      for (TVTriangulation::Face f = 0; f < tvT.T.nbFaces(); f++)
      {
        TVTriangulation::VertexRange V = tvT.T.verticesAroundFace(f);
        std::vector<std::pair<TVTriangulation::RealPoint, bool>> tr;
        tr.push_back(std::make_pair(tvT.T.position(V[0]), false));
        tr.push_back(std::make_pair(tvT.T.position(V[1]), false));
        tr.push_back(std::make_pair(tvT.T.position(V[2]), false));
        tr.push_back(std::make_pair(tvT.T.position(V[0]), false));

        exp.addContour(tr, DGtal::Color(0, 200, 200), 0.01);
      }
      std::vector<std::vector<std::pair<TVTriangulation::RealPoint, bool>>> contour = trackBorders(tvT, numColor);
      for (auto c : contour)
      {
        DGtal::Color col;
        /*if (ContourHelper::isCounterClockWise(c))
        {
          col = DGtal::Color(200, 20, 200);
        }
        else
        {*/
        col = DGtal::Color(200, 200, 20);
        //}
        exp.addContour(c, col, 0.1);
      }
    }
    exp.closeFigure();
  }

  void exportVectContoursDual(TVTriangulation &tvT, const std::vector<std::string> &vectName, unsigned int width,
                              unsigned int height, double scale)
  {
    std::vector<TVTriangulation::ColorContours> contourCol = trackAllBorders(tvT, width, height);
    for (auto name : vectName)
    {
      BasicVectoImageExporter exp(name, width, height, false, scale);
      for (auto c : contourCol)
      {
        DGtal::Color col = c.first;
        exp.addRegions(c.second, col);
      }
      exp.closeFigure();
    }
  }

} // namespace DGtal


int main(int argc, char **argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
//  std::string inputFileName;
//  std::string outputFileName {"result.raw"};
//  DGtal::int64_t rescaleInputMin {0};
//  DGtal::int64_t rescaleInputMax {255};
//
//
  double dt {0.248};
  double lambda {0.0};
  double tol {0.01};
  int N  {20}; //= vm["tv-max-iter"].as<int>();
  int quant  {256}; //= vm["quantify"].as<int>();
  std::string img_fname;
  std::string bname {"after-tv-opt"};

  double p {0.5};
  double sim {0.0};
  std::string conn {"Order"};
  bool debug {false};

  
  int display {0};//   = vm["display-tv"].as<int>();
  double b {2.0};
  double st  {0.9}; //= vm["stiffness"].as<double>();

  
  double disc {0.1};
  double am {1.0};
  int miter {200};
  int strat {4};
  int nbAlt {1};
  double vectScale;
  
  int displayFlip {4};
  std::vector<std::string> vectMeshName;
  std::vector<std::string> vectMeshDualName;
  std::vector<std::string> vectContoursDualName;
  int numColorExportVectDual {0};
  int nbItRegContour {20};
  bool versionDispl {false};
  bool displayMesh {false};
  double zip {1.0};
  std::string z_method {"Laplacian"};

  using namespace DGtal;
  app.description("Generate the best piecewise linear TV from a color image. It works by simply flipping edges of"
                  "some initial triangulation to optimize the given energy."
                  "Basic usage: \n \ttv-triangulation-color [options] -i <image.ppm> -b 4");
  app.add_option("-i,--input,1", img_fname, "Specifies the input shape as a 2D image PPM filename.")
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o", bname, "Specifies the basename of output displayed images. See option -D, file name is completed by -g, -lg, -pu, -2nd, -2nd-with-slices depending on display choice.", true);

  app.add_option("--lambda,-l", lambda, "The data fidelity term in TV denoising (if lambda <= 0, then the data fidelity is exact.", true);
  app.add_option("--dt", dt, "The time step in TV denoising (should be lower than 0.25).", true);
  app.add_option("--tolerance,-t", tol, "The tolerance to stop the TV denoising.", true);
  app.add_option("--tv-max-iter,-N", N, "The maximum number of iteration in TV's algorithm.", true);
  app.add_option("--quantify", quant, "The quantification for colors (number of levels, q=2 means binary.", true);

  app.add_option("--tv-power,-p", p, "The power coefficient used to compute the gradient ie |Grad I|^{2p}. ", true);
  app.add_option("--similarity", sim, "Tells when two colors are considered identical for connectedness.", true);
  app.add_option("--connectivity",conn,  "Indicates the strategy for determining the connectivity of ambiguous pixels: Size | Order. Size is best for 1-bit images, Order is best for color images.", true);
  app.add_flag("--debug", debug, "specifies the DEBUG mode: output some intermediate results in cc.ppm and order.ppm.");
  app.add_option("--display-tv,-d",display, "Tells the display mode after standard TV regularization per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles.", true);
  app.add_option("--bitmap,-b", b, "Rasterization magnification factor [arg] for PNG export.", true);
  app.add_option("--stiffness",st, "Tells how to stiff the gradient around discontinuities (amplitude value is changed at stiffness * middle).", true);
  app.add_option("--discontinuities,-S",disc,  "Tells to display a % of the TV discontinuities (the triangles with greatest energy).", true);
  app.add_option("--amplitude", am,"Tells the amplitude of the stiffness for the gradient around discontinuities.", true);
  app.add_option("--limit,-L",miter, "Gives the maximum number of passes (default is 200).", true );
  app.add_option("--strategy,-s", strat, "Strategy for quadrilatera with equal energy: 0: do nothing, 1: subdivide, 2: flip all, 3: flip all when #flipped normal = 0, 4: flip approximately half when #flipped normal = 0, 5: flip approximately half when #flipped normal = 0, subdivide if everything fails.", true);
  app.add_option("--nb-alt-iter,-A",nbAlt, "The number of iteration for alternating TV and TV-flip.", true);
  auto vectScaleOpt = app.add_option("vectScale", vectScale, "Change the default vectorial scale to ajust default display size. By default the scale is given from the Rasterization magnification factor (--bitmap,b option) (using 10 will display easely small images while 1.0 is more adapted to bigger images)." );
  app.add_option("--display-flip,-D",displayFlip,"Tells the display mode aftergeometric TV flips per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles (quite nice and flat), 0x8: output partition of unity Bezier triangles (slow and not good), 0x10: output 2nd order reconstruction with Bezier curve (nicest but slower than 0x4), 0x20: same as 0x8 but displays also the similarity/discontinuity graph (for debug/illustrations), 0x40: use Poisson reconstruction to reconstruct linear gradient with discontinuities.", true);
  app.add_option("--exportVectMesh,-e", vectMeshName,"Export the triangle mesh.", true );
  app.add_option("--exportVectMeshDual,-E", vectMeshDualName,"Export the triangle mesh.", true );
  app.add_option("--exportVectContoursDual,-C", vectContoursDualName,"Export the image regions filled. Use 3rd bezier curves for the edges.", true );
  app.add_option("--numColorExportVectDual", numColorExportVectDual, "num of the color of the map.", true );

  app.add_flag("--displayMesh", displayMesh, "display mesh of the vectorial display.");
  app.add_option("--regularizeContour,-R", nbItRegContour,"regularizes the dual contours for <nb> iterations.", true );
  app.add_option("--zip,-z", zip, "Compresses the triangulation to keep only the given proportion of vertices.", true);
  app.add_option("--zip-method,-Z",z_method, "zip method in Laplacian | Merge.", true);
  app.add_flag("--version", versionDispl, "Display the version number.");

  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
  
// //
//  // parse command line ----------------------------------------------
//  po::options_description general_opt("Allowed options are: ");
//  general_opt.add_options()("help,h", "display this message")
//    ("output,o", po::value<std::string>()->default_value("after-tv-opt"), "Specifies the basename of output displayed images. See option -D, file name is completed by -g, -lg, -pu, -2nd, -2nd-with-slices depending on display choice.")
//    ("limit,L", po::value<int>()->default_value(200), "Gives the maximum number of passes (default is 200).")
//    ("bitmap,b", po::value<double>()->default_value(2.0), "Rasterization magnification factor [arg] for PNG export.")
//    ("strategy,s", po::value<int>()->default_value(4), "Strategy for quadrilatera with equal energy: 0: do nothing, 1: subdivide, 2: flip all, 3: flip all when #flipped normal = 0, 4: flip approximately half when #flipped normal = 0, 5: flip approximately half when #flipped normal = 0, subdivide if everything fails.")
//    ("tv-power,p", po::value<double>()->default_value(0.5), "The power coefficient used to compute the gradient ie |Grad I|^{2p}. ")
//    ("lambda,l", po::value<double>()->default_value(0.0), "The data fidelity term in TV denoising (if lambda <= 0, then the data fidelity is exact")
//    ("dt", po::value<double>()->default_value(0.248), "The time step in TV denoising (should be lower than 0.25)")
//    ("tolerance,t", po::value<double>()->default_value(0.01), "The tolerance to stop the TV denoising.")
//    ("quantify,q", po::value<int>()->default_value(256), "The quantification for colors (number of levels, q=2 means binary.")
//    ("tv-max-iter,N", po::value<int>()->default_value(20), "The maximum number of iteration in TV's algorithm.")
//    ("nb-alt-iter,A", po::value<int>()->default_value(1), "The number of iteration for alternating TV and TV-flip.")
//    ("display-tv,d", po::value<int>()->default_value(0), "Tells the display mode after standard TV regularization per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles.")
//    ("display-flip,D", po::value<int>()->default_value(4), "Tells the display mode aftergeometric TV flips per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles (quite nice and flat), 0x8: output partition of unity Bezier triangles (slow and not good), 0x10: output 2nd order reconstruction with Bezier curve (nicest but slower than 0x4), 0x20: same as 0x8 but displays also the similarity/discontinuity graph (for debug/illustrations), 0x40: use Poisson reconstruction to reconstruct linear gradient with discontinuities.")
//    ("discontinuities,S", po::value<double>()->default_value(0.1), "Tells to display a % of the TV discontinuities (the triangles with greatest energy).")
    
//    ("amplitude", po::value<double>()->default_value(1.0), "Tells the amplitude of the stiffness for the gradient around discontinuities.")
//    ("similarity", po::value<double>()->default_value(0.0), "Tells when two colors are considered identical for connectedness.")
//    ("connectivity", po::value<std::string>()->default_value("Order"), "Indicates the strategy for determining the connectivity of ambiguous pixels: Size | Order. Size is best for 1-bit images, Order is best for color images.")
    //("debug", "specifies the DEBUG mode: output some intermediate results in cc.ppm and order.ppm.")
//    ("displayMesh", "display mesh of the vectorial display.")
//    ("exportVectMesh,e", po::value<std::vector<std::string>>()->multitoken(), "Export the triangle mesh.")
//    ("exportVectMeshDual,E", po::value<std::vector<std::string>>()->multitoken(), "Export the triangle mesh.")
//    ("exportVectContoursDual,C", po::value<std::vector<std::string>>()->multitoken(), "Export the image regions filled. Use 3rd bezier curves for the edges")
//    ("vectScale", po::value<double>(), "Change the default vectorial scale to ajust default display size. By default the scale is given from the Rasterization magnification factor (--bitmap,b option) (using 10 will display easely small images while 1.0 is more adapted to bigger images).")
//    ("numColorExportVectDual", po::value<unsigned int>()->default_value(0), "num of the color of the map.")
//    ("regularizeContour,R", po::value<int>()->default_value(20), "regularizes the dual contours for <nb> iterations.")
//    ("zip,z", po::value<double>()->default_value(1.0), "Compresses the triangulation to keep only the given proportion of vertices.")
//    ("zip-method,Z", po::value<std::string>()->default_value("Laplacian"), "zip method in Laplacian | Merge.")("version", "Display the version number.")
      // ("nb-zip-geometry,Z", po::value<int>()->default_value( 100 ), "Maximum number of iterations to optimize the TV geometry of the zipped triangulation." )
      ;
 
  

  if (versionDispl)
  {
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Using " << argv[0] << " version " << version << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
  }
  
 
  
  
  
  
  // Useful types
  using namespace std;
  using namespace DGtal;

  typedef Z2i::Space Space;
  typedef Z2i::Domain Domain;
  typedef ImageSelector<Z2i::Domain, unsigned int>::Type Image;
  typedef ImageSelector<Z2i::Domain, Color>::Type ColorImage;
  Clock c;
  trace.beginBlock("Loading image");

  Image image = GenericReader<Image>::import(img_fname);
  std::string extension = img_fname.substr(img_fname.find_last_of(".") + 1);
  bool color = false;
  if (extension == "ppm")
    color = true;
  trace.info() << "Image <" << img_fname
               << "> size=" << image.extent()[0]
               << "x" << image.extent()[1]
               << " color=" << (color ? "True" : "False") << std::endl;
  trace.info() << std::fixed;
  trace.endBlock();


  double timeTV, timeTriangulation, timeDisplayTriangulation,
      timeOptimisationTr, timeRegulContours, timeExport = 0.0;
  c.startClock();
  if (false)
  {
    trace.beginBlock("Image usual TV regularization");
    bool out_color = color;
    if (color)
    {
      typedef ImageTVRegularization<Space, 3> ColorTV;
      ColorTV tv;
      tv.init(image, ColorTV::Color2ValueFunctor());
      tv.optimize(lambda, dt, tol, N);
      if (out_color)
        tv.outputU(image, ColorTV::Value2ColorFunctor(quant));
      else
        tv.outputU(image, ColorTV::Value2GrayLevelFunctor(quant));
    }
    else
    {
      typedef ImageTVRegularization<Space, 1> GrayLevelTV;
      GrayLevelTV tv;
      tv.init(image, GrayLevelTV::GrayLevel2ValueFunctor());
      tv.optimize(lambda, dt, tol, N);
      if (out_color)
        tv.outputU(image, GrayLevelTV::Value2ColorFunctor(quant));
      else
        tv.outputU(image, GrayLevelTV::Value2GrayLevelFunctor(quant));
    }
    trace.endBlock();
  }
  else if (quant != 256)
  {
    trace.beginBlock("Image quantification");
    bool out_color = color;
    if (color)
    {
      typedef ImageTVRegularization<Space, 3> ColorTV;
      ColorTV tv;
      tv.init(image, ColorTV::Color2ValueFunctor());
      if (out_color)
        tv.outputU(image, ColorTV::Value2ColorFunctor(quant));
      else
        tv.outputU(image, ColorTV::Value2GrayLevelFunctor(quant));
    }
    else
    {
      typedef ImageTVRegularization<Space, 1> GrayLevelTV;
      GrayLevelTV tv;
      tv.init(image, GrayLevelTV::GrayLevel2ValueFunctor());
      if (out_color)
        tv.outputU(image, GrayLevelTV::Value2ColorFunctor(quant));
      else
        tv.outputU(image, GrayLevelTV::Value2GrayLevelFunctor(quant));
    }
    trace.endBlock();
  }
  timeTV = c.stopClock();
  c.startClock();
  trace.beginBlock("Construction of the triangulation");
 
 
  int conn_s = conn == "Size" ? 0 : 1;
  TVTriangulation TVT(image, color, p, sim, conn_s, debug);
  trace.info() << TVT.T << std::endl;
  trace.endBlock();
  // by default we take the scale of the scale factor of the zoom image.
  double scaleVisuVect = b;
  if (vectScaleOpt->count() > 0)
  {
    scaleVisuVect = vectScale;
  }
  trace.beginBlock("Output TV image (possibly quantified)");
  Image J(image.domain());
  bool ok = TVT.outputU(J);
  struct UnsignedInt2Color
  {
    Color operator()(unsigned int val) const { return Color(val); }
  };
  struct UnsignedInt2GrayLevel
  {
    unsigned char operator()(unsigned int val) const
    {
      return static_cast<unsigned char>(val);
    }
  };
  PPMWriter<Image, UnsignedInt2Color>::exportPPM("output-tv.ppm", J);
  trace.endBlock();
  timeTriangulation = c.stopClock();
  c.startClock();
  trace.beginBlock("Displaying triangulation");
  {
   
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = (double)image.domain().upperBound()[0];
    double y1 = (double)image.domain().upperBound()[1];
    viewTVTriangulationAll(TVT, b, x0, y0, x1, y1, color, "after-tv",
                           display, disc, st, am);
  }
  trace.endBlock();
  timeDisplayTriangulation = c.stopClock();
  c.startClock();
  trace.beginBlock("Optimizing the triangulation");
 

  if (quant != 256 && nbAlt != 1)
  {
    nbAlt = 1;
    trace.warning() << "Quantification is not compatible with alternating TV + flips" << std::endl;
  }
  std::pair<int, int> nbs;
  double last_energy = -1;
  for (int n = 0; n < nbAlt; ++n)
  {
    trace.info() << "******************** PASS " << n << " *******************"
                 << endl;
    if (n >= 0 && lambda > 0.0)
    {
      trace.info() << "------------- optimize u ----------------" << std::endl;
      TVTriangulation::ValueForm save_u = TVT._u;
      TVT.tvPass(lambda, dt, tol, N);
      if (last_energy >= 0.0 && TVT.getEnergyTV() >= last_energy && n >= nbAlt / 2)
      {
        TVT._u = save_u;
        TVT.computeEnergyTV();
        trace.info() << "TV( u ) = " << TVT.getEnergyTV()
                     << " (unable to do better)" << std::endl;
        break;
      }
    }
    if (n == 0)
    {
      Image L(image.domain());
      TVT.outputLaplacian(L, 1.0);
      PGMWriter<Image, UnsignedInt2GrayLevel>::exportPGM("output-laplacian-before.pgm", L);
    }
    int iter = 0;
    int last = 1;
    bool subdivide = false;
    trace.info() << "------------- optimize geometry --------------" << std::endl;
    trace.info() << "TV( u ) = " << TVT.getEnergyTV() << std::endl;
    while (true)
    {
      if (iter++ >= miter)
        break;
      double energy = 0.0;
      nbs = TVT.onePass(energy, strat);
      if ((last == 0) && (nbs.first == 0))
      {
        if (subdivide || strat != 5)
          break;
        subdivide = true;
        nbs = TVT.onePass(energy, 1);
      }
      last = nbs.first;
    }
    last_energy = TVT.getEnergyTV();
  }
  trace.endBlock();
  timeOptimisationTr = c.stopClock();
  c.startClock();

  trace.beginBlock("regularizing contours");
  {
    TVT.regularizeContours(0.001, nbItRegContour);
  }
  trace.endBlock();
  timeRegulContours = c.stopClock();
  c.startClock();

  trace.beginBlock("Displaying triangulation");
  {
   
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = (double)image.domain().upperBound()[0];
    double y1 = (double)image.domain().upperBound()[1];

    viewTVTriangulationAll(TVT, b, x0, y0, x1, y1, color, bname,
                           displayFlip, disc, st, am);
    // if ( ( display & 64 ) == 64 ) {
    //   typedef ImageContainerBySTLVector< Domain, Color > ColorImage;
    //   typedef ColorImage::Point Point;
    //   typedef ColorImage::Domain Domain;
    //   ColorImage poisson_image( Domain( Point::zero, Point::zero ) );
    //   image::Zoom Z( image.domain(), (int) round( b ) );
    //   TVT.outputPoissonImage( poisson_image, Z );
    //   struct IdFunctor {
    //     typedef Color Argument;
    //     typedef Color Value;
    //     Value operator()( Argument c ) const { return c; }
    //   };
    //   std::string fpoisson = bname + "-poisson.ppm";
    //   PPMWriter< ColorImage, IdFunctor >::exportPPM( fpoisson, poisson_image );
    // }
  }
  trace.endBlock();
  timeDisplayTriangulation = c.stopClock();

  if (z_method == "Merge")
  {
    trace.beginBlock("Compressing triangulation");

    if (zip < 1.0)
    {
      TVTriangulation TVTzip(image.domain(), TVT._u, TVT.T,
                             color, p);
      TVTzip.simplify(zip);
      trace.info() << TVTzip.T << std::endl;

      trace.info() << "------------- optimize geometry --------------" << std::endl;
      trace.info() << "TV( u ) = " << TVTzip.getEnergyTV() << std::endl;
      int iter = 0;
      int last = 1;
      // int miter = vm[ "nb-zip-geometry" ].as<int>();
      while (true)
      {
        if (iter++ >= miter)
          break;
        double energy = 0.0;
        nbs = TVTzip.onePass(energy, strat);
        if ((last == 0) && (nbs.first == 0))
          break;
        last = nbs.first;
      }

     
      double x0 = 0.0;
      double y0 = 0.0;
      double x1 = (double)image.domain().upperBound()[0];
      double y1 = (double)image.domain().upperBound()[1];
      std::string bname = "zip";
      viewTVTriangulationAll(TVTzip, b, x0, y0, x1, y1, color, bname,
                             displayFlip, disc, st, am);
      unsigned int w = image.extent()[0];
      unsigned int h = image.extent()[1];
      std::string pname = "zip-primal.eps";
      // std::string dname = "zip-dual.eps";

      exportVectMesh(TVTzip, pname, w, h, true, scaleVisuVect);
      //      exportVectMeshDual(TVTzip, dname, w, h , true, 0, scaleVisuVect);
    }
    trace.endBlock();
  } //  else if ( z_method == "Laplacian" ) {
  //   trace.beginBlock("Compressing triangulation");
  //   {
  //     Image L( image.domain() );
  //     TVT.outputLaplacian( L, 1.0 );
  //     PGMWriter<Image,UnsignedInt2GrayLevel>::exportPGM( "output-laplacian-after.pgm", L );
  //     double zip = vm[ "zip" ].as<double>();
  //     if ( zip < 1.0 ) {
  //       ImageTriangulation< Z2i::Space, 3 > IT;
  //       TVTriangulation::ScalarForm lp( TVT.T.nbVertices() );
  //       TVTriangulation::VertexIndex v = 0;
  //       for ( auto val : L ) lp[ v++ ] = val;
  //       IT.init( image.domain(), TVT._I, lp, zip );
  //       TVTriangulation TVTzip( IT._domain, IT._I, IT._T,
  //       			color, p );
  //       trace.info() << TVTzip.T << std::endl;

  //       trace.info() << "------------- optimize geometry --------------" << std::endl;
  //       trace.info() << "TV( u ) = " << TVTzip.getEnergyTV() << std::endl;
  //       int iter = 0;
  //       int last = 1;
  //       int miter = vm[ "limit" ].as<int>();
  //       //int miter = vm[ "nb-zip-geometry" ].as<int>();
  //       while ( true ) {
  //         if ( iter++ >= miter ) break;
  //         double energy = 0.0;
  //         nbs = TVTzip.onePass( energy, strat );
  //         if ( ( last == 0 ) && ( nbs.first == 0 ) ) break;
  //         last = nbs.first;
  //       }

  //       int  display = vm[ "display-flip" ].as<int>();
  //       double     b = vm[ "bitmap" ].as<double>();
  //       double  disc = vm[ "discontinuities" ].as<double>();
  //       double    st = vm[ "stiffness" ].as<double>();
  //       double    am = vm[ "amplitude" ].as<double>();
  //       double    x0 = 0.0;
  //       double    y0 = 0.0;
  //       double    x1 = (double) image.domain().upperBound()[ 0 ];
  //       double    y1 = (double) image.domain().upperBound()[ 1 ];
  //       std::string bname = "zip";
  //       viewTVTriangulationAll( TVTzip, b, x0, y0, x1, y1, color, bname,
  //       			display, disc, st, am );
  //       unsigned int w = image.extent()[ 0 ];
  //       unsigned int h = image.extent()[ 1 ];
  //       std::string pname = "zip-primal.eps";
  //       // std::string dname = "zip-dual.eps";
  //       double vectScale = vm["vectScale"].as<double>();
  //       exportVectMesh(TVTzip, pname, w, h , true, vectScale);
  //       //      exportVectMeshDual(TVTzip, dname, w, h , true, 0, vectScale);
  //     }
  //   }
  //   trace.endBlock();
  // }

  c.startClock();

  trace.beginBlock("Export base triangulation");

  if (vectMeshName.size() != 0)
  {
    unsigned int w = image.extent()[0];
    unsigned int h = image.extent()[1];
    for (auto name : vectMeshName)
    {
      exportVectMesh(TVT, name, w, h, displayMesh, scaleVisuVect);
    }
  }
  if (vectMeshDualName.size() != 0)
  {
    unsigned int w = image.extent()[0];
    unsigned int h = image.extent()[1];
    for (auto name : vectMeshDualName)
    {
      unsigned int numColor = numColorExportVectDual;
      exportVectMeshDual(TVT, name, w, h, displayMesh, numColor, scaleVisuVect);
    }
  }
  if (vectContoursDualName.size() != 0)
  {
    unsigned int w = image.extent()[0];
    unsigned int h = image.extent()[1];
    exportVectContoursDual(TVT, vectContoursDualName, w, h, scaleVisuVect);
  }
  trace.endBlock();
  timeExport = c.stopClock();
  std::cout << "----------------------\n";
  std::cout << "Time measure: (ms) \n"
            << "----------------------\n"
            << "TV:" << timeTV << std::endl
            << "Triangulation:" << timeTriangulation << std::endl
            << "Triangulation display:" << timeDisplayTriangulation << std::endl
            << "Optimisation triangulation:" << timeOptimisationTr << std::endl
            << "Contour regularisation: " << timeRegulContours << std::endl
            << "Export vectorial representation:" << timeExport << std::endl
            << "----------------------\n"
            << "Total time: " << timeTV + timeTriangulation + timeDisplayTriangulation + timeOptimisationTr + timeRegulContours + timeExport << std::endl
            << "----------------------\n";

  // trace.beginBlock("Merging triangles");
  // {

  // }
  // trace.endBlock();
  return 0;
}
