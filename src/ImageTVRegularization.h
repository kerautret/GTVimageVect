/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file ImageTVRegularization.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module ImageTVRegularization.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ImageTVRegularization_RECURSES)
#error Recursive header files inclusion detected in ImageTVRegularization.h
#else // defined(ImageTVRegularization_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImageTVRegularization_RECURSES

#if !defined ImageTVRegularization_h
/** Prevents repeated inclusion of headers. */
#define ImageTVRegularization_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/kernel/domains/Linearizer.h>
#include <DGtal/helpers/StdDefs.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class ImageTVRegularization
  /**
     Description of class 'ImageTVRegularization' <p> \brief Aim: This
     class regularizes an input image with respect to its total
     variation. More precisely, for $\Omega$ a domain, the total variation of \f$ u \in L^1(\Omega,\mathbb{R})\f$  is defined as

     \f[
     TV( u ) := \sup \left\{ \int_\Omega u \div \xi \text{d}x : \xi \in C_c^1(\Omega, \mathbb{R}), \|\xi\|_\infty \le 1 \right\}.
     \f]

     Note that \f$ TV(u) \f$ is finite iff its distributional
     derivative \f$ Du \f$ is a finite Radon measure, and in this case
     \f$ TV(u)=|Du|(\Omega) \f$. Also, if \f$ u \f$ has a gradient \f$
     \nabla u \in L^1(\Omega,\R^2) \f$, then $TV(u) = \int_\Omega
     |\nabla u | \f$. This expains why \f$ TV(u) \f$ is often
     "defined" as such, by abuse of the notation. The space of
     functions of \a bounded \a variations is thus defined as

     \f[
     BV( \Omega ) := \left\{ u \in L^1(\Omega) : \int_\Omega | \nabla u | < \infty \right\}.
     \f]

     We solve here the following unconstrained minimization problem,
     where $\lamba$ is a parameter that tunes the fit to original data
     \f$ f \f$:

     \f[
     \min_{u \in BV(\Omega)} \int_\Omega |\nabla u | +
     \frac{\lambda}{2}\| u - f \|_2^2.
     \f]

     Assuming \f$ f \f$ is an input grey-level image, and \f$ g \f$ is
     the output TV denoised image, the following code does the job:

     \code
     template <typename Image>
     void TV( Image& g, const Image& f, double lambda ) {
       typedef typename Image::Domain::Space Space;
       typedef ImageTVRegularization<Space, 1> GrayLevelTV;
       GrayLevelTV tv;
       tv.init( f, GrayLevelTV::GrayLevel2ValueFunctor() );
       tv.optimize( lambda );
       tv.outputU( g, GrayLevelTV::Value2GrayLevelFunctor() );
     }
     \endcode

     @see tv-image.cpp
     
     @tparam TSpace the digital space for images (choose Z2i::Space or
     Z3i::Space according to the dimension of images you are
     processing).
     
     @tparam M the number of scalar per pixel/voxel, generally 1 for
     gray-level images and 3 for color images, but other choices are
     possible.
  */
  template <typename TSpace, int M>
  class ImageTVRegularization
  {
  public:
    typedef TSpace                            Space;
    typedef ImageTVRegularization< TSpace, M> Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension <= 3 ));
    BOOST_STATIC_ASSERT (( M >= 1 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef HyperRectDomain<Space>                Domain;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< M, Scalar >              Value;
    typedef std::array< Value, Space::dimension > VectorValue;
    ///static constants to store the dimension.
    static const Dimension N = Space::dimension;
    typedef std::vector<Scalar>                   ScalarForm;
    typedef std::vector<Value>                    ValueForm;
    typedef std::vector<VectorValue>              VectorValueForm;

    /// The image values at each vertex
    ValueForm            _I;
    /// The domain of the image and of the computations.
    Domain               _domain;
    /// The extent of the domain of the image and of the computations.
    Vector               _extent;
    /// The regularized values at each vertex
    ValueForm            _u;
    /// The TV-regularized vectors
    VectorValueForm      _p;
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /**
     * Destructor.
     */
    ~ImageTVRegularization() {}

    /// Default constructor. The object is invalid.
    ImageTVRegularization() : _domain( Point(), Point() ) {}

    /// Functor used to feed the TV with a color image (M should be 3).
    struct Color2ValueFunctor {
      Scalar operator()( unsigned int color, unsigned int m ) const {
	switch( m ) {
	case 0: return (color >> 16) & 0xff;
	case 1: return (color >> 8) & 0xff;
	default: return color & 0xff;
	}
      }
    };
    /// Functor used to feed the TV with a gray-level image (M should be 1).
    struct GrayLevel2ValueFunctor {
      Scalar operator()( unsigned int color, unsigned int m ) const {
	return color & 0xff;
      }
    };
    
    /// @tparam Functor the type of function: Image Value x int --> Scalar,
    /// where int is the dimension in the Value.
    /// @see Color2ValueFunctor
    /// @see GrayLevel2ValueFunctor
    template <typename Image, typename Functor>
    void init( const Image& I, Functor f )
    {
      _I.reserve( I.size() );
      for ( auto val_I : I ) {
	Value v;
	for ( unsigned int m = 0; m < M; ++m )
	  v[ m ] = f( val_I, m );
	_I.push_back( v );
      }
      _u = _I;                  // u = image at initialization
      _p.resize( _I.size() );   // p = 0     at initialization
      _domain = I.domain();
      _extent = I.extent();
    }

    /// Bounds the scalar value in [0,255] and rounds it to the nearest integer.
    static unsigned int dig( Scalar v )
    {
      return (unsigned int ) round( std::min( 255.0, std::max( 0.0, v ) ) );
    }

    /// Functor for outputing the result into a color image.
    struct Value2ColorFunctor {
      const double factor;
      Value2ColorFunctor( int level = 256 )
	: factor ( 255.0 / ( level - 1 ) ) {}

      double q( double v ) const  {
      	return round( v / factor ) * factor;
      }

      unsigned int operator()( const Value& v ) const
      {
	if ( M <= 2 ) 
	  return ( dig( q( v[ 0 ] ) ) << 16 )
	    +    ( dig( q( v[ 0 ] ) ) << 8 )
	    +    ( dig( q( v[ 0 ] ) ) );
	else 
	  return ( dig( q( v[ 0 ] ) ) << 16 )
	    +    ( dig( q( v[ 1 ] ) ) << 8 )
	    +    ( dig( q( v[ 2 ] ) ) );
      }
    };
    
    /// Functor for outputing the result into a gray-level image.
    struct Value2GrayLevelFunctor {
      const double factor;
      
      Value2GrayLevelFunctor( int level = 256 )
	: factor ( 255.0 / ( level - 1 ) ) {}

      double q( double v ) const  {
      	return round( v / factor ) * factor;
      }
      unsigned int operator()( const Value& v ) const
      {
	unsigned int gl = 0;
	for ( unsigned int m = 0; m < M; ++m )
	  gl += dig( q( v[ m ] ) );
	return gl / M;
      }
    };
    
    /// @see Value2ColorFunctor
    /// @see Value2GrayLevelFunctor
    template <typename Image, typename Functor>
    bool outputU( Image& J, Functor f ) const
    {
      J      = Image( _domain ); 
      Size i = 0;
      for ( auto & val_J : J ) {
	val_J = f( _u[ i ] );
	i++;
      }
      return i == _u.size();
    }

    
    /// Linearize the point to an index.
    Size index( const Point &aPoint ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getIndex( aPoint, _extent );
    }
    /// Delinearize the index to a point.
    Point point( const Size index ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getPoint( index, _extent );
    }

    // ---------------- Calculus norms and operators ------------------
    
    /// Computes the square of x.
    static Scalar square( Scalar x ) { return x*x; }

    /// The norm for Value (Euclidean)
    Scalar normX( const Value& v ) const
    {
      Scalar x = square( v[ 0 ] );
      for ( unsigned int m = 1; m < M; ++m )
	x += square( v[ m ] );
      return sqrt( x );
    }

    /// The norm for VectorValue (TV and VTV)
    Scalar normY( const VectorValue& vv ) const
    {
      Scalar xx = 0.0;
      for ( unsigned int n = 0; n < N; ++n )
	for ( unsigned int m = 0; m < M; ++m )
	  xx += square( vv[ n ][ m ] );
      return sqrt( xx );
    }

    /// @return the vector of the norms of each vector in p (one per vertex). 
    ScalarForm norm( const VectorValueForm& p ) const
    {
      ScalarForm S( p.size() );
      for ( Size i = 0; i < p.size(); i++ )
	S[ i ] = normY( p[ i ] );
      return S;
    }

    // Definition of a global gradient operator that assigns vectors to triangles
    VectorValueForm grad( const ValueForm& u ) const
    { // it suffices to traverse all (valid) triangles.
      VectorValueForm G( u.size() );
      Size i = 0;
      for ( auto p : _domain ) {
	G[ i ] = grad( u, p, i );
	i++;
      }
      return G;
    }

    inline
    VectorValue grad( const ValueForm& u, const Point& p ) {
      return grad( u, p, index( p ) );
    } 

    inline
    VectorValue grad( const ValueForm& u, Size i ) {
      return grad( u, point( i ), i );
    } 

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad( const ValueForm& u, const Point& p,
		      const Size i ) const
    {
      VectorValue G;
      for ( unsigned int n = 0; n < N; ++n )
	if ( p[ n ] < _extent[ n ] - 1 )
	  for ( unsigned int m = 0; m < M; ++m )
	    G[ n ][ m ] = u[ index( p + Point::base( n ) ) ][ m ] - u[ i ][ m ];
      // otherwise zero.
      return G;
    }

    // Definition of a (global) divergence operator that assigns
    // scalars to vertices from a vector field.
    ValueForm div( const VectorValueForm& G ) const
    {
      ValueForm S( G.size() );
      for ( Size i = 0; i < G.size(); ++i )
	S[ i ] = div( G, point( i ), i );
      return S;
    }

    Value div( const VectorValueForm& G, const Point& p ) const {
      return div( G, p, index( p ) );
    }
    Value div( const VectorValueForm& G, Size i ) const {
      return div( G, point( i ), i );
    }

    // Definition of a (local) divergence operator that assigns
    // scalars to vertices from a vector field.
    Value div( const VectorValueForm& G, const Point& p,
	       const Size i ) const
    {
      Value v;
      for ( unsigned int n = 0; n < N; ++n )
	if ( p[ n ] < _extent[ n ] - 1 ) {
	  for ( unsigned int m = 0; m < M; ++m )
	    v[ m ] += G[ i ][ n ][ m ];
	} else if ( p[ n ] > 0 ) {
	  Size i_1 = index( p - Point::base( n ) );
	  for ( unsigned int m = 0; m < M; ++m )
	    v[ m ] -= G[ i_1 ][ n ][ m ];
	}
      return v;
    }

    /// @return the scalar form lambda.u
    ValueForm multiplication( Scalar lambda, const ValueForm& u ) const
    {
      ValueForm S = u;
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] *= lambda;
      return S;
    }
    
    /// @return the scalar form u - v
    ValueForm subtraction( const ValueForm& u, const ValueForm& v ) const
    {
      ValueForm S = u;
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] -= v[ i ];
      return S;
    }

    /// u -= v
    void subtract( ValueForm& u, const ValueForm& v ) const
    {
      for ( Size i = 0; i < u.size(); ++i )
	u[ i ] -= v[ i ];
    }

    /// @return the scalar form a.u + b.v
    ValueForm combination( const Scalar a, const ValueForm& u,
			   const Scalar b, const ValueForm& v ) const
    {
      ValueForm S( u.size() );
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] = a * u[ i ] + b * v[ i ];
      return S;
    }

    /// @return the Total Variation of _u (current approximation of image _I).
    Scalar energyTV() const
    {
      ScalarForm energies = norm( grad( _u ) );
      Scalar            E = 0.0;
      for ( Scalar e_vtx : energies ) E += e_vtx;
      return E;
    }
    
    /// Does one pass of TV regularization (u, p and I must have the
    /// meaning of the previous iteration).
    /// @note Chambolle, Pock primal-dual algorithm 1
    Scalar optimize( Scalar lambda,
		     Scalar dt = 0.248, Scalar tol = 0.01, int max_iter = 15 )
    {
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      //trace.info() << "lambda.f" << std::endl;
      VectorValueForm p( _p.size() ); // this is p^{n+1}
      ValueForm      lf = multiplication( lambda, _I );
      Scalar     diff_p = 0.0;
      int          iter = 0; // iteration number
      do {
	// trace.info() << "div( p ) - lambda.f" << std::endl;
	ValueForm        dp_lf = subtraction( div( _p ), lf );
	// trace.info() << "G := grad( div( p ) - lambda.f)" << std::endl;
	VectorValueForm gdp_lf = grad( dp_lf );
	// trace.info() << "N := | G |" << std::endl;
	ScalarForm     ngdp_lf = norm( gdp_lf );
	// trace.info() << "p^n+1 := ( p + dt * G ) / ( 1 + dt | G | )" << std::endl;
	diff_p = 0.0;
	for ( Size i = 0; i < p.size(); i++ ) {
	  Scalar alpha = 1.0 / ( 1.0 + dt * ngdp_lf[ i ] );
	  if ( alpha <= 0.0 )
	    trace.warning() << "Face " << i << " alpha=" << alpha << std::endl;
	  VectorValue delta;
	  for ( unsigned int n = 0; n < N; ++n ) {
	    p[ i ][ n ] = alpha * ( _p[ i ][ n ] + dt * gdp_lf[ i ][ n ] );
	    delta[ n ]  = p[ i ][ n ] - _p[ i ][ n ];
	  }
	  diff_p = std::max( diff_p, normY( delta ) );
	}
	trace.info() << "Iter n=" << (iter++) << " diff_p=" << diff_p
		     << " tol=" << tol << std::endl;
	std::swap( p, _p );
      } while ( ( diff_p > tol ) && ( iter < max_iter ) );
      _u = combination( 1.0, _I, -1.0/lambda, div( _p ) );
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      return diff_p;
    }

    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ImageTVRegularization



  template <typename TSpace, int M>
  class ImageTVZoom : public ImageTVRegularization< TSpace, M >
  {
  public:
    typedef TSpace                            Space;
    typedef ImageTVZoom< TSpace, M>           Self;
    typedef ImageTVRegularization< TSpace, M> Base;
    using typename Base::Point;
    using typename Base::Vector;
    using typename Base::RealPoint;
    using typename Base::RealVector;
    using typename Base::Size;
    using typename Base::Domain;
    using typename Base::Scalar;
    using typename Base::Value;
    using typename Base::VectorValue;
    using typename Base::ScalarForm;
    using typename Base::ValueForm;
    using typename Base::VectorValueForm;
    using Base::N;
    using Base::_I;
    using Base::_domain;
    using Base::_extent;
    using Base::_u;
    using Base::_p;
    using Base::grad;
    using Base::div;
    using Base::norm;
    using Base::normX;
    using Base::normY;
    using Base::multiplication;
    using Base::subtraction;
    using Base::subtract;
    using Base::combination;
    using Base::energyTV;
    using Base::index;
    using Base::point;
    
    /// The zoom factor (2,3,...)
    int        _zoom;
    /// The unzoomed domain.
    Domain     _uz_domain;
    /// The unzoomed extent.
    Vector     _uz_extent;
    /// The scalar form v
    ValueForm  _v;
  public:
    ImageTVZoom() : Base() {}

    /// @tparam Functor the type of function: Image Value x int --> Scalar,
    /// where int is the dimension in the Value.
    /// @see Color2ValueFunctor
    /// @see GrayLevel2ValueFunctor
    template <typename Image, typename Functor>
    void init( const Image& I, Functor f, int zoom = 2 )
    {
      _zoom = zoom;
      ASSERT( _zoom >= 2 );
      _I.reserve( I.size() );
      _uz_domain = I.domain();
      _uz_extent = I.extent();
      _domain    = Domain( Point::zero, I.domain().upperBound() * zoom );
      _extent    = I.domain().upperBound() * zoom + Point::diagonal( 1 );
      const Size z_size = _domain.size();
      _p.resize( z_size ); // p = vec(0)
      _v.resize( z_size ); // v = 0  at initialization
      _u.resize( z_size ); // u = f at sampled points.
      Value v;
      Domain zd( Point::zero, Point::diagonal( zoom - 1 ) );
      for ( auto p : I.domain() ) {
	auto val_I = I( p );
	for ( unsigned int m = 0; m < M; ++m )
	  v[ m ] = f( val_I, m );
	_I.push_back( v );
	Point b = zoom * p;
	for ( Point q : zd ) {
	  Point r = b + q;
	  if ( r.sup( _domain.upperBound() ) == _domain.upperBound() )
	    _u[ index( b + q ) ] = v;
	}
      }
    }

    /// Does one pass of TV regularization (u, p and I must have the
    /// meaning of the previous iteration).
    ///
    /// @note Li, Bao, Liu, and Zhang algorithm, adaptation of Pock
    /// primal-dual algorithm 1
    Scalar optimize( Scalar lambda, Scalar theta,
		     Scalar dt = 0.248, Scalar tol = 0.01, int max_iter = 15 )
    {
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      //trace.info() << "lambda.f" << std::endl;
      VectorValueForm p( _p.size() ); // this is p^{n+1}
      Scalar          lt = lambda * theta;
      Scalar      diff_p = 0.0;
      int           iter = 0; // iteration number
      Scalar           h = 1.0 / (Scalar) _zoom;
      ValueForm       dp = multiplication( 1.0/h, div( _p ) );
      do {
	// Computing p^{k+1}
	ValueForm        dp_ut = combination( 1.0, dp,
					      1.0/(theta), _u );
	VectorValueForm gdp_ut = grad( dp_ut );
	ScalarForm     ngdp_ut = norm( gdp_ut );
	diff_p = 0.0;
	for ( Size i = 0; i < p.size(); i++ ) {
	  Scalar alpha = 1.0 / ( 1.0 + dt * ngdp_ut[ i ] );
	  if ( alpha <= 0.0 )
	    trace.warning() << "Face " << i << " alpha=" << alpha << std::endl;
	  VectorValue delta;
	  for ( unsigned int n = 0; n < N; ++n ) {
	    p[ i ][ n ] = alpha * ( _p[ i ][ n ] + dt * gdp_ut[ i ][ n ] );
	    delta[ n ]  = p[ i ][ n ] - _p[ i ][ n ];
	  }
	  diff_p = std::max( diff_p, normY( delta ) );
	}
	trace.info() << "Iter n=" << (iter++) << " diff_p=" << diff_p
		     << " tol=" << tol << std::endl;
	dp = div( p );
	// Computing v^{k+1}
	_v = combination( 1.0, _u, -theta, dp );
	// Computing u^{k+1}
	Size i = 0;
	for ( Point p : _domain ) {
	  Point q = p / _zoom;
	  if ( ( q * _zoom ) == p )
	    // Sample point
	    _u[ i ] = ( lt * _I[ uzindex( q ) ] + _v[ i ] ) / ( 1.0 + lt );
	  else
	    _u[ i ] = _v[ i ];
	  i++;
	}
	// Preparing next iteration
	std::swap( p, _p );
      } while ( ( diff_p > tol ) && ( iter < max_iter ) );
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      return diff_p;
    }

    
    /// @see Value2ColorFunctor
    /// @see Value2GrayLevelFunctor
    template <typename Image, typename Functor>
    bool outputU( Image& J, Functor f ) const
    {
      J      = Image( _domain ); 
      Size i = 0;
      for ( auto & val_J : J ) {
	val_J = f( _u[ i ] );
	i++;
      }
      return i == _u.size();
    }

    /// Linearize the point to an index.
    Size uzindex( const Point &aPoint ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getIndex( aPoint, _uz_extent );
    }
    /// Delinearize the index to a point.
    Point uzpoint( const Size index ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getPoint( index, _uz_extent );
    }

    
    
  };
} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ImageTVRegularization_h

#undef ImageTVRegularization_RECURSES
#endif // else defined(ImageTVRegularization_RECURSES)
