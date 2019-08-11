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
 * @file BezierCurve.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module BezierCurve.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(BezierCurve_RECURSES)
#error Recursive header files inclusion detected in BezierCurve.h
#else // defined(BezierCurve_RECURSES)
/** Prevents recursive inclusion of headers. */
#define BezierCurve_RECURSES

#if !defined BezierCurve_h
/** Prevents repeated inclusion of headers. */
#define BezierCurve_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <functional>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/linalg/SimpleMatrix.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class BezierCurve
  /**
     Description of class 'BezierCurve' <p> \brief Aim: This class
     represents a Bezier triangle of degree 2, as a
     polynomial R^2 -> R. 

     @tparam TSpace the digital space for images (choose Z2i::Space).
  */
  template <typename TSpace>
  class BezierCurve
  {
  public:
    typedef TSpace                            Space;
    typedef BezierCurve< TSpace> Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef typename Vector::Component            Integer;
    typedef typename RealVector::Component        Scalar;
    typedef std::function< Point( RealPoint p ) > Digitizer;

    
    /// The Bezier points of the curve
    std::vector<RealPoint>   _bpoints;

    
    // ----------------------- Standard services ------------------------------
  public:
  
    /// Destructor.
    ~BezierCurve() = default;
    /// Default constructor. The object is invalid.
    BezierCurve() = default;
    /// Default copy constructor. 
    BezierCurve( const BezierCurve& ) = default;
    /// Default move constructor. 
    BezierCurve( BezierCurve&& ) = default;
    /// Default assignment.
    BezierCurve& operator=( BezierCurve&& ) = default;

    /// Constructeur from Bezier points and values in order b200,
    /// b020, b002, b110, b011, b101
    BezierCurve( const std::vector<RealPoint>& bpoints )
      : _bpoints( bpoints ) {}

    /// @return the i-th Bezier control point (in order b200, b020,
    /// b002, b110, b011, b101).
    const RealPoint& b( int i ) const { return _bpoints[ i ]; }

    /// @param[in] t the parameter (between 0 and 1) where you want to
    /// evaluate Bezier curve.
    ///
    /// @return the corresponding point B(t)
    RealPoint operator()( Scalar t ) const
    {
      if ( _bpoints.size() == 1 ) return _bpoints[ 0 ];
      const Size n = _bpoints.size() - 1;
      std::vector<RealPoint> tP( n );
      for ( Size i = 0; i < n; i++ )
	tP[ i ] = (1.0-t) * b( i ) + t*b( i+1 );
      return Self( tP ) ( t );
    }

      
    /// @param[out] lo the lowest digital point of the digital bounding box
    /// @param[out] hi the highet digital point of the digital bounding box
    void digitalBoundingBox( const Digitizer& dig,
			     Point& lo, Point& hi ) const
    {
      hi = lo = dig( b( 0 ) );
      for ( Size i = 1; i < _bpoints.size(); ++i ) {
	lo = lo.inf( dig( b( i ) ) );
	hi = hi.sup( dig( b( i ) ) );
      }
    }
    
    /// @param[out] lo the lowest point of the bounding box
    /// @param[out] hi the highet point of the bounding box
    void boundingBox( RealPoint& rlo, RealPoint& rhi ) const
    {
      rlo = b( 0 );
      rhi = b( 0 );
      for ( Size i = 1; i < _bpoints.size(); ++i ) {
	rlo = rlo.inf( b( i ) );
	rhi = rhi.sup( b( i ) );
      }
    }

    /// @param[out] p a digital point traced by the Bezier curve.
    /// @return 'true' iff the Bezier curve is reduced to one digital point.
    bool isOnePoint( const Digitizer& dig, Point& p ) const
    {
      Point hi;
      digitalBoundingBox( dig, p, hi );
      return p == hi;
    }

    /// @param[out] p a digital point traced by the Bezier curve.
    /// @param[out] q a digital point traced by the Bezier curve.
    ///
    /// @return 'true' iff the Bezier curve is reduced to a group of
    /// at most four digital points.
    bool isTwoPoints( const Digitizer& dig, Point& p, Point& q ) const
    {
      digitalBoundingBox( dig, p, q );
      return ( ( p[ 0 ] + 1 ) >= q[ 0 ] )
	&&   ( ( p[ 1 ] + 1 ) >= q[ 1 ] );
    }

    /// Traces the digital points covered by the Bezier curve.
    ///
    /// @param[out] dp the vector of digital points traced by the Bezier curve.
    /// @param[out] rp the vector of associated points traced by the Bezier curve.
    /// @param[out] dt the associated vector of parameters t, ie `dp[i]=B(t[i])`.
    void trace( const Digitizer& dig,
		std::vector<Point>& dp,
		std::vector<RealPoint>& rp,
		std::vector<Scalar>& dt ) const
    {
      std::vector<Scalar> dd;
      trace( dig, dp, rp, dt, dd, 0.0, 1.0 );
    }

    /// Traces the digital points covered by the Bezier curve.
    ///
    /// @param[out] dp the vector of digital points traced by the Bezier curve.
    /// @param[out] rp the vector of associated points traced by the Bezier curve.
    /// @param[out] dt the associated vector of parameters t, ie `dp[i]=B(t[i])`.
    void traceDirect( const Digitizer& dig,
    		      std::vector<Point>& dp,
    		      std::vector<RealPoint>& rp,
    		      std::vector<Scalar>& dt ) const
    {
      std::vector<Scalar> dd;
      traceDirect( dig, dp, rp, dt, dd );
    }
    
    /// Traces the digital points covered by the Bezier curve.
    ///
    /// @param[out] dp the vector of digital points traced by the Bezier curve.
    /// @param[out] rp the vector of associated points traced by the Bezier curve.
    /// @param[out] dt the associated vector of parameters t, ie `dp[i]=B(t[i])`.
    /// @param[out] dd the associated error in the digitization.
    ///
    /// @param[in] t0, t1 the current interval [t0,t1] in the
    /// recursion (should be called with [0,1]).
    void traceDirect( const Digitizer& dig,
		      std::vector<Point>& dp,
		      std::vector<RealPoint>& rp,
		      std::vector<Scalar>& dt,
		      std::vector<Scalar>& dd ) const
    {
      // Add first point
      const Scalar t0 = 0.0;
      const Scalar t1 = 1.0;
      const RealPoint rp0 = b( 0 ); 
      const RealPoint rp1 = b( _bpoints.size() - 1 ); 
      addPoint( dig, rp0, t0, dp, rp, dt, dd );
      Scalar delta_t = 1.0 / ( (dig( rp1 ) - dig( rp0 )).norm() + 1.0 );
      Scalar t = t0 - 1.e-4;
      while ( t + delta_t < t1 ) {
	RealPoint r = (*this)( t + delta_t );
	Point    dr = dig( r );
	auto      d = ( dr - dp.back() ).normInfinity();
	if ( d <= 1 ) {
	  addPoint( dig, r, t + delta_t, dp, rp, dt, dd );
	  t += delta_t;
	  // if ( d == 0 ) delta_t *= 1.5;
	} else {
	  delta_t *= 0.5;
	}
      }
      addPoint( dig, rp1, t1, dp, rp, dt, dd );
    }
    
    bool addPoint( const Digitizer& dig,
		   const RealPoint& r,
		   Scalar t,
		   std::vector<Point>& dp,
		   std::vector<RealPoint>& rp,
		   std::vector<Scalar>& dt,
		   std::vector<Scalar>& dd ) const
    {
      Point      dr = dig( r );
      RealPoint rdr = RealPoint( dr[ 0 ], dr[ 1 ] );
      if ( dp.empty() || ( dr != dp.back() ) ) {
	dp.push_back( dr );
	rp.push_back( r );
	dd.push_back( (rdr - r).norm() );
	dt.push_back( t );
	return true;
      } else if ( (rdr - r).norm() < dd.back() ) {
	rp.back() = r;
	dd.back() = (rdr - r).norm();
	dt.back() = t;
      }
      return false;
    }
    
    /// Traces the digital points covered by the Bezier curve.
    ///
    /// @param[out] dp the vector of digital points traced by the Bezier curve.
    /// @param[out] rp the vector of associated points traced by the Bezier curve.
    /// @param[out] dt the associated vector of parameters t, ie `dp[i]=B(t[i])`.
    /// @param[out] dd the associated error in the digitization.
    ///
    /// @param[in] t0, t1 the current interval [t0,t1] in the
    /// recursion (should be called with [0,1]).
    void trace( const Digitizer& dig,
		std::vector<Point>& dp,
		std::vector<RealPoint>& rp,
		std::vector<Scalar>& dt,
		std::vector<Scalar>& dd,
		Scalar t0 = 0.0, Scalar t1 = 1.0 ) const
    {
      const Size    n = _bpoints.size();
      Point p, q;
      if ( isTwoPoints( dig, p, q ) ) {
	Scalar        t = t0;
	const Scalar ti = (t1-t0) / ( _bpoints.size() - 1 );
	for ( Size i = 0; i < n; i += n-1 ) {
	  RealPoint   r = _bpoints[ i ];
	  Point      dr = dig( r );
	  RealPoint rdr = RealPoint( dr[ 0 ], dr[ 1 ] );
	  if ( dp.empty() || ( dr != dp.back() ) ) {
	    dp.push_back( dr );
	    rp.push_back( r );
	    dd.push_back( (rdr - r).norm() );
	    dt.push_back( t );
	  } else if ( (rdr - r).norm() < dd.back() ) {
	    rp.back() = r;
	    dd.back() = (rdr - r).norm();
	    dt.back() = t;
	  }
	  t += ti * (n-1);
	}
      } else {
	const Scalar tm = (t0+t1)/2.0;
	std::vector<RealPoint> left ( n );
	std::vector<RealPoint> right( n );
	std::vector<RealPoint> mid = _bpoints;
	for ( Size i = 0; i < n; i++ ) {
	  left [ i ]     = mid[ 0 ];
	  right[ n-1-i ] = mid[ n-1-i ];
	  for ( Size j = 0; j < (n-i-1); j++ )
	    mid[ j ]     = 0.5 * ( mid[ j ] + mid[ j+1 ] );
	}
	Self bleft ( left  );
	Self bright( right );
	bleft .trace( dig, dp, rp, dt, dd, t0, tm );
	bright.trace( dig, dp, rp, dt, dd, tm, t1 );
      }
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class BezierCurve

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined BezierCurve_h

#undef BezierCurve_RECURSES
#endif // else defined(BezierCurve_RECURSES)
