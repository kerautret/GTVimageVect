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
 * @file BezierTriangle2.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module BezierTriangle2.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(BezierTriangle2_RECURSES)
#error Recursive header files inclusion detected in BezierTriangle2.h
#else // defined(BezierTriangle2_RECURSES)
/** Prevents recursive inclusion of headers. */
#define BezierTriangle2_RECURSES

#if !defined BezierTriangle2_h
/** Prevents repeated inclusion of headers. */
#define BezierTriangle2_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/linalg/SimpleMatrix.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class BezierTriangle2
  /**
     Description of class 'BezierTriangle2' <p> \brief Aim: This class
     represents a Bezier triangle of degree 2, as a
     polynomial R^2 -> R^M. 

     @tparam TSpace the digital space for images (choose Z2i::Space).

     @tparam M the number of scalar per value, generally 1 for
     gray-level images and 3 for color images.
  */
  template <typename TSpace, int M>
  class BezierTriangle2
  {
  public:
    typedef TSpace                            Space;
    typedef BezierTriangle2< TSpace, M> Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));
    BOOST_STATIC_ASSERT (( M >= 1 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< M, Scalar >              Value;
    typedef PointVector< 3, Scalar >              RealPoint3;

    typedef std::function< Scalar( Scalar, Scalar, Scalar ) > BPolynomial;
    
    /// The six bezier points, in order b200, b020, b002, b110, b011, b101
    std::array<RealPoint,6>   _bpoints;
    /// The values at the six bezier points, in order b200, b020,
    /// b002, b110, b011, b101
    std::array<Value,6>       _bvalues;

    static inline
    Scalar B200( Scalar r, Scalar  , Scalar   ) { return r*r; }
    static inline
    Scalar B020( Scalar  , Scalar s, Scalar   ) { return s*s; }
    static inline
    Scalar B002( Scalar  , Scalar  , Scalar t ) { return t*t; }
    static inline
    Scalar B110( Scalar r, Scalar s, Scalar   ) { return r*s; }
    static inline
    Scalar B011( Scalar  , Scalar s, Scalar t ) { return s*t; }
    static inline
    Scalar B101( Scalar r, Scalar  , Scalar t ) { return r*t; }

    static inline
    Scalar det( const RealVector& v, const RealVector& w )
    {
      return v[ 0 ] * w[ 1 ] - v[ 1 ] * w[ 0 ];
    }
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /// Destructor.
    ~BezierTriangle2() = default;
    /// Default constructor. The object is invalid.
    BezierTriangle2() = default;
    /// Default copy constructor. 
    BezierTriangle2( const BezierTriangle2& ) = default;
    /// Default move constructor. 
    BezierTriangle2( BezierTriangle2&& ) = default;
    /// Default assignment.
    BezierTriangle2& operator=( BezierTriangle2&& ) = default;

    /// Constructeur from Bezier points and values in order b200,
    /// b020, b002, b110, b011, b101
    BezierTriangle2( const std::array<RealPoint,6>& bpoints,
		     const std::array<Value,6>&     bvalues )
      : _bpoints( bpoints ), _bvalues( bvalues ) {}
    
    /// Constructeur from Bezier points b200, b020, b002, and
    /// orthogonality conditions at b110, b011, b101 and values at
    /// b200, b020, b002
    ///
    /// @note The orthogonality conditions is that at b110 (resp. b011
    /// and b101), the gradient of the Bezier polynomial should be
    /// orthogonal to orthvecs[ 0 ] (resp. orthvecs[ 1 ] and orthvecs[
    /// 2 ]).
    BezierTriangle2( const std::array<RealPoint,3> & bpoints,
		     const std::array<RealVector,3>& orthvecs,
		     const std::array<Value,3>     & bvalues )
    {
      _bpoints[ 0 ] = bpoints[ 0 ];
      _bpoints[ 1 ] = bpoints[ 1 ];
      _bpoints[ 2 ] = bpoints[ 2 ];
      _bvalues[ 0 ] = bvalues[ 0 ];
      _bvalues[ 1 ] = bvalues[ 1 ];
      _bvalues[ 2 ] = bvalues[ 2 ];
      RealVector de[ 3 ] = { decompose( orthvecs[ 0 ], 0 ),
			     decompose( orthvecs[ 1 ], 1 ),
			     decompose( orthvecs[ 2 ], 2 ) };
      _bpoints[ 3 ] = 0.5 * ( _bpoints[ 0 ] + _bpoints[ 1 ] );
      _bpoints[ 4 ] = 0.5 * ( _bpoints[ 1 ] + _bpoints[ 2 ] );
      _bpoints[ 5 ] = 0.5 * ( _bpoints[ 2 ] + _bpoints[ 0 ] );
      for ( int i = 0; i < 3; ++i ) {
	std::cout << "de[" << i << "] = " << de[ i ] << std::endl;
	if ( de[ i ][ 1 ] == 0.0 )
	  trace.error() << "Invalid base for BezoutTriangle2." << std::endl;
      }
      Scalar d_over_e[ 3 ] = { de[ 0 ][ 0 ] / de[ 0 ][ 1 ],
				de[ 1 ][ 0 ] / de[ 1 ][ 1 ],
				de[ 2 ][ 0 ] / de[ 2 ][ 1 ] };
      // 0.5*( -a*(d2/e2) + b*(d1/e1 + 1) + c*(d2/e2 - d1/e1 + 1) )
      _bvalues[ 3 ] =
	0.5 * ( -d_over_e[ 2 ] * bvalues[ 0 ]
		+ ( d_over_e[ 1 ] + 1 ) * bvalues[ 1 ]
		+ ( d_over_e[ 2 ] - d_over_e[ 1 ] + 1 ) * bvalues[ 2 ] );
	
      // 0.5*( a*(d0/e0 - d2/e2 + 1) - b*d0/e0 + c*(d2/e2 + 1) )
      _bvalues[ 4 ] =
	0.5 * ( -d_over_e[ 0 ] * bvalues[ 1 ]
		+ ( d_over_e[ 2 ] + 1 ) * bvalues[ 2 ]
		+ ( d_over_e[ 0 ] - d_over_e[ 2 ] + 1 ) * bvalues[ 0 ] );

      // 0.5*( a*(d0/e0 + 1) + b*(d1/e1 - d0/e0 + 1) - c*d1/e1 )
      _bvalues[ 5 ] =
	0.5 * ( -d_over_e[ 1 ] * bvalues[ 2 ]
		+ ( d_over_e[ 0 ] + 1 ) * bvalues[ 0 ]
		+ ( d_over_e[ 1 ] - d_over_e[ 0 ] + 1 ) * bvalues[ 1 ] );
      for ( int i = 0; i < 6; ++i )
	std::cout << "b(" << i << ")=" << b(i) << " v(i)=" << v(i) << std::endl;
    }

    /// @return the i-th Bezier control point (in order b200, b020,
    /// b002, b110, b011, b101).
    const RealPoint& b( int i ) const { return _bpoints[ i ]; }

    /// @return the i-th Bezier value point (in order b200, b020,
    /// b002, b110, b011, b101).
    const Value& v( int i ) const { return _bvalues[ i ]; }

    /// Decomposes a vector p onto the base b[i]b[i+1], b[i]b[i+1],
    /// indices taken modulo 3.
    RealVector decompose( const RealVector& p, int i ) const
    {
      typedef SimpleMatrix<Scalar,2,2>      Matrix;
      typedef typename Matrix::ColumnVector ColumnVector;
      int       j = (i+1)%3;
      int       k = (i+2)%3;
      Matrix    P = { b(j)[ 0 ] - b(i)[ 0 ], b(j)[ 1 ] - b(i)[ 1 ],
		      b(k)[ 0 ] - b(i)[ 0 ], b(k)[ 1 ] - b(i)[ 1 ] };
      Matrix    Q = P.inverse().transpose();
      ColumnVector U = { p[ 0 ], p[ 1 ] };
      ColumnVector V = Q * U;
      return RealVector( V[ 0 ], V[ 1 ] );
    }

    /// @param[in] p the point where you want to evaluate the
    /// BezierTriangle.
    ///
    /// @return the value at this point
    Value operator()( const RealPoint& p ) const
    {
      return operator()( barycentric( p ) );
    }
    
    /// @param[in] b the barycentric coordinates where you want to
    /// evaluate the BezierTriangle: (1,0,0) is b200, (0,1,0) is b020,
    /// (0,0,1) is b002, (a,b,c) with a+b+c=1, 0<=a<=1, 0<=b<=1,
    /// 0<=c<=1 is anywhere between these points.
    ///
    /// @return the value at this point
    Value operator()( const RealPoint3& bc ) const
    {
      return B200( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 0 )
	+    B020( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 1 )
	+    B002( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 2 )
	+    B110( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 3 )
	+    B011( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 4 )
	+    B101( bc[ 0 ], bc[ 1 ], bc[ 2 ] ) * v( 5 );
    }

    /// @param[in] r,s,t the barycentric coordinates where you want to
    /// evaluate the BezierTriangle: (1,0,0) is b200, (0,1,0) is b020,
    /// (0,0,1) is b002, (a,b,c) with a+b+c=1, 0<=a<=1, 0<=b<=1,
    /// 0<=c<=1 is anywhere between these points.
    ///
    /// @return the value at this point
    Value operator()( Scalar r, Scalar s, Scalar t ) const
    {
      return B200( r, s, t ) * v( 0 )
	+    B020( r, s, t ) * v( 1 )
	+    B002( r, s, t ) * v( 2 )
	+    B110( r, s, t ) * v( 3 )
	+    B011( r, s, t ) * v( 4 )
	+    B101( r, s, t ) * v( 5 );
    }

    /// @return the barycentric coordinates of point \a p in this triangle.
    RealPoint3 barycentric( const RealPoint& p ) const
    {
      RealPoint3 bc = { det( b( 1 ) - p, b( 2 ) - p ),
			det( b( 2 ) - p, b( 0 ) - p ),
			det( b( 0 ) - p, b( 1 ) - p ) };
      return bc / ( bc[ 0 ] + bc[ 1 ] + bc[ 2 ] );
    }

    /// @param bc any barycentric coordinates (sum is 1).
    ///
    /// @return 'true' iff the given barycentric coordinates \a bc
    /// indicates a point inside or on the boundary of the triangle.
    bool isInTriangle( const RealPoint3& bc ) const
    {
      return ( 0 <= bc[ 0 ] ) && ( bc[ 0 ] <= 1 )
	&&   ( 0 <= bc[ 1 ] ) && ( bc[ 1 ] <= 1 )
	&&   ( 0 <= bc[ 2 ] ) && ( bc[ 2 ] <= 1 );
    }

      
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class BezierTriangle2

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined BezierTriangle2_h

#undef BezierTriangle2_RECURSES
#endif // else defined(BezierTriangle2_RECURSES)
