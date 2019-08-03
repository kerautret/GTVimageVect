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
 * @file LeastSquares.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module LeastSquares.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(LeastSquares_RECURSES)
#error Recursive header files inclusion detected in LeastSquares.h
#else // defined(LeastSquares_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LeastSquares_RECURSES

#if !defined LeastSquares_h
/** Prevents repeated inclusion of headers. */
#define LeastSquares_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include "BezierTriangle2.h"
#include "cairo.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class LeastSquares
  /**
     Description of class 'LeastSquares' <p> \brief Aim: This class
     gives functions to draw triangles in images.

     In this class, colors are represented as a 3-vector with
     components in 0..255.

     @tparam TSpace the digital space for images (choose Z2i::Space).
  */
  template <typename TPoint>
  class LeastSquares {
  public:
    typedef TPoint                        Point;
    typedef typename Point::Coordinate    Scalar;
    BOOST_STATIC_ASSERT (( Point::dimension == 2 ));
    
    static 
    std::array<Scalar,6>
    linearFit( const std::vector< Point >&  X,
	       const std::vector< Scalar >& V )
    {
      bool ok;
      std::array<Scalar,6> S = { V[0], 0, 0, 0, 0, 0 };
      S[ 0 ] = V[ 0 ];
      auto n = std::min( X.size(), V.size() );
      if ( n == 3 ) {
	typedef SimpleMatrix<Scalar,3,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 4 ) {
	typedef SimpleMatrix<Scalar,4,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 5 ) {
	typedef SimpleMatrix<Scalar,5,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 6 ) {
	typedef SimpleMatrix<Scalar,6,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 7 ) {
	typedef SimpleMatrix<Scalar,7,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 8 ) {
	typedef SimpleMatrix<Scalar,8,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 9 ) {
	typedef SimpleMatrix<Scalar,9,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 10 ) {
	typedef SimpleMatrix<Scalar,10,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 11 ) {
	typedef SimpleMatrix<Scalar,11,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 12 ) {
	typedef SimpleMatrix<Scalar,12,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 13 ) {
	typedef SimpleMatrix<Scalar,13,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n >= 14 ) {
	typedef SimpleMatrix<Scalar,14,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } 
      return S;
    }

    static 
    std::array<Scalar,6>
    quadraticFit( const std::vector< Point >&  X,
		  const std::vector< Scalar >& V )
    {
      bool ok;
      std::array<Scalar,6> S = { V[0], 0, 0, 0, 0, 0 };
      auto n = std::min( X.size(), V.size() );
      if ( n == 3 ) {
	typedef SimpleMatrix<Scalar,3,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 4 ) {
	typedef SimpleMatrix<Scalar,4,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareLinearFit( M, Y, X, V );
	ok = computeFit( S, M, Y );
      } else if ( n == 5 ) {
	typedef SimpleMatrix<Scalar,5,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 6 ) {
	typedef SimpleMatrix<Scalar,6,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 7 ) {
	typedef SimpleMatrix<Scalar,7,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 8 ) {
	typedef SimpleMatrix<Scalar,8,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 9 ) {
	typedef SimpleMatrix<Scalar,9,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 10 ) {
	typedef SimpleMatrix<Scalar,10,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 11 ) {
	typedef SimpleMatrix<Scalar,11,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 12 ) {
	typedef SimpleMatrix<Scalar,12,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n == 13 ) {
	typedef SimpleMatrix<Scalar,13,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      } else if ( n >= 14 ) {
	typedef SimpleMatrix<Scalar,14,2>      Matrix;
	typedef typename Matrix::ColumnVector ColumnVector;
	Matrix M; ColumnVector Y;
	prepareQuadraticFit( M, Y, X, V );
	if ( ( ok = computeFit( S, M, Y ) ) == false )
	  return linearFit( X, V ); 
      }
      return S;
    }
    
    template <typename Matrix, typename ColumnVector>
    static
    bool computeFit( std::array<Scalar,6>& S,
		     const Matrix& M, const ColumnVector& Y )
    {
      const auto  tM = M.transpose();
      const auto tMM = tM * M;
      if ( tMM.determinant() != 0 ) {
	const auto R = tMM.inverse() * ( tM * Y );
	S[ 1 ] = R[ 0 ]; S[ 2 ] = R[ 1 ];
	if ( R.dimension > 2 ) {
	  S[ 3 ] = R[ 2 ];
	  S[ 4 ] = R[ 3 ];
	  S[ 5 ] = R[ 4 ];
	}
	return true;
      }
      return false;
    }
    
    template <typename Matrix, typename ColumnVector>
    static
    void prepareLinearFit( Matrix& M, ColumnVector& Y,
			   const std::vector< Point >&  X,
			   const std::vector< Scalar >& V )
    {
      for ( Dimension i = 0; i < M.rows(); ++i ) {
	M.setComponent( i, 0, ( X[i][0] - X[0][0] ) );
	M.setComponent( i, 1, ( X[i][1] - X[0][1] ) );
	Y[ i ] = V[ i ] - V[ 0 ];
      }	
    }

    template <typename Matrix, typename ColumnVector>
    static
    void prepareQuadraticFit( Matrix& M, ColumnVector& Y,
			   const std::vector< Point >&  X,
			   const std::vector< Scalar >& V )
    {
      for ( Dimension i = 0; i < M.rows(); ++i ) {
	M.setComponent( i, 0, ( X[i][0] - X[0][0] ) );
	M.setComponent( i, 1, ( X[i][1] - X[0][1] ) );
	M.setComponent( i, 2, ( X[i][0] - X[0][0] ) * ( X[i][0] - X[0][0] ) );
	M.setComponent( i, 3, ( X[i][0] - X[0][0] ) * ( X[i][1] - X[0][1] ) );
	M.setComponent( i, 4, ( X[i][1] - X[0][1] ) * ( X[i][1] - X[0][1] ) );
	Y[ i ] = V[ i ] - V[ 0 ];
      }	
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class LeastSquares

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LeastSquares_h

#undef LeastSquares_RECURSES
#endif // else defined(LeastSquares_RECURSES)
