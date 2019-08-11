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
 * @file CairoViewer.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module CairoViewer.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CairoViewer_RECURSES)
#error Recursive header files inclusion detected in CairoViewer.h
#else // defined(CairoViewer_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CairoViewer_RECURSES

#if !defined CairoViewer_h
/** Prevents repeated inclusion of headers. */
#define CairoViewer_h

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
  // class CairoViewer
  /**
     Description of class 'CairoViewer' <p> \brief Aim: This class
     gives functions to draw triangles in images.

     In this class, colors are represented as a 3-vector with
     components in 0..255.

     @tparam TSpace the digital space for images (choose Z2i::Space).
  */
  template <typename TSpace>
  class CairoViewer {
  public:
    typedef TSpace                            Space;
    typedef CairoViewer< TSpace > Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef HyperRectDomain<Space>                Domain;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< 3, Scalar>               RealPoint3;
    typedef PointVector< 3, Scalar>               RealVector3;
    typedef PointVector< 3, Scalar >              Value; // Color
    typedef PointVector< 3, Scalar >              Color;
    
  protected:
    const double _redf, _greenf, _bluef;
    int _x0, _y0;
    int _width, _height;
    Domain _draw_domain;
    double _xf, _yf;
    int _shading; // 0: flat, 1: Gouraud, 2: LinearGradient, 3: Discontinuities
    bool _color;
    cairo_surface_t* _surface;
    cairo_t* _cr;
    double _st;        ///< discontinuity stiffness.
    double _am;        ///< discontinuity amplitude.
    double _s0, _sm, _s1; ///< precomputed abscissae from stiffness.

  public:
    /**
       Constructor. 
    */
    CairoViewer( double x0, double y0, double x1, double y1,
		 double xfactor = 1.0, double yfactor = 1.0,
		 int shading = 0,
		 bool color = true,
		 double disc_stiffness = 0.5,
		 double disc_amplitude = 0.75 )
      : _redf( 1.0/255.0f ), _greenf( 1.0/255.0f ), _bluef( 1.0/255.0f ),
      _x0( round( x0 ) ), _y0( round( y0 ) ),
      _width( round( (x1-x0) * xfactor + 1 ) ),
      _height(round( (y1-y0) * xfactor + 1 ) ),
      _draw_domain( Point( _x0, _y0 ),
		    Point( _x0 + _width - 1, _y0 + _height - 1 ) ),
      /* _draw_domain( Point( _x0, _y0 ), */
      /* 		    Point( _x0 + _width - 2, _y0 + _height - 2 ) ), */
      _xf( xfactor ), _yf( yfactor ), _shading( shading ),
      _color( color ), _st( disc_stiffness ), _am( disc_amplitude )
    {
      _surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
					     _width, _height );
      _cr = cairo_create ( _surface );
      // Fill the background with black
      cairo_set_source_rgba( _cr, 0.0, 0.0, 0.0, 1.0 );
      cairo_rectangle ( _cr, 0, 0, _width, _height );
      cairo_fill( _cr );
      _s0 = _st * 0.5;
      _sm = 0.5;
      _s1 = 1.0 - _st * 0.5;
      startAdd();
    }
    
    /// Destructor.
    ~CairoViewer()
    {
      cairo_destroy( _cr );
      cairo_surface_destroy( _surface );
    }

    void startAdd()
    {
      cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
      cairo_set_line_width( _cr, 0.0 ); 
      cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
      cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
    }

    void save( const char* file_name ) const
    {
      cairo_surface_write_to_png( _surface, file_name );
    }

    inline double i( double x ) const
    {
      //return ( (x+0.5) * _xf ) - _x0;
      // Avoids bad approximations around 1/(_xf*_yf) pixels
      return x * _xf - _x0 + 0.5;
    }
    inline double x( double i ) const
    {
      return ( i - 0.5 + _x0 ) / _xf;
    }
    
    inline double j( double y ) const
    {
      //return _height - (( (y+0.5) * _yf ) - _y0) - 1;
      // Avoids bad approximations around 1/(_xf*_yf) pixels
      return _height - ( y * _yf - _y0) - 0.5;
    }
    inline double y( double j ) const
    {
      return ( _height + _y0 - j - 0.5 ) / _yf;
    }

    inline RealPoint ij( RealPoint a_xy ) const
    {
      return RealPoint( i( a_xy[ 0 ] ), j( a_xy[ 1 ] ) );
    }

    inline RealPoint xy( RealPoint a_ij ) const
    {
      return RealPoint( x( a_ij[ 0 ] ), y( a_ij[ 1 ] ) );
    }
    
    bool computeLinearGradient( RealPoint a, RealPoint b, RealPoint c,
				RealVector3 V,
				RealPoint& s, RealPoint& mid, RealPoint& e,
				Scalar& gs, Scalar& gmid, Scalar& ge )
    {
      const RealVector3 One = RealVector3::diagonal( 1 );
      const RealVector3   X = RealVector3( a[ 0 ], b[ 0 ], c[ 0 ] );
      const RealVector3   Y = RealVector3( a[ 1 ], b[ 1 ], c[ 1 ] );
      const RealPoint    Gr = RealPoint( V.crossProduct( Y ).dot( One ),
					 X.crossProduct( V ).dot( One ) );
      if ( Gr == RealPoint::zero ) {
	s   = ij( a );	mid = ij( b );	e   = ij( c );
	gs  = V[ 0 ];	gmid= V[ 1 ];	ge  = V[ 2 ];
	return false;
      }
      const RealPoint Ur = Gr.getNormalized();
      const Scalar    da = Ur.dot( a );
      const Scalar    db = Ur.dot( b );
      const Scalar    dc = Ur.dot( c );
      Scalar td[ 3 ] = { da, db, dc };
      const int middle[3][3] = { { -1, 2, 1 }, { 2, -1, 0 }, { 1, 0, -1 } }; 
      int m = ( da < db ) ? ( ( da < dc ) ? 0 : 2 ) : ( db < dc ? 1 : 2 );
      int M = ( da >= db ) ? ( ( da >= dc ) ? 0 : 2 ) : ( db >= dc ? 1 : 2 );
      int k = middle[ m ][ M ];
      if ( ( k == -1 ) || ( td[ m ] == td[ M ] ) )
	trace.error() << "Invalid mid m=" << m << " k=" << k << " M=" << M
		      << " da=" << da << " db=" << db << " dc=" << dc
		      << " X=" << X << "Y=" << Y << " V=" << V
		      << " Gr=" << Gr << " Ur=" << Ur
		      << std::endl;
      const RealPoint pts[ 3 ] = { a, b, c };
      const RealPoint d1( pts[ m ] );
      s    = RealPoint( i( d1[ 0 ] ), j( d1[ 1 ] ) );
      gs   = V[ m ];
      const RealPoint d2( pts[ m ] + ( pts[ k ] - pts[ m ] ).dot( Ur ) * Ur );
      mid  = RealPoint( i( d2[ 0 ] ), j( d2[ 1 ] ) );
      gmid = V[ k ];
      const RealPoint d3( pts[ m ] + ( pts[ M ] - pts[ m ] ).dot( Ur ) * Ur );
      e    = RealPoint( i( d3[ 0 ] ), j( d3[ 1 ] ) );
      ge   = V[ M ];
      return true;
    }

    void drawLinearGradientTriangle( RealPoint a, RealPoint b, RealPoint c, 
				     Value val_a, Value val_b, Value val_c ) 
    {
      RealPoint s, m, e;
      Scalar  gs, gm, ge;
      Scalar  t;
      cairo_pattern_t *pat;
      if ( _color ) {
	const RealVector3 Vr = Value( val_a[ 0 ], val_b[ 0 ], val_c[ 0 ] );
	const RealVector3 Vg = Value( val_a[ 1 ], val_b[ 1 ], val_c[ 1 ] );
	const RealVector3 Vb = Value( val_a[ 2 ], val_b[ 2 ], val_c[ 2 ] );
	// Draw path
	cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
	cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
	cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
	cairo_close_path( _cr );
	// Draw red
	if ( computeLinearGradient( a, b, c, Vr, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  t = (m-s).norm() / (e-s).norm();
	  cairo_pattern_add_color_stop_rgb (pat, 0, gs * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, t, gm * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, 1, ge * _redf, 0, 0);
	  cairo_set_source( _cr, pat );
	  cairo_fill_preserve( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, gs * _redf, 0, 0 );
	  cairo_fill_preserve( _cr );
	}
	// Draw green
	if ( computeLinearGradient( a, b, c, Vg, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  t = (m-s).norm() / (e-s).norm();
	  cairo_pattern_add_color_stop_rgb (pat, 0, 0, gs * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, t, 0, gm * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, 1, 0, ge * _greenf, 0);
	  cairo_set_source( _cr, pat );
	  cairo_fill_preserve( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, 0, gs * _greenf, 0 );
	  cairo_fill_preserve( _cr );
	}
	// Draw blue
	if ( computeLinearGradient( a, b, c, Vb, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  t = (m-s).norm() / (e-s).norm();
	  cairo_pattern_add_color_stop_rgb (pat, 0, 0, 0, gs * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, t, 0, 0, gm * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, 1, 0, 0, ge * _bluef );
	  cairo_set_source( _cr, pat );
	  cairo_fill( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, 0, 0, gs * _bluef );
	  cairo_fill( _cr );
	}
      } else { // monochrome
	const Value  Vm = Value( val_a[ 0 ], val_b[ 0 ], val_c[ 0 ] );
	// Draw path
	cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
	cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
	cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
	cairo_close_path( _cr );
	// Draw gray-level
	if ( computeLinearGradient( a, b, c, Vm, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  t = (m-s).norm() / (e-s).norm();
	  cairo_pattern_add_color_stop_rgb (pat, 0, gs * _redf, gs * _greenf, gs * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, t, gm * _redf, gm * _greenf, gm * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, 1, ge * _redf, ge * _greenf, ge * _bluef );
	  cairo_set_source( _cr, pat );
	  cairo_fill( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, gs * _redf, gs * _greenf, gs * _bluef );
	  cairo_fill( _cr );
	}
      }
    }

    double disY0( double gs, double ge ) const {
      return std::max( 0.0, std::min( 255.0, _am * gs + (1.0 - _am ) * ge ) );
    }
    double disYm( double gs, double ge ) const {
      return 0.5 * ( gs + ge );
    }
    double disY1( double gs, double ge ) const {
      return std::max( 0.0, std::min( 255.0, _am * ge + (1.0 - _am ) * gs ) );
    }
      
    void drawNonLinearGradientTriangle( RealPoint a, RealPoint b, RealPoint c, 
					Value val_a, Value val_b, Value val_c )
    {
      RealPoint s, m, e;
      Scalar  gs, gm, ge;
      cairo_pattern_t *pat;
      if ( _color ) {
	const Value  Vr = Value( val_a[ 0 ], val_b[ 0 ], val_c[ 0 ] );
	const Value  Vg = Value( val_a[ 1 ], val_b[ 1 ], val_c[ 1 ] );
	const Value  Vb = Value( val_a[ 2 ], val_b[ 2 ], val_c[ 2 ] );
	// Draw path
	cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
	cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
	cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
	cairo_close_path( _cr );
	// Draw red
	if ( computeLinearGradient( a, b, c, Vr, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  cairo_pattern_add_color_stop_rgb (pat, 0.0, gs           * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _s0, disY0(gs,ge) * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _sm, disYm(gs,ge) * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _s1, disY1(gs,ge) * _redf, 0, 0);
	  cairo_pattern_add_color_stop_rgb (pat, 1.0, ge           * _redf, 0, 0);
	  cairo_set_source( _cr, pat );
	  cairo_fill_preserve( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, gs * _redf, 0, 0 );
	  cairo_fill_preserve( _cr );
	}
	// Draw green
	if ( computeLinearGradient( a, b, c, Vg, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  cairo_pattern_add_color_stop_rgb (pat, 0.0, 0, gs           * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _s0, 0, disY0(gs,ge) * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _sm, 0, disYm(gs,ge) * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, _s1, 0, disY1(gs,ge) * _greenf, 0);
	  cairo_pattern_add_color_stop_rgb (pat, 1.0, 0, ge           * _greenf, 0);
	  cairo_set_source( _cr, pat );
	  cairo_fill_preserve( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, 0, gs * _greenf, 0 );
	  cairo_fill_preserve( _cr );
	}
	// Draw blue
	if ( computeLinearGradient( a, b, c, Vb, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  cairo_pattern_add_color_stop_rgb (pat, 0.0, 0, 0, gs           * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _s0, 0, 0, disY0(gs,ge) * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _sm, 0, 0, disYm(gs,ge) * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _s1, 0, 0, disY1(gs,ge) * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, 1.0, 0, 0, ge           * _bluef );
	  cairo_set_source( _cr, pat );
	  cairo_fill( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, 0, 0, gs * _bluef );
	  cairo_fill( _cr );
	}
      } else { // monochrome
	const Value  Vm = Value( val_a[ 0 ], val_b[ 0 ], val_c[ 0 ] );
	// Draw path
	cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
	cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
	cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
	cairo_close_path( _cr );
	// Draw gray-level
	if ( computeLinearGradient( a, b, c, Vm, s, m, e, gs, gm, ge ) ) {
	  pat = cairo_pattern_create_linear(s[0],s[1],e[0],e[1]);
	  cairo_pattern_add_color_stop_rgb(pat, 0.0, gs * _redf, gs * _greenf, gs * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _s0, disY0(gs,ge) * _redf, disY0(gs,ge) * _greenf, disY0(gs,ge) * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _sm, disYm(gs,ge) * _redf, disYm(gs,ge) * _greenf, disYm(gs,ge) * _bluef );
	  cairo_pattern_add_color_stop_rgb (pat, _s1, disY1(gs,ge) * _redf, disY1(gs,ge) * _greenf, disY1(gs,ge) * _bluef);
	  cairo_pattern_add_color_stop_rgb (pat, 1.0, ge * _redf, ge * _greenf, ge * _bluef );
	  cairo_set_source( _cr, pat );
	  cairo_fill( _cr );
	  cairo_pattern_destroy( pat );
	} else {
	  cairo_set_source_rgb( _cr, gs * _redf, gs * _greenf, gs * _bluef );
	  cairo_fill( _cr );
	}
      }
    }
    void drawGouraudTriangle( RealPoint a, RealPoint b, RealPoint c, 
			      Value val_a, Value val_b, Value val_c ) 
    {
      cairo_pattern_t * pattern = cairo_pattern_create_mesh();
      /* Add a Gouraud-shaded triangle */
      cairo_mesh_pattern_begin_patch (pattern);
      cairo_mesh_pattern_move_to (pattern, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_mesh_pattern_line_to (pattern, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_mesh_pattern_line_to (pattern, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 0,
					       val_a[ 0 ] * _redf,
					       val_a[ 1 ] * _greenf,
					       val_a[ 2 ] * _bluef );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 1,
					       val_b[ 0 ] * _redf,
					       val_b[ 1 ] * _greenf,
					       val_b[ 2 ] * _bluef );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 2,
					       val_c[ 0 ] * _redf,
					       val_c[ 1 ] * _greenf,
					       val_c[ 2 ] * _bluef );
      cairo_mesh_pattern_end_patch (pattern);
      cairo_set_source( _cr, pattern );
      cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_close_path( _cr );
      cairo_fill( _cr );
      cairo_pattern_destroy( pattern );
    }

    void drawFlatTriangle( RealPoint a, RealPoint b, RealPoint c, 
			   Value val )
    {
      cairo_set_source_rgb( _cr,
			    val[ 0 ] * _redf,
			    val[ 1 ] * _greenf,
			    val[ 2 ] * _bluef );
      cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_close_path( _cr );
      cairo_fill( _cr );
    }

    void drawFlatLine( RealPoint a, RealPoint b,
		       Value val )
    {
      cairo_set_source_rgb( _cr,
			    val[ 0 ] * _redf,
			    val[ 1 ] * _greenf,
			    val[ 2 ] * _bluef );
      cairo_set_line_width( _cr, 1.0 ); 
      cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
      //cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_close_path( _cr );
      cairo_stroke( _cr );
      cairo_set_line_width( _cr, 0.0 ); 
    }

    
    template <typename BezierTriangle>
    void drawColorBezierTriangle( const BezierTriangle& BT )
    {
      // Scans all pixels in triangle
      RealPoint low = ij( BT.vertex(0) ).inf( ij( BT.vertex(1) ) ).inf( ij( BT.vertex(2) ) ); 
      RealPoint sup = ij( BT.vertex(0) ).sup( ij( BT.vertex(1) ) ).sup( ij( BT.vertex(2) ) ); 
      Point   ijlow = Point( (int) low[ 0 ], (int) low[ 1 ] );
      Point   ijsup = Point( (int) sup[ 0 ], (int) sup[ 1 ] );
      Domain domain( ijlow, ijsup );
      for ( Point p : domain ) {
	RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	auto bary_q = BT.barycentric( q );
	std::cout << "p=" << p << " q=" << q << " bary_q=" << bary_q << std::endl;
	if ( BT.isInTriangle( bary_q ) ) {
	  Value val = BT( bary_q );
	  std::cout << "   Inside val=" << val << std::endl;
          // cairo_set_source_rgb( _cr, 255.0 * _redf,
	  // 			0.0 * _greenf, 0.0 * _bluef );
          cairo_set_source_rgb( _cr, val[ 0 ] * _redf,
	  			val[ 1 ] * _greenf, val[ 2 ] * _bluef );
	  cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	  cairo_fill( _cr );
	}
      }
    }

    template <typename BezierTriangle>
    void drawMonochromeBezierTriangle( const BezierTriangle& BT, int channel )
    {
      Scalar   red = ( channel == 0 || channel == 3 ) ? 1.0 : 0.0;
      Scalar green = ( channel == 1 || channel == 3 ) ? 1.0 : 0.0;
      Scalar  blue = ( channel == 2 || channel == 3 ) ? 1.0 : 0.0;
      // Scans all pixels in triangle
      RealPoint low = ij( BT.b(0) ).inf( ij( BT.b(1) ) ).inf( ij( BT.b(2) ) ); 
      RealPoint sup = ij( BT.b(0) ).sup( ij( BT.b(1) ) ).sup( ij( BT.b(2) ) ); 
      Point   ijlow = Point( (int) low[ 0 ], (int) low[ 1 ] );
      Point   ijsup = Point( (int) sup[ 0 ], (int) sup[ 1 ] );
      Domain domain( ijlow, ijsup );
      for ( Point p : domain ) {
	RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	auto bary_q = BT.barycentric( q );
	std::cout << "p=" << p << " q=" << q << " bary_q=" << bary_q << std::endl;
	if ( BT.isInTriangle( bary_q ) ) {
	  Value val = BT( bary_q );
	  std::cout << "   Inside val=" << val << std::endl;
          cairo_set_source_rgb( _cr,
				red   * val[ 0 ] * _redf,
				green * val[ 0 ] * _greenf,
				blue  * val[ 0 ] * _bluef );
	  cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	  cairo_fill( _cr );
	}
      }
    }

    static inline
    Scalar det( const RealVector& v, const RealVector& w )
    {
      return v[ 0 ] * w[ 1 ] - v[ 1 ] * w[ 0 ];
    }
    
    /// @return the barycentric coordinates of point \a p in triangle (abc).
    static
    RealPoint3 barycentric( const RealPoint& p,
			    const RealPoint& a,
			    const RealPoint& b,
			    const RealPoint& c )
    {
      RealPoint3 bc = { det( b - p, c - p ),
			det( c - p, a - p ),
			det( a - p, b - p ) };
      return bc / ( bc[ 0 ] + bc[ 1 ] + bc[ 2 ] );
    }


    /// @param bc any barycentric coordinates (sum is 1).
    ///
    /// @return 'true' iff the given barycentric coordinates \a bc
    /// indicates a point inside or on the boundary of the triangle.
    static
    bool isInTriangle( const RealPoint3& bc )
    {
      return ( 0 <= bc[ 0 ] ) && ( bc[ 0 ] <= 1 )
	&&   ( 0 <= bc[ 1 ] ) && ( bc[ 1 ] <= 1 )
	&&   ( 0 <= bc[ 2 ] ) && ( bc[ 2 ] <= 1 );
    }

    /// bc[ 0 ] + bc[ 1 ] + bc[ 2 ] = 1
    static
    Scalar linearPOU( Scalar va, Scalar vb, Scalar vc, const RealPoint3& bc )
    {
      return bc[ 0 ] * va + bc[ 1 ] * vb + bc[ 2 ] * vc;
    }

    static
    Scalar sqr( Scalar x )
    { return x*x; }

    static
    Scalar gaussianPOU( Scalar va, Scalar vb, Scalar vc, const RealPoint3& bc )
    {
      Scalar ca = bc[ 0 ] > 0.0 ? exp( 1.0 - 1.0 / ( bc[ 0 ] ) ) : 0.0; 
      Scalar cb = bc[ 1 ] > 0.0 ? exp( 1.0 - 1.0 / ( bc[ 1 ] ) ) : 0.0; 
      Scalar cc = bc[ 2 ] > 0.0 ? exp( 1.0 - 1.0 / ( bc[ 2 ] ) ) : 0.0; 
      return ( ca * va + cb * vb + cc * vc ) / ( ca + cb + cc );
    }
    
    template <typename Fct3>
    void drawPartitionOfUnityTriangle
    ( RealPoint a, RealPoint b, RealPoint c,
      Fct3     fa, Fct3     fb, Fct3     fc )
    {
      if ( _color ) drawColorPartitionOfUnityTriangle( a, b, c, fa, fb, fc );
      else  drawGrayLevelPartitionOfUnityTriangle( a, b, c, fa, fb, fc );
    }
    
    template <typename Fct3>
    void drawGrayLevelPartitionOfUnityTriangle
    ( RealPoint a, RealPoint b, RealPoint c,
      Fct3     fa, Fct3     fb, Fct3     fc )
    {
      // Scans all pixels in triangle
      RealPoint low = ij( a ).inf( ij( b ) ).inf( ij( c ) ); 
      RealPoint sup = ij( a ).sup( ij( b ) ).sup( ij( c ) ); 
      Point   ijlow = Point( (int) low[ 0 ], (int) low[ 1 ] );
      Point   ijsup = Point( (int) sup[ 0 ], (int) sup[ 1 ] );
      Domain domain( ijlow, ijsup );
      for ( Point p : domain ) {
	const RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	const auto     bq = barycentric( q, a, b, c );
	if ( isInTriangle( bq ) ) {
	  Scalar v = gaussianPOU( fa[ 0 ]( q ), fb[ 0 ]( q ), fc[ 0 ]( q ), bq );
          cairo_set_source_rgb( _cr, v * _redf, v * _greenf, v * _bluef );
	  cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	  cairo_fill( _cr );
	}
      }
    }
    template <typename Fct3>
    void drawColorPartitionOfUnityTriangle
    ( RealPoint a, RealPoint b, RealPoint c,
      Fct3     fa, Fct3     fb, Fct3     fc )
    {
      // Scans all pixels in triangle
      RealPoint low = ij( a ).inf( ij( b ) ).inf( ij( c ) ); 
      RealPoint sup = ij( a ).sup( ij( b ) ).sup( ij( c ) ); 
      Point   ijlow = Point( (int) low[ 0 ], (int) low[ 1 ] );
      Point   ijsup = Point( (int) sup[ 0 ], (int) sup[ 1 ] );
      Domain domain( ijlow, ijsup );
      for ( Point p : domain ) {
	const RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	const auto     bq = barycentric( q, a, b, c );
	if ( isInTriangle( bq ) ) {
	  Scalar vr = gaussianPOU( fa[ 0 ]( q ), fb[ 0 ]( q ), fc[ 0 ]( q ), bq );
	  Scalar vg = gaussianPOU( fa[ 1 ]( q ), fb[ 1 ]( q ), fc[ 1 ]( q ), bq );
	  Scalar vb = gaussianPOU( fa[ 2 ]( q ), fb[ 2 ]( q ), fc[ 2 ]( q ), bq );
          cairo_set_source_rgb( _cr, vr * _redf, vg * _greenf, vb * _bluef );
	  cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	  cairo_fill( _cr );
	}
      }
    }
    template <typename EvalMetrics>
    void drawPartitionOfUnity( Domain    domain, const EvalMetrics& f )
    {
      if ( _color ) drawColorPartitionOfUnity( domain, f );
      else          drawGrayLevelPartitionOfUnity( domain, f );
    }
    
    template <typename EvalMetrics>
    void drawColorPartitionOfUnity( Domain    domain, const EvalMetrics&   f )
    {
      // Scans all pixels in domain
      std::cout << "Domain=" << domain << std::endl;
      std::vector<Point> D;
      for ( Point p : domain ) D.push_back( p );
      for ( Point p : D ) {
	const RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	//std::cout << "p=" << p << " q=" << q << std::endl;
	const Value     V = f.evalMetrics( q );
	cairo_set_source_rgb( _cr, V[0] * _redf, V[1] * _greenf, V[2] * _bluef );
	cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	cairo_fill( _cr );
      }
    }
    template <typename EvalMetrics>
    void drawGrayLevelPartitionOfUnity( Domain    domain, const EvalMetrics&   f )
    {
      // Scans all pixels in domain
      for ( Point p : domain ) {
	const RealPoint q = { x( p[ 0 ] ), y( p[ 1 ] ) };
	const auto      V = f.evalMetrics( q );
	cairo_set_source_rgb( _cr, V[0] * _redf, V[0] * _greenf, V[0] * _bluef );
	cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
	cairo_fill( _cr );
      }
    }

    void drawPixel( Point p, Value V )
    {
      if ( _color )
	cairo_set_source_rgb( _cr, V[0] * _redf, V[1] * _greenf, V[2] * _bluef );
      else
	cairo_set_source_rgb( _cr, V[0] * _redf, V[0] * _greenf, V[0] * _bluef );
      cairo_rectangle( _cr, p[ 0 ], p[ 1 ], 1, 1 );
      cairo_fill( _cr );
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class CairoViewer

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CairoViewer_h

#undef CairoViewer_RECURSES
#endif // else defined(CairoViewer_RECURSES)
