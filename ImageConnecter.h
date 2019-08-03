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
 * @file ImageConnecter.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module ImageConnecter.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ImageConnecter_RECURSES)
#error Recursive header files inclusion detected in ImageConnecter.h
#else // defined(ImageConnecter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImageConnecter_RECURSES

#if !defined ImageConnecter_h
/** Prevents repeated inclusion of headers. */
#define ImageConnecter_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/CImage.h>
#include <DGtal/io/writers/PPMWriter.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class ImageConnecter
  /**
     Description of class 'ImageConnecter' <p> \brief Aim: This class 
     optimizes the initial connectedness (4 or 8) of a given image.

     In this class, colors are represented as a 3-vector with
     components in 0..255.

     @tparam TImage the type of 2D image.
  */
  template <typename TImage>
  class ImageConnecter {
  public:
    typedef TImage                            Image;
    typedef ImageConnecter< TImage >          Self;
    BOOST_CONCEPT_ASSERT(( concepts::CImage< Image > ));

    typedef typename Image::Domain                Domain;
    typedef typename Image::Value                 Value;
    typedef typename Domain::Space                Space;
    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef typename Vector::Component            Integer;
    typedef typename RealVector::Component        Scalar;
    typedef std::function<Scalar( Value, Value )> Comparator; // 0 is equal

    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));
    
    struct Quad {
      Point base;
      Size  min;
      Size  max;
      bool operator<( const Quad& other ) const
      {
	// Gives better results.
	return ( min < other.min )
      	  || ( ( min == other.min )
      	       && ( max > other.max ) ); 
	// return ( max > other.max )
	// 	    || ( ( max == other.max )
	// 		 && ( min < other.min ) );
	// return ( ( max - min ) > ( other.max - other.min ) )
	//     || ( ( ( max - min ) == ( other.max - other.min ) )
	// 	 && ( min < other.min ) );
      }
    };
    
    // Union-find data structure.
    struct Element {
      Point    point;
      Size     rank;
      Element* father;
      Size     nb4;
      Size     nb8;
      int      order;
      Element( Element* self = nullptr, Point p = Point( 0, 0 ) )
	: point( p ), rank( 0 ), father( self ), nb4 ( 1 ), nb8( 1 ),
	  order( -1 ) {}
    };

    // return the root of the tree containing e
    static Element* find( Element* e )
    {
      if ( e != e->father )
	e->father = find( e->father );
      return e->father;
    }
    // link the two trees of roots x and y.
    static void link4( Element* x, Element* y )
    {
      if ( x->rank > y->rank ) {
	y->father = x;
	x->nb4   += y->nb4; 
	x->nb8   += y->nb8;
	x->order  = ( y->order == -1 ) ? x->order
	  : ( x->order == - 1 ) ? y->order
	  : std::min( x->order, y->order );
      }
      else {
	x->father = y;
	y->nb4   += x->nb4; 
	y->nb8   += x->nb8; 
	y->order  = ( y->order == -1 ) ? x->order
	  : ( x->order == - 1 ) ? y->order
	  : std::min( x->order, y->order );
	if ( x->rank == y->rank ) y->rank += 1;
      }
    }
    // link the two trees of roots x and y.
    static void link8( Element* x, Element* y )
    {
      if ( x->rank > y->rank ) {
	y->father = x;
	x->nb8   += y->nb8; 
      }
      else {
	x->father = y;
	y->nb8   += x->nb8; 
	if ( x->rank == y->rank ) y->rank += 1;
      }
    }
    // Given two elements x and y in different trees, makes the union
    // of the two trees.
    static void merge4( Element* x, Element* y )
    {
      link4( find( x ), find( y ) );
    }
    // Given two elements x and y in different trees, makes the union
    // of the two trees.
    static void merge8( Element* x, Element* y )
    {
      link8( find( x ), find( y ) );
    }

    /// The possible strategies for solving connectivity ambiguities.
    /// SizeOfConnectedComponents:  best for 1-bit images.
    /// OrderOfConnectedComponents: best in the general case.
    enum Strategy { SizeOfConnectedComponents,
		    OrderOfConnectedComponents };
    /// The three possible configurations between 4 points.
    enum DiagonalConfiguration { Default = 0, Diagonal00_11, Diagonal10_01 };

    /// Encodes connections between a group of four pixels sharing a pointel.
    struct QuadConfiguration {
      DiagonalConfiguration diagonal;   //< choice of diagonal edge
      bool                  horizontal; //< force horizontal edge
      bool                  vertical;   //< force vertical edge
      /// Default constructor.
      QuadConfiguration()
	: diagonal( Default ), horizontal( false ), vertical( false )
      {}
    };

    /// Stores how edges should be connected within a group of four
    /// pixels sharing a pointel.
    std::vector< QuadConfiguration > connections;

    /// Useful for union-find algorithm
    std::vector< Element >           labels;

    /// The width of the image
    Integer                          width;

    /// Tells is in debug mode.
    bool                             _debug;
  public:

    /// Constructor.
    ///
    /// @param debug if true, outputs images of connected components
    /// (cc.ppm) or orders (order.ppm)
    ImageConnecter( bool debug = false )
      : _debug( debug )
    {}

    /// Given a point \a p representing the base corner of a group of
    /// 4 pixels sharing a vertex, returns how edges should be
    /// connected.
    QuadConfiguration howConnected( Point p ) const
    {
      return connections[ p[ 1 ] * width + p[ 0 ] ];
    }

    /// Compute connectivities in the image \a I, according to the
    /// chosen strategy \a s and a comparator such that 'comp( I(p1),
    /// I(p2) ) <= same' indicates that the pixels \a p1 and \a p2
    /// have the same value.
    void init
    ( Strategy s, const Image& I, Comparator comp, Scalar same )
    {
      if ( s == SizeOfConnectedComponents )
	initWithSizeOfConnectedComponents ( I, comp, same );
      else
	initWithOrderOnConnectedComponents ( I, comp, same );
    }

    /// Compute connectivities using the respective of 4 and 8
    /// connected components.
    void initWithSizeOfConnectedComponents
    ( const Image& I, Comparator comp, Scalar same )
    {
      // Creates 4-connected components
      const Domain& domain = I.domain();
      const Domain rdomain ( domain.lowerBound(),
			     domain.upperBound() - Point::diagonal( 1 ) );
      const Vector& extent = I.extent();
      width  = extent[ 0 ];
      labels.resize( domain.size() );
      connections.resize( domain.size() );
      Size i = 0;
      for ( auto p : domain ) {
	Element& self = labels[ i++ ];
	self = Element( &self, p );
      }
      // Scan and merge
      Point up = domain.upperBound();
      for ( auto p : domain ) {
	Value v = I( p );
	if ( p[ 0 ] < up[ 0 ] ) {
	  Point q = p + Vector( 1, 0 );
	  if ( comp( v, I( q ) ) <= same ) {
	    Element* e1 = & labels[ p[ 1 ] * width + p[ 0 ] ];
	    Element* e2 = & labels[ q[ 1 ] * width + q[ 0 ] ];
	    if ( find( e1 ) != find( e2 ) ) {
	      merge4( e1, e2 );
	      // std::cout << "merge " << e1->nb << " " << e2->nb << std::endl;
	    }
	  }
	}
	if ( p[ 1 ] < up[ 1 ] ) {
	  Point q = p + Vector( 0, 1 );
	  if ( comp( v, I( q ) ) <= same ) {
	    Element* e1 = & labels[ p[ 1 ] * width + p[ 0 ] ];
	    Element* e2 = & labels[ q[ 1 ] * width + q[ 0 ] ];
	    if ( find( e1 ) != find( e2 ) ) {
	      merge4( e1, e2 );
	      // std::cout << "merge " << e1->nb << " " << e2->nb << std::endl;
	      //std::cout << "merge " << p << " " << q << std::endl;
	    }
	  }
	}
      } // end scan for 4-connectedness
      // Connect diagonals around big regions.
      std::vector<Quad> quads( rdomain.size() );
      for ( auto p : rdomain ) {
	const Point p00 = p;
	const Point p10 = p + Vector( 1, 0 );
	const Point p01 = p + Vector( 0, 1 );
	const Point p11 = p + Vector( 1, 1 );
	Element*    e00 = find( & labels[ p00[ 1 ] * width + p00[ 0 ] ] );
	Element*    e10 = find( & labels[ p10[ 1 ] * width + p10[ 0 ] ] );
	Element*    e01 = find( & labels[ p01[ 1 ] * width + p01[ 0 ] ] );
	Element*    e11 = find( & labels[ p11[ 1 ] * width + p11[ 0 ] ] );
	Quad          q = { p00,
			    std::min( std::min( e00->nb4, e10->nb4 ),
				      std::min( e01->nb4, e11->nb4 ) ),
			    std::max( std::max( e00->nb4, e10->nb4 ),
				      std::max( e01->nb4, e11->nb4 ) ) };
	quads.push_back( q );
      }
      std::sort( quads.begin(), quads.end() );

      if ( _debug ) {
	typedef ImageContainerBySTLVector< Domain, Color > ColorImage;
	ColorImage debug( domain );
	std::vector< Color > cc_color( labels.size() );
	for ( int i = 0; i < cc_color.size(); ++i )
	  cc_color[ i ] = Color( (unsigned int) rand() % 0xffffff );
	for ( int i = 0; i < cc_color.size(); ++i ) {
	  Element* e = & labels[ i ];
	  if ( e != find( e ) ) 
	    cc_color[ i ] = cc_color[ find( e ) - & labels[ 0 ] ];
	}
	for ( auto p : domain )
	  debug.setValue( p, cc_color[ p[ 1 ] * width + p[ 0 ] ] );
	PPMWriter<ColorImage, std::function< Color(Color) > >
	  ::exportPPM("cc.ppm", debug, [] (Color c) { return c; } );
      }
	
      // Connect diagonals around big regions.
      Size nb00_11 = 0;
      Size nb10_01 = 0;
      for ( auto q : quads ) {
	//	std::cout << q.base << " " << q.min << " " << q.max << std::endl;
	const Point p00 = q.base;
	const Point p10 = p00 + Vector( 1, 0 );
	const Point p01 = p00 + Vector( 0, 1 );
	const Point p11 = p00 + Vector( 1, 1 );
	Element*    e00 = find( & labels[ p00[ 1 ] * width + p00[ 0 ] ] );
	Element*    e10 = find( & labels[ p10[ 1 ] * width + p10[ 0 ] ] );
	Element*    e01 = find( & labels[ p01[ 1 ] * width + p01[ 0 ] ] );
	Element*    e11 = find( & labels[ p11[ 1 ] * width + p11[ 0 ] ] );
	const Value v00 = I( p00 );
	const Value v10 = I( p10 );
	const Value v01 = I( p01 );
	const Value v11 = I( p11 );
	Scalar   s00_11 = comp( v00, v11 );
	Scalar   s10_01 = comp( v10, v01 );
	
	QuadConfiguration c;
	if ( comp( v00, v10 ) <= same ) c.horizontal = true;
	if ( comp( v00, v01 ) <= same ) c.vertical   = true;
	if ( ( same   <  s00_11 ) && ( same < s10_01 ) )
	  c.diagonal = Default;
	else if ( ( s00_11 <= same   ) && ( same < s10_01 ) )
	  c.diagonal = Diagonal00_11;
	else if ( ( s10_01 <= same ) && ( same < s00_11 ) )
	  c.diagonal = Diagonal10_01;
	// simpler solutions like this are better on 1-bit pixel art
	// than complex ones, which introduces weird dithering.
	else if ( ( e00 == e11 ) && ( e10 != e01 ) ) c.diagonal = Diagonal10_01;
	else if ( ( e10 == e01 ) && ( e00 != e11 ) ) c.diagonal = Diagonal00_11;
	else c.diagonal = Diagonal00_11;
	// else if ( ( e00 == e11 ) && ( e10 != e01 ) ) {
	//   if ( ( ( e00->nb4 + e11->nb4 ) < ( e10->nb4 + e01->nb4 ) ) )
	//     //&&( ( e00->nb8 + e11->nb8 ) >= 2 * ( e00->nb4 + e11->nb4 ) ) )
	//     c.diagonal = Diagonal00_11;
	//   else c.diagonal = Diagonal10_01;
	// }
	// else if ( ( e10 == e01 ) && ( e00 != e11 ) ) {
	//   if ( ( ( e10->nb4 + e01->nb4 ) < ( e00->nb4 + e11->nb4 ) ) )
	//     // &&( ( e10->nb8 + e01->nb8 ) >= 2 * ( e10->nb4 + e01->nb4 ) ) )
	//     c.diagonal = Diagonal10_01;
	//   else c.diagonal = Diagonal00_11;
	// }
	// else c.diagonal = Diagonal00_11;
	// else { 
	//   int  r00_11 = ( e00->nb8 + e11->nb8 ) / ( e00->nb4 + e11->nb4 );
	//   int  r10_01 = ( e10->nb8 + e01->nb8 ) / ( e10->nb4 + e01->nb4 );
	//   if ( r00_11 < r10_01 ) c.diagonal = Diagonal10_01;
	//   else if ( r00_11 > r10_01 ) c.diagonal = Diagonal00_11;
	//   else c.diagonal = Diagonal00_11;
	// }
	if ( c.diagonal == Diagonal00_11 ) {
	  if ( e00 != e11 ) merge8( e00, e11 );
	  nb00_11++;
	} else if ( c.diagonal == Diagonal10_01 ) {
	  if ( e10 != e01 ) merge8( e10, e01 );
	  nb10_01++;
	}

	connections[ p00[ 1 ] * width + p00[ 0 ] ] = c;
      }
      trace.info() << "nb00_11=" << nb00_11
		   << " nb10_01=" << nb10_01;
    }

    void fillRegionWithOrder( const Image& I, Comparator comp, Scalar same,
			      Point p, int& order,
			      std::set<Point>& M )
    {
      const Point lo = I.domain().lowerBound();
      const Point up = I.domain().upperBound();
      std::queue<Point> Q;
      Q.push( p );
      while ( ! Q.empty() ) {
	Point    pp = Q.front(); Q.pop();
	if ( M.find( pp ) != M.end() ) continue;
	M.insert( pp );
	Value     v = I( pp );
	Element* e1 = & labels[ linearize( pp ) ];
	e1->order   = ( e1->order == -1 )
	  ? order : std::min( e1->order, order );
	Point T[ 4 ] = { pp + Vector(  1,  0 ),
			 pp + Vector( -1,  0 ),
			 pp + Vector(  0,  1 ),
			 pp + Vector(  0, -1 ) };
	for ( int i = 0; i < 4; ++i ) {
	  Point q = T[ i ];
	  if ( ( q.inf( lo ) == lo ) && ( q.sup( up ) == up ) ) {
	    Element* e2 = & labels[ linearize( q ) ];
	    if ( comp( v, I( q ) ) <= same )
	      Q.push( q );
	    else
	      e2->order   = ( e2->order == -1 )
		? order+1 : std::min( e2->order, order+1 );
	  }
	} // 	for ( int i = 0; i < 4; ++i ) {
      }	 // while (! Q.empty() ) {
      ++order;
    }
    
    void processPoint( const Image& I, Comparator comp, Scalar same,
		       Point p, Size& order )
    {
      const Point lo = I.domain().lowerBound();
      const Point up = I.domain().upperBound();
      Value     v = I( p );
      Element* e1 = & labels[ linearize( p ) ];
      if ( p[ 0 ] < up[ 0 ] ) {
	Point q = p + Vector( 1, 0 );
	Element* e2 = & labels[ linearize( q ) ];
	if ( comp( v, I( q ) ) <= same ) {
	  if ( find( e1 ) != find( e2 ) ) {
	    merge4( e1, e2 );
	  }
	} else {
	  Element* fe2 = find( e2 );
	  if ( fe2->order == -1 ) fe2->order = ++order;
	}
      }
      if ( lo[ 0 ] < p[ 0 ] ) {
	Point q = p + Vector( -1, 0 );
	Element* e2 = & labels[ linearize( q ) ];
	if ( comp( v, I( q ) ) <= same ) {
	  if ( find( e1 ) != find( e2 ) ) {
	    merge4( e1, e2 );
	  }
	} else {
	  Element* fe2 = find( e2 );
	  if ( fe2->order == -1 ) fe2->order = ++order;
	}
      }
      if ( p[ 1 ] < up[ 1 ] ) {
	Point q = p + Vector( 0, 1 );
	Element* e2 = & labels[ linearize( q ) ];
	if ( comp( v, I( q ) ) <= same ) {
	  if ( find( e1 ) != find( e2 ) ) {
	    merge4( e1, e2 );
	  }
	} else {
	  Element* fe2 = find( e2 );
	  if ( fe2->order == -1 ) fe2->order = ++order;
	}
      }
      if ( lo[ 1 ] < p[ 1 ] ) {
	Point q = p + Vector( 0, -1 );
	Element* e2 = & labels[ linearize( q ) ];
	if ( comp( v, I( q ) ) <= same ) {
	  if ( find( e1 ) != find( e2 ) ) {
	    merge4( e1, e2 );
	  }
	} else {
	  Element* fe2 = find( e2 );
	  if ( fe2->order == -1 ) fe2->order = ++order;
	}
      }
    }
    // The strategy is to label the components from the outside toward
    // the inside, in order to get an order on regions.
    void initWithOrderOnConnectedComponents
    ( const Image& I, Comparator comp, Scalar same )
    {
      // Creates 4-connected components
      const Domain& domain = I.domain();
      const Domain rdomain ( domain.lowerBound(),
			     domain.upperBound() - Point::diagonal( 1 ) );
      const Vector& extent = I.extent();
      width  = extent[ 0 ];
      // Initializes all regions as singletons.
      labels.resize( domain.size() );
      connections.resize( domain.size() );
      Size i = 0;
      for ( auto p : domain ) {
	Element& self = labels[ i++ ];
	self = Element( &self, p );
      }
      // We scan the domain from the outside toward the inside
      const Point lo = domain.lowerBound();
      const Point up = domain.upperBound();
      Point      clo = lo;
      Point      cup = up;
      int      order = 0;
      std::set<Point> M;
      while ( ( clo[ 0 ] <= cup[ 0 ] ) &&
	      ( clo[ 1 ] <= cup[ 1 ] ) ) {
	for ( int x = clo[ 0 ]; x <= cup[ 0 ]; ++x ) {
	  Point p1( x, clo[ 1 ] );
	  if ( M.find( p1 ) == M.end() )
	    fillRegionWithOrder( I, comp, same, p1, order, M );
	  Point p2( x, cup[ 1 ] );
	  if ( M.find( p2 ) == M.end() )
	    fillRegionWithOrder( I, comp, same, p2, order, M );
	}
	for ( int y = clo[ 1 ]+1; y < cup[ 1 ]; ++y ) {
	  Point p1( clo[ 0 ], y );
	  if ( M.find( p1 ) == M.end() )
	    fillRegionWithOrder( I, comp, same, p1, order, M );
	  Point p2( cup[ 0 ], y );
	  if ( M.find( p2 ) == M.end() )
	    fillRegionWithOrder( I, comp, same, p2, order, M );
	}
	clo += Point::diagonal( 1 );
	cup -= Point::diagonal( 1 );
      }
      // end scan for 4-connectedness

      if ( _debug ) {
	typedef ImageContainerBySTLVector< Domain, Color > ColorImage;
	ColorImage debug( domain );
	std::vector< Color > cc_color( labels.size() );
	for ( int i = 0; i < cc_color.size(); ++i )
	  cc_color[ i ] = Color( (unsigned int) rand() % 0xffffff );
	for ( int i = 0; i < cc_color.size(); ++i ) {
	  Element* e = & labels[ i ];
	  if ( e != find( e ) ) 
	    cc_color[ i ] = cc_color[ find( e ) - & labels[ 0 ] ];
	}
	for ( auto p : domain )
	  debug.setValue( p, cc_color[ p[ 1 ] * width + p[ 0 ] ] );
	PPMWriter<ColorImage, std::function< Color(Color) > >
	  ::exportPPM("cc.ppm", debug, [] (Color c) { return c; } );
	for ( auto p : domain ) {
	  int order = find( & labels[ linearize( p ) ] )->order;
	  debug.setValue( p, Color( ( ( order % 16 ) * 16 ) % 256, ( order / 16 ) % 256, order % 256 ) );
	}
	PPMWriter<ColorImage, std::function< Color(Color) > >
	  ::exportPPM("order.ppm", debug, [] (Color c) { return c; } );
      }
      
      // Connect diagonals around big regions.
      Size nb00_11 = 0;
      Size nb10_01 = 0;
      for ( auto p : rdomain ) {
	const Point p00 = p;
	const Point p10 = p + Vector( 1, 0 );
	const Point p01 = p + Vector( 0, 1 );
	const Point p11 = p + Vector( 1, 1 );
	Element*    e00 = find( & labels[ linearize( p00 ) ] );
	Element*    e10 = find( & labels[ linearize( p10 ) ] );
	Element*    e01 = find( & labels[ linearize( p01 ) ] );
	Element*    e11 = find( & labels[ linearize( p11 ) ] );
	const Value v00 = I( p00 );
	const Value v10 = I( p10 );
	const Value v01 = I( p01 );
	const Value v11 = I( p11 );
	Scalar   s00_11 = comp( v00, v11 );
	Scalar   s10_01 = comp( v10, v01 );
	QuadConfiguration c;
	if ( comp( v00, v10 ) <= same ) c.horizontal = true;
	if ( comp( v00, v01 ) <= same ) c.vertical   = true;
	if ( ( same   <  s00_11 ) && ( same < s10_01 ) )
	  c.diagonal = Default;
	else if ( ( s00_11 <= same   ) && ( same < s10_01 ) )
	  c.diagonal = Diagonal00_11;
	else if ( ( s10_01 <= same ) && ( same < s00_11 ) )
	  c.diagonal = Diagonal10_01;
	else {
	  Size m00_11 = std::min( e00->order, e11->order ); 
	  Size m10_01 = std::min( e10->order, e01->order );
	  if ( m00_11 < m10_01 ) c.diagonal = Diagonal10_01;
	  else c.diagonal = Diagonal00_11;
	}
	if ( c.diagonal == Diagonal00_11 ) {
	  //if ( e00 != e11 ) merge8( e00, e11 );
	  nb00_11++;
	} else if ( c.diagonal == Diagonal10_01 ) {
	  //if ( e10 != e01 ) merge8( e10, e01 );
	  nb10_01++;
	}
	connections[ p00[ 1 ] * width + p00[ 0 ] ] = c;
      }
      trace.info() << "nb00_11=" << nb00_11
		   << " nb10_01=" << nb10_01;
    }
    
    /// Destructor.
    ~ImageConnecter()
    {
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:
    /// @param p any point
    /// @return the corresponding index in the 1-dimensional array.
    Size linearize( Point p ) const
    {
      return p[ 1 ] * width + p[ 0 ];
    }

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ImageConnecter

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ImageConnecter_h

#undef ImageConnecter_RECURSES
#endif // else defined(ImageConnecter_RECURSES)
