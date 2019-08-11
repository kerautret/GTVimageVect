#include <cmath>
#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/images/CImage.h>
#include <DGtal/io/Color.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/io/writers/PGMWriter.h>
//#include <DGtal/math/RealFFT.h>

namespace DGtal {

  namespace image {
    struct ColorToRedFunctor {
      int operator()( int color ) const
      { return (color >> 16) & 0xff; }
    };
    struct ColorToGreenFunctor {
      int operator()( int color ) const
      { return (color >> 8) & 0xff; }
    };
    struct ColorToBlueFunctor {
      int operator()( int color ) const
      { return color & 0xff; }
    };
    struct ColorToGrayFunctor {
      int operator()( int color ) const
      {
	double val = 0.2126 * ( (double) ( ( color >> 16 ) & 0xff ) )
	  +          0.7152 * ( (double) ( ( color >> 8  ) & 0xff ) )
	  +          0.0722 * ( (double) ( ( color >> 0  ) & 0xff ) );
	return ( (int) round( val ) ) & 0xff;
      }
    };
    struct DoubleToGrayFunctor {
      unsigned char operator()( double scalar ) const
      { return (unsigned char) std::min( 255.0, std::max( 0.0, round( scalar ) ) ); }
    };
    struct SignedDoubleToGrayFunctor {
      unsigned char operator()( double scalar ) const
      { return (unsigned char) std::min( 255.0, std::max( 0.0, 128.0 + round( scalar ) ) ); }
    };
    struct GrayToGrayFunctor {
      int operator()( int color ) const
      { return color & 0xff; }
    };
    struct GrayToRedGreen {
      DGtal::Color operator()( int value ) const
      { 
	if ( value >= 128 ) return DGtal::Color( 0, (value-128)*2+1, 0 );
      else                return DGtal::Color( 0, 0, 255-value*2 );
      }
    };

    struct Zoom {
      typedef Z2i::Integer   Integer;
      typedef Z2i::Point     Point;
      typedef Z2i::RealPoint RealPoint;
      typedef Z2i::Domain    Domain;
      typedef RealPoint::Component Scalar;

      Zoom( Domain aDomain, int aZoom )
	: _domain( aDomain ), _zoom( aZoom ) {
	_zDomain = getZoomedDomain();
	_stride  =  _domain.upperBound()[ 0 ] -  _domain.lowerBound()[ 0 ] + 1;
	_zStride = _zDomain.upperBound()[ 0 ] - _zDomain.lowerBound()[ 0 ] + 1;
      }

      const Domain& domain() const { return _domain; }
      const Domain& zoomedDomain() const { return _zDomain; }
      Domain getZoomedDomain() const {
	return Domain( zoom( _domain.lowerBound() ), zoom( _domain.upperBound() ) );
      }

      Integer index( Point p ) const
      {
	return p[ 1 ] * _stride + p[ 0 ];
      }

      Integer zoomedIndex( Point p ) const
      {
	return p[ 1 ] * _zStride + p[ 0 ];
      }
      
      Point zoom( Point p ) const
      {
	return p * _zoom;
      }
      
      RealPoint zoom( RealPoint p ) const
      {
	return p * ( Scalar ) _zoom;
      }

      Point unzoom( Point p ) const
      {
	return p / _zoom;
      }

      RealPoint real( Point p ) const
      { return RealPoint( p[ 0 ], p[ 1 ] ); }

      Point closest( RealPoint p ) const
      { return Point( (Integer) round( p[ 0 ] ), (Integer) round( p[ 1 ] ) ); }
      
      RealPoint project( Point p ) const
      {
	return real( p ) / (Scalar) _zoom;
      }
      
      RealPoint unzoom( RealPoint p ) const
      {
	return p / ( Scalar ) _zoom;
      }
      
      Domain  _domain;
      Domain  _zDomain;
      Integer _zoom;
      Integer _stride;
      Integer _zStride;
    };
    
    struct Utils {
      typedef double                     Scalar;
      typedef PointVector< 3, Scalar >   Value;
      // typedef std::array< Value, 2 >     VectorValue;
      struct VectorValue {
	Value x;
	Value y;
      };
      typedef std::vector<Scalar>        ScalarForm;
      typedef std::vector<Value>         ValueForm;
      typedef std::vector<VectorValue>   VectorValueForm;
      
      template <typename Image>
      ScalarForm getScalarForm( const Image& image, bool color )
      {
	unsigned int v = 0;
	ScalarForm I( image.size() );
	if ( color ) {
	  auto F = ColorToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = (Scalar) F( val );
	} else {
	  auto F = GrayToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = (Scalar) F( val );
	}
	return I;
      }
      template <typename Image>
      ValueForm getValueForm( const Image& image, bool color )
      {
	unsigned int v = 0;
	ValueForm I( image.size() );
	if ( color ) {
	  auto R = ColorToRedFunctor();
	  auto G = ColorToGreenFunctor();
	  auto B = ColorToBlueFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = Value( (Scalar) R( val ), (Scalar) G( val ), (Scalar) B( val ) );
	} else {
	  auto F = GrayToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = Value( (Scalar) F( val ), (Scalar) F( val ), (Scalar) F( val ) );
	}
	return I;
      }
    };

    // ----------------------------------------------------------------------
    namespace functions {
      template <typename TImageDst, typename TImageSrc,
		typename TFunctor>
      void transform( TImageDst& dst, const TImageSrc& src,
		      TFunctor F )
      {
	std::transform( src.begin(), src.end(), dst.begin(), F );
      }

      template <typename ScalarImage>
      bool exportPGM( const char* name, const ScalarImage& image )
      {
	return PGMWriter<ScalarImage,DoubleToGrayFunctor>::exportPGM( name, image );
      }
      template <typename ScalarImage>
      bool exportGradPGM( const char* name, const ScalarImage& image )
      {
	return PGMWriter<ScalarImage,SignedDoubleToGrayFunctor>::exportPGM( name, image );
      }
      
      // template <typename TImage>
      // void bandFilter( TImage& image,
      //   	       double low_freq = 0.0, double high_freq = 0.5,
      //   	       double shift = 0.0,
      //   	       double output_min = 0.0,
      //   	       double output_max = 255.0 )
      // {
      //   typedef TImage Image;
      //   BOOST_CONCEPT_ASSERT(( DGtal::concepts::CImage< Image > ));

      //   typedef typename Image::Domain         Domain;
      //   typedef typename Domain::Space         Space;
      //   typedef typename Space::Point          Point;
      //   typedef typename Space::RealVector     RealVector;
      //   typedef typename RealVector::Component Scalar;
      //   typedef typename Image::Value          Value;

      //   Domain domain( image.domain() );
      //   RealFFT< Domain, Scalar > F ( domain );
      //   auto IF = F.getSpatialImage();
      //   auto FF = F.getFreqImage();
      //   // trace.info() << "spatial   domain = " << F.getSpatialDomain() << std::endl;
      //   // trace.info() << "frequency domain = " << F.getFreqDomain() << std::endl;
      //   auto IF_it = IF.begin();
      //   for ( auto I_it  = image.begin(), I_itE = image.end(); I_it != I_itE; I_it++ )
      //     *IF_it++ = (Scalar) *I_it;
      //   F.forwardFFT();
      //   Scalar l2 = (Scalar) ( low_freq  * low_freq  );
      //   Scalar h2 = (Scalar) ( high_freq * high_freq );
      //   for ( auto&& fp : F.getFreqDomain() )
      //     {
      //       auto  freq = F.calcScaledFreqCoords( fp );
      //       auto nfreq = freq.dot( freq );
      //       if ( nfreq < l2 || nfreq > h2 )
      //         FF.setValue( fp, (Scalar) 0 );
      //     }
      //   F.backwardFFT();
      //   auto I_it  = image.begin();
      //   for ( auto IF_it = IF.begin(), IF_itE = IF.end(); IF_it != IF_itE; IF_it++ ) {
      //     Value val = (Value) std::max( output_min,
      //   				std::min( output_max, shift + (*IF_it) ) );
      //     *I_it++ = val;
      //   }
      // }      

      // template <typename ScalarImage>
      // void poisson( ScalarImage& I, const ScalarImage& Gx, const ScalarImage& Gy )
      // {
      //   typedef typename ScalarImage::Point Point;
      //   typedef typename ScalarImage::Domain Domain;
      //   typedef typename ScalarImage::Value  Scalar;
      //   Domain domain( I.domain() );
      //   Domain sym_domain( Domain( -domain.upperBound() - Point::diagonal( 1 ),
      //   			   domain.upperBound() ) );
      //   RealFFT< Domain, double > U ( sym_domain );
      //   // Calcule l'image gradient.
      //   RealFFT< Domain, double > Vx( sym_domain );
      //   RealFFT< Domain, double > Vy( sym_domain );
      //   auto IU  = U.getSpatialImage();
      //   auto IVx = Vx.getSpatialImage();
      //   auto IVy = Vy.getSpatialImage();
      //   for ( auto&& p : domain ) 
      //     {
      //       // Extend by symetry
      //       Point sym_px ( -1 - p[ 0 ], p[ 1 ] );
      //       Point sym_py ( p[ 0 ], -1 - p[ 1 ] );
      //       Point sym_pxy( -1 - p[ 0 ], -1 - p[ 1 ] );
      //       std::array<Point,4> all_p = { p, sym_px, sym_py, sym_pxy };
      //       double  vI = (double)  I( p );
      //       double vGx = (double) Gx( p );
      //       double vGy = (double) Gy( p );
      //       for ( auto&&q : all_p ) IU.setValue(  q, vI );
      //       IVx.setValue(       p,  vGx );
      //       IVx.setValue(  sym_px, -vGx );
      //       IVx.setValue(  sym_py,  vGx );
      //       IVx.setValue( sym_pxy, -vGx );
      //       IVy.setValue(       p,  vGy );
      //       IVy.setValue(  sym_px,  vGy );
      //       IVy.setValue(  sym_py, -vGy );
      //       IVy.setValue( sym_pxy, -vGy );
      //     }
      //   trace.info() << "Compute FFT[V]" << std::endl;
      //   trace.info() << "spatial   domain = " << Vx.getSpatialDomain() << std::endl;
      //   trace.info() << "frequency domain = " << Vx.getFreqDomain() << std::endl;
      //   U.forwardFFT();
      //   Vx.forwardFFT();
      //   Vy.forwardFFT();
      //   Domain freq_domain = Vx.getFreqDomain();
      //   auto FU  = U.getFreqImage();
      //   auto FVx = Vx.getFreqImage();
      //   auto FVy = Vy.getFreqImage();
      //   typedef typename RealFFT< Domain, double >::Complex Complex;
      //   const double two_pi = 2.0 * M_PI;
      //   const Complex     i = Complex( 0.0, 1.0 );
      //   const double      J = ( sym_domain.upperBound()[0] - sym_domain.lowerBound()[0] ) + 1;
      //   const double      L = ( sym_domain.upperBound()[1] - sym_domain.lowerBound()[1] ) + 1;
      //   for ( auto&& p : freq_domain )
      //     {
      //       auto  freq = U.calcScaledFreqCoords( p );
      //       const double  mj = freq[ 0 ];
      //       const double  nl = freq[ 1 ];
      //       Complex val = i * ( mj * FVx( p ) + nl * FVy( p ) )
      //         / ( two_pi * ( mj * mj + nl * nl ) );
      //       if ( p == Point::zero ) continue;
      //       else FU.setValue( p, -val );
      //     }
      //   U.backwardFFT();
      //   for ( auto&& p : domain ) 
      //     I.setValue( p, (Scalar) std::max( 0.0, std::min( 255.0, IU( p ) ) ) );
      // } // poisson
    } // namespace functions
  } // namespace image
} // namespace DGtal
