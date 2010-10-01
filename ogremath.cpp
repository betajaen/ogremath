


#include "ogremath.h"

namespace OgreMath
{
 namespace Private
 {
  
  const Real Math::POS_INFINITY = std::numeric_limits<Real>::infinity();
  const Real Math::NEG_INFINITY = -std::numeric_limits<Real>::infinity();
  const Real Math::PI = Real( 4.0 * atan( 1.0 ) );
  const Real Math::TWO_PI = Real( 2.0 * PI );
  const Real Math::HALF_PI = Real( 0.5 * PI );
  const Real Math::fDeg2Rad = PI / Real(180.0);
  const Real Math::fRad2Deg = Real(180.0) / PI;
  const Real Math::LOG2 = log(Real(2.0));
  
  int Math::mTrigTableSize;
  Math::AngleUnit Math::msAngleUnit;
  
  Real  Math::mTrigTableFactor;
  Real *Math::mSinTable = NULL;
  Real *Math::mTanTable = NULL;
  

 }
}