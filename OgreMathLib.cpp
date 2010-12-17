/**

    This file is part of OgreMathLib.

    Copyright (c) 2000-2009 Torus Knot Software Ltd
    Copyright (c) 2010 Robin Southern

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files(the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.

    Parts of this file may also contain code from:

     Wild Magic 0.2 Matrix math (free source code), http://www.geometrictools.com

*/

#include "OgreMathLib.h"

#define OGREMATHLIB_PLATFORM_WIN32 1
#define OGREMATHLIB_PLATFORM_LINUX 2
#define OGREMATHLIB_PLATFORM_APPLE 3
#define OGREMATHLIB_PLATFORM_SYMBIAN 4
#define OGREMATHLIB_PLATFORM_IPHONE 5
#define OGREMATHLIB_PLATFORM_ANDROID 6
#define OGREMATHLIB_PLATFORM_TEGRA2 7

#define OGREMATHLIB_COMPILER_MSVC 1
#define OGREMATHLIB_COMPILER_GNUC 2
#define OGREMATHLIB_COMPILER_BORL 3
#define OGREMATHLIB_COMPILER_WINSCW 4
#define OGREMATHLIB_COMPILER_GCCE 5

#define OGREMATHLIB_ENDIAN_LITTLE 1
#define OGREMATHLIB_ENDIAN_BIG 2

#define OGREMATHLIB_ARCHITECTURE_32 1
#define OGREMATHLIB_ARCHITECTURE_64 2

/* Finds the compiler type and version.
*/
#if defined( __GCCE__ )
#   define OGREMATHLIB_COMPILER OGREMATHLIB_COMPILER_GCCE
#   define OGREMATHLIB_COMP_VER _MSC_VER
//#	include <staticlibinit_gcce.h> // This is a GCCE toolchain workaround needed when compiling with GCCE
#elif defined( __WINSCW__ )
#   define OGREMATHLIB_COMPILER OGREMATHLIB_COMPILER_WINSCW
#   define OGREMATHLIB_COMP_VER _MSC_VER
#elif defined( _MSC_VER )
#   define OGREMATHLIB_COMPILER OGREMATHLIB_COMPILER_MSVC
#   define OGREMATHLIB_COMP_VER _MSC_VER
#elif defined( __GNUC__ )
#   define OGREMATHLIB_COMPILER OGREMATHLIB_COMPILER_GNUC
#   define OGREMATHLIB_COMP_VER (((__GNUC__)*100) + \
        (__GNUC_MINOR__*10) + \
        __GNUC_PATCHLEVEL__)

#elif defined( __BORLANDC__ )
#   define OGREMATHLIB_COMPILER OGREMATHLIB_COMPILER_BORL
#   define OGREMATHLIB_COMP_VER __BCPLUSPLUS__
#   define __FUNCTION__ __FUNC__
#else
#   pragma error "No known compiler. Abort! Abort!"

#endif

/* See if we can use __forceinline or if we need to use __inline instead */
#if OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC
#   if OGREMATHLIB_COMP_VER >= 1200
#       define OGREMATHLIB_FORCEINLINE __forceinline
#   endif
#elif defined(__MINGW32__)
#   if !defined(OGREMATHLIB_FORCEINLINE)
#       define OGREMATHLIB_FORCEINLINE __inline
#   endif
#else
#   define OGREMATHLIB_FORCEINLINE __inline
#endif

/* Finds the current platform */

#if defined( __SYMBIAN32__ )
#   define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_SYMBIAN
#elif defined( __WIN32__ ) || defined( _WIN32 )
#   define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_WIN32
#elif defined( __APPLE_CC__)
// Device                                                     Simulator
// Both requiring OS version 3.0 or greater
#   if __ENVIRONMENT_IPHONE_OS_VERSION_MIN_REQUIRED__ >= 30000 || __IPHONE_OS_VERSION_MIN_REQUIRED >= 30000
#       define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_IPHONE
#   else
#       define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_APPLE
#   endif
#elif defined(linux) && defined(__arm__)
// TODO: This is NOT the correct way to detect the Tegra 2 platform but it works for now.
// It doesn't appear that GCC defines any platform specific macros.
#   define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_TEGRA2
#elif defined(__ANDROID__)
#	define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_ANDROID
#else
#   define OGREMATHLIB_PLATFORM OGREMATHLIB_PLATFORM_LINUX
#endif

/* Find the arch type */
#if defined(__x86_64__) || defined(_M_X64) || defined(__powerpc64__) || defined(__alpha__) || defined(__ia64__) || defined(__s390__) || defined(__s390x__)
#   define OGREMATHLIB_ARCH_TYPE OGREMATHLIB_ARCHITECTURE_64
#else
#   define OGREMATHLIB_ARCH_TYPE OGREMATHLIB_ARCHITECTURE_32
#endif


#define OGRELIBMATH_ALLOC_T(TYPE, SIZE, CATEGORY) (TYPE*) malloc(sizeof(TYPE) * SIZE)
#define OGREMATHLIB_FREE(MEM, CATEGORY) free(MEM)
#if OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC
#   if OGREMATHLIB_COMP_VER >= 1200
#       define OGREMATHLIB_FORCEINLINE __forceinline
#   endif
#elif defined(__MINGW32__)
#   if !defined(OGREMATHLIB_FORCEINLINE)
#       define OGREMATHLIB_FORCEINLINE __inline
#   endif
#else
#   define OGREMATHLIB_FORCEINLINE __inline
#endif

namespace OgreMathLib
{

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC
#  pragma warning (push)
// disable "instruction may be inaccurate on some Pentiums"
#  pragma warning (disable : 4725)
#endif

/*=============================================================================
 ASM math routines posted by davepermen et al on flipcode forums
=============================================================================*/
const float pi = 4.0f * atan( 1.0f );
const float half_pi = 0.5f * pi;

/*=============================================================================
	NO EXPLICIT RETURN REQUIRED FROM THESE METHODS!!
=============================================================================*/
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC && OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
#	pragma warning( push )
#	pragma warning( disable: 4035 )
#endif

float asm_arccos( float r ) {
  // return half_pi + arctan( r / -sqr( 1.f - r * r ) );
  
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
  
  float asm_one = 1.f;
  float asm_half_pi = half_pi;
  __asm {
    fld r // r0 = r
    fld r // r1 = r0, r0 = r
    fmul r // r0 = r0 * r
    fsubr asm_one // r0 = r0 - 1.f
    fsqrt // r0 = sqrtf( r0 )
    fchs // r0 = - r0
    fdiv // r0 = r1 / r0
    fld1 // {{ r0 = atan( r0 )
    fpatan // }}
    fadd asm_half_pi // r0 = r0 + pi / 2
  } // returns r0
  
#else
  
  return float( acos( r ) );
  
#endif
}

float asm_arcsin( float r ) {
  // return arctan( r / sqr( 1.f - r * r ) );
  
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
  
  const float asm_one = 1.f;
  __asm {
    fld r // r0 = r
    fld r // r1 = r0, r0 = r
    fmul r // r0 = r0 * r
    fsubr asm_one // r0 = r0 - 1.f
    fsqrt // r0 = sqrtf( r0 )
    fdiv // r0 = r1 / r0
    fld1 // {{ r0 = atan( r0 )
    fpatan // }}
  } // returns r0
  
#else
  
  return float( asin( r ) );
  
#endif
  
}

float asm_arctan( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  __asm {
    fld r // r0 = r
    fld1 // {{ r0 = atan( r0 )
    fpatan // }}
  } // returns r0
  
#else
  
  return float( atan( r ) );
  
#endif
  
}

float asm_sin( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  __asm {
    fld r // r0 = r
    fsin // r0 = sinf( r0 )
  } // returns r0
  
#else
  
  return sin( r );
  
#endif
  
}

float asm_cos( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  __asm {
    fld r // r0 = r
    fcos // r0 = cosf( r0 )
  } // returns r0
  
#else
  
  return cos( r );
  
#endif
}

float asm_tan( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  // return sin( r ) / cos( r );
  __asm {
    fld r // r0 = r
    fsin // r0 = sinf( r0 )
    fld r // r1 = r0, r0 = r
    fcos // r0 = cosf( r0 )
    fdiv // r0 = r1 / r0
  } // returns r0
  
#else
  
  return tan( r );
  
#endif
}

// returns a for a * a = r
float asm_sqrt( float r )
{
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  __asm {
    fld r // r0 = r
    fsqrt // r0 = sqrtf( r0 )
  } // returns r0
  
#else
  
  return sqrt( r );
  
#endif
}

// returns 1 / a for a * a = r
// -- Use this for Vector normalisation!!!
float asm_rsq( float r )
{
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  __asm {
    fld1 // r0 = 1.f
    fld r // r1 = r0, r0 = r
    fsqrt // r0 = sqrtf( r0 )
    fdiv // r0 = r1 / r0
  } // returns r0
  
#else
  
  return 1. / sqrt( r );
  
#endif
}

// returns 1 / a for a * a = r
// Another version
float apx_rsq( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  const float asm_dot5 = 0.5f;
  const float asm_1dot5 = 1.5f;
  
  __asm {
    fld r // r0 = r
    fmul asm_dot5 // r0 = r0 * .5f
    mov eax, r // eax = r
    shr eax, 0x1 // eax = eax >> 1
    neg eax // eax = -eax
    add eax, 0x5F400000 // eax = eax & MAGICAL NUMBER
    mov r, eax // r = eax
    fmul r // r0 = r0 * r
    fmul r // r0 = r0 * r
    fsubr asm_1dot5 // r0 = 1.5f - r0
    fmul r // r0 = r0 * r
  } // returns r0
  
#else
  
  return 1. / sqrt( r );
  
#endif
}

/* very MS-specific, commented out for now
   Finally the best InvSqrt implementation?
   Use for vector normalisation instead of 1/length() * x,y,z
*/
#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

__declspec(naked) float __fastcall InvSqrt(float fValue)
{
  __asm
  {
    mov        eax, 0be6eb508h
    mov        dword ptr[esp-12],03fc00000h
    sub        eax, dword ptr[esp + 4]
    sub        dword ptr[esp+4], 800000h
    shr        eax, 1
    mov        dword ptr[esp -  8], eax
    
    fld        dword ptr[esp -  8]
    fmul    st, st
    fld        dword ptr[esp -  8]
    fxch    st(1)
    fmul    dword ptr[esp +  4]
    fld        dword ptr[esp - 12]
    fld        st(0)
    fsub    st,st(2)
    
    fld        st(1)
    fxch    st(1)
    fmul    st(3),st
    fmul    st(3),st
    fmulp    st(4),st
    fsub    st,st(2)
    
    fmul    st(2),st
    fmul    st(3),st
    fmulp    st(2),st
    fxch    st(1)
    fsubp    st(1),st
    
    fmulp    st(1), st
    ret 4
  }
}

#endif

// returns a random number
OGREMATHLIB_FORCEINLINE float asm_rand()
{

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
#if 0
#if OGREMATHLIB_COMP_VER >= 1300

  static unsigned __int64 q = time( NULL );
  
  _asm {
    movq mm0, q
    
    // do the magic MMX thing
    pshufw mm1, mm0, 0x1E
    paddd mm0, mm1
    
    // move to integer memory location and free MMX
    movq q, mm0
    emms
  }
  
  return float( q );
#endif
#else
  // VC6 does not support pshufw
  return float( rand() );
#endif
#else
  // GCC etc
  
  return float( rand() );
  
#endif
}

// returns the maximum random number
OGREMATHLIB_FORCEINLINE float asm_rand_max()
{

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
#if 0
#if OGREMATHLIB_COMP_VER >= 1300

  return (std::numeric_limits< unsigned __int64 >::max)();
  return 9223372036854775807.0f;
#endif
#else
  // VC6 does not support unsigned __int64
  return float( RAND_MAX );
#endif
  
#else
  // GCC etc
  return float( RAND_MAX );
  
#endif
}

// returns log2( r ) / log2( e )
float asm_ln( float r ) {

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32

  const float asm_1_div_log2_e = .693147180559f;
  const float asm_neg1_div_3 = -.33333333333333333333333333333f;
  const float asm_neg2_div_3 = -.66666666666666666666666666667f;
  const float asm_2 = 2.f;
  
  int log_2 = 0;
  
  __asm {
    // log_2 = ( ( r >> 0x17 ) & 0xFF ) - 0x80;
    mov eax, r
    sar eax, 0x17
    and eax, 0xFF
    sub eax, 0x80
    mov log_2, eax
    
    // r = ( r & 0x807fffff ) + 0x3f800000;
    mov ebx, r
    and ebx, 0x807FFFFF
    add ebx, 0x3F800000
    mov r, ebx
    
    // r = ( asm_neg1_div_3 * r + asm_2 ) * r + asm_neg2_div_3;   // (1)
    fld r
    fmul asm_neg1_div_3
    fadd asm_2
    fmul r
    fadd asm_neg2_div_3
    fild log_2
    fadd
    fmul asm_1_div_log2_e
  }
  
#else
  
  return log( r );
  
#endif
}

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC &&  OGREMATHLIB_ARCH_TYPE == OGREMATHLIB_ARCHITECTURE_32
#	pragma warning( pop )
#endif

#if  OGREMATHLIB_COMPILER == OGREMATHLIB_COMPILER_MSVC
#  pragma warning (pop)
#endif

const Real Math::POS_INFINITY = std::numeric_limits<Real>::infinity();
const Real Math::NEG_INFINITY = -std::numeric_limits<Real>::infinity();
const Real Math::PI = Real( 4.0 * atan( 1.0 ) );
const Real Math::TWO_PI = Real( 2.0 * PI );
const Real Math::HALF_PI = Real( 0.5 * PI );
const Real Math::fDeg2Rad = PI / Real(180.0);
const Real Math::fRad2Deg = Real(180.0) / PI;
const Real Math::LOG2 = log(Real(2.0));

int Math::mTrigTableSize;

Real  Math::mTrigTableFactor;
Real *Math::mSinTable = NULL;
Real *Math::mTanTable = NULL;

// these functions must be defined down here, because they rely on the
// angle unit conversion functions in class Math:

inline Real Radian::valueDegrees() const
{
  return Math::RadiansToDegrees ( mRad );
}

inline Real Degree::valueRadians() const
{
  return Math::DegreesToRadians ( mDeg );
}
//-----------------------------------------------------------------------
Math::Math( unsigned int trigTableSize )
{

  mTrigTableSize = trigTableSize;
  mTrigTableFactor = mTrigTableSize / Math::TWO_PI;
  
  mSinTable = OGRELIBMATH_ALLOC_T(Real, mTrigTableSize, MEMCATEGORY_GENERAL);
  mTanTable = OGRELIBMATH_ALLOC_T(Real, mTrigTableSize, MEMCATEGORY_GENERAL);
  
  buildTrigTables();
}

//-----------------------------------------------------------------------
Math::~Math()
{
  OGREMATHLIB_FREE(mSinTable, MEMCATEGORY_GENERAL);
  OGREMATHLIB_FREE(mTanTable, MEMCATEGORY_GENERAL);
}

//-----------------------------------------------------------------------
void Math::buildTrigTables(void)
{
  // Build trig lookup tables
  // Could get away with building only PI sized Sin table but simpler this
  // way. Who cares, it'll ony use an extra 8k of memory anyway and I like
  // simplicity.
  Real angle;
  for (int i = 0; i < mTrigTableSize; ++i)
  {
    angle = Math::TWO_PI * i / mTrigTableSize;
    mSinTable[i] = sin(angle);
    mTanTable[i] = tan(angle);
  }
}
//-----------------------------------------------------------------------
Real Math::SinTable (Real fValue)
{
  // Convert range to index values, wrap if required
  int idx;
  if (fValue >= 0)
  {
    idx = int(fValue * mTrigTableFactor) % mTrigTableSize;
  }
  else
  {
    idx = mTrigTableSize - (int(-fValue * mTrigTableFactor) % mTrigTableSize) - 1;
  }
  
  return mSinTable[idx];
}
//-----------------------------------------------------------------------
Real Math::TanTable (Real fValue)
{
  // Convert range to index values, wrap if required
  int idx = int(fValue *= mTrigTableFactor) % mTrigTableSize;
  return mTanTable[idx];
}
//-----------------------------------------------------------------------
int Math::ISign (int iValue)
{
  return ( iValue > 0 ? +1 : ( iValue < 0 ? -1 : 0 ) );
}
//-----------------------------------------------------------------------
Radian Math::ACos (Real fValue)
{
  if ( -1.0 < fValue )
  {
    if ( fValue < 1.0 )
      return Radian(acos(fValue));
    else
      return Radian(0.0);
  }
  else
  {
    return Radian(PI);
  }
}
//-----------------------------------------------------------------------
Radian Math::ASin (Real fValue)
{
  if ( -1.0 < fValue )
  {
    if ( fValue < 1.0 )
      return Radian(asin(fValue));
    else
      return Radian(HALF_PI);
  }
  else
  {
    return Radian(-HALF_PI);
  }
}
//-----------------------------------------------------------------------
Real Math::Sign (Real fValue)
{
  if ( fValue > 0.0 )
    return 1.0;
    
  if ( fValue < 0.0 )
    return -1.0;
    
  return 0.0;
}
//-----------------------------------------------------------------------
Real Math::InvSqrt(Real fValue)
{
  return Real(asm_rsq(fValue));
}
//-----------------------------------------------------------------------
Real Math::UnitRandom ()
{
  return asm_rand() / asm_rand_max();
}

//-----------------------------------------------------------------------
Real Math::RangeRandom (Real fLow, Real fHigh)
{
  return (fHigh-fLow)*UnitRandom() + fLow;
}

//-----------------------------------------------------------------------
Real Math::SymmetricRandom ()
{
  return 2.0f * UnitRandom() - 1.0f;
}

//-----------------------------------------------------------------------
bool Math::RealEqual( Real a, Real b, Real tolerance )
{
  if (fabs(b-a) <= tolerance)
    return true;
  else
    return false;
}



inline Vector3 Vector3::randomDeviant(
  const Radian& angle,
  const Vector3& up  ) const
{
  Vector3 newUp;
  
  if (up == Vector3::ZERO)
  {
    // Generate an up vector
    newUp = this->perpendicular();
  }
  else
  {
    newUp = up;
  }
  
  // Rotate up vector by random amount around this
  Quaternion q;
  q.FromAngleAxis( Radian(Math::UnitRandom() * Math::TWO_PI), *this );
  newUp = q * newUp;
  
  // Finally rotate this by given angle around randomised up
  q.FromAngleAxis( angle, newUp );
  return q * (*this);
}

inline Radian Vector3::angleBetween(const Vector3& dest)
{
  Real lenProduct = length() * dest.length();
  
  // Divide by zero check
  if (lenProduct < 1e-6f)
    lenProduct = 1e-6f;
    
  Real f = dotProduct(dest) / lenProduct;
  
  f = Math::Clamp(f, (Real)-1.0, (Real)1.0);
  return Math::ACos(f);
  
}

const Vector3 Vector3::ZERO( 0, 0, 0 );

const Vector3 Vector3::UNIT_X( 1, 0, 0 );
const Vector3 Vector3::UNIT_Y( 0, 1, 0 );
const Vector3 Vector3::UNIT_Z( 0, 0, 1 );
const Vector3 Vector3::NEGATIVE_UNIT_X( -1,  0,  0 );
const Vector3 Vector3::NEGATIVE_UNIT_Y(  0, -1,  0 );
const Vector3 Vector3::NEGATIVE_UNIT_Z(  0,  0, -1 );
const Vector3 Vector3::UNIT_SCALE(1, 1, 1);

Quaternion Vector3::getRotationTo(const Vector3& dest,
                                  const Vector3& fallbackAxis) const
{
  // Based on Stan Melax's article in Game Programming Gems
  Quaternion q;
  // Copy, since cannot modify local
  Vector3 v0 = *this;
  Vector3 v1 = dest;
  v0.normalise();
  v1.normalise();
  
  Real d = v0.dotProduct(v1);
  // If dot == 1, vectors are the same
  if (d >= 1.0f)
  {
    return Quaternion::IDENTITY;
  }
  if (d < (1e-6f - 1.0f))
  {
    if (fallbackAxis != Vector3::ZERO)
    {
      // rotate 180 degrees about the fallback axis
      q.FromAngleAxis(Radian(Math::PI), fallbackAxis);
    }
    else
    {
      // Generate an axis
      Vector3 axis = Vector3::UNIT_X.crossProduct(*this);
      if (axis.isZeroLength()) // pick another if colinear
        axis = Vector3::UNIT_Y.crossProduct(*this);
      axis.normalise();
      q.FromAngleAxis(Radian(Math::PI), axis);
    }
  }
  else
  {
    Real s = Math::Sqrt( (1+d)*2 );
    Real invs = 1 / s;
    
    Vector3 c = v0.crossProduct(v1);
    
    q.x = c.x * invs;
    q.y = c.y * invs;
    q.z = c.z * invs;
    q.w = s * 0.5f;
    q.normalise();
  }
  return q;
}

const Real Quaternion::ms_fEpsilon = 1e-03f;
const Quaternion Quaternion::ZERO(0,0,0,0);
const Quaternion Quaternion::IDENTITY(1,0,0,0);

//-----------------------------------------------------------------------
void Quaternion::FromRotationMatrix (const Matrix3& kRot)
{
  // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
  // article "Quaternion Calculus and Fast Animation".
  
  Real fTrace = kRot[0][0]+kRot[1][1]+kRot[2][2];
  Real fRoot;
  
  if ( fTrace > 0.0 )
  {
    // |w| > 1/2, may as well choose w > 1/2
    fRoot = Math::Sqrt(fTrace + 1.0f);  // 2w
    w = 0.5f*fRoot;
    fRoot = 0.5f/fRoot;  // 1/(4w)
    x = (kRot[2][1]-kRot[1][2])*fRoot;
    y = (kRot[0][2]-kRot[2][0])*fRoot;
    z = (kRot[1][0]-kRot[0][1])*fRoot;
  }
  else
  {
    // |w| <= 1/2
    static size_t s_iNext[3] = { 1, 2, 0 };
    size_t i = 0;
    if ( kRot[1][1] > kRot[0][0] )
      i = 1;
    if ( kRot[2][2] > kRot[i][i] )
      i = 2;
    size_t j = s_iNext[i];
    size_t k = s_iNext[j];
    
    fRoot = Math::Sqrt(kRot[i][i]-kRot[j][j]-kRot[k][k] + 1.0f);
    Real* apkQuat[3] = { &x, &y, &z };
    *apkQuat[i] = 0.5f*fRoot;
    fRoot = 0.5f/fRoot;
    w = (kRot[k][j]-kRot[j][k])*fRoot;
    *apkQuat[j] = (kRot[j][i]+kRot[i][j])*fRoot;
    *apkQuat[k] = (kRot[k][i]+kRot[i][k])*fRoot;
  }
}
//-----------------------------------------------------------------------
void Quaternion::ToRotationMatrix (Matrix3& kRot) const
{
  Real fTx  = x+x;
  Real fTy  = y+y;
  Real fTz  = z+z;
  Real fTwx = fTx*w;
  Real fTwy = fTy*w;
  Real fTwz = fTz*w;
  Real fTxx = fTx*x;
  Real fTxy = fTy*x;
  Real fTxz = fTz*x;
  Real fTyy = fTy*y;
  Real fTyz = fTz*y;
  Real fTzz = fTz*z;
  
  kRot[0][0] = 1.0f-(fTyy+fTzz);
  kRot[0][1] = fTxy-fTwz;
  kRot[0][2] = fTxz+fTwy;
  kRot[1][0] = fTxy+fTwz;
  kRot[1][1] = 1.0f-(fTxx+fTzz);
  kRot[1][2] = fTyz-fTwx;
  kRot[2][0] = fTxz-fTwy;
  kRot[2][1] = fTyz+fTwx;
  kRot[2][2] = 1.0f-(fTxx+fTyy);
}
//-----------------------------------------------------------------------
void Quaternion::FromAngleAxis (const Radian& rfAngle,
                                const Vector3& rkAxis)
{
  // assert:  axis[] is unit length
  //
  // The quaternion representing the rotation is
  //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
  
  Radian fHalfAngle ( 0.5f * rfAngle );
  Real fSin = Math::Sin(fHalfAngle);
  w = Math::Cos(fHalfAngle);
  x = fSin*rkAxis.x;
  y = fSin*rkAxis.y;
  z = fSin*rkAxis.z;
}
//-----------------------------------------------------------------------
void Quaternion::ToAngleAxis (Radian& rfAngle, Vector3& rkAxis) const
{
  // The quaternion representing the rotation is
  //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
  
  Real fSqrLength = x*x+y*y+z*z;
  if ( fSqrLength > 0.0 )
  {
    rfAngle = 2.0*Math::ACos(w);
    Real fInvLength = Math::InvSqrt(fSqrLength);
    rkAxis.x = x*fInvLength;
    rkAxis.y = y*fInvLength;
    rkAxis.z = z*fInvLength;
  }
  else
  {
    // angle is 0 (mod 2*pi), so any axis will do
    rfAngle = Radian(0.0);
    rkAxis.x = 1.0;
    rkAxis.y = 0.0;
    rkAxis.z = 0.0;
  }
}
//-----------------------------------------------------------------------
void Quaternion::FromAxes (const Vector3* akAxis)
{
  Matrix3 kRot;
  
  for (size_t iCol = 0; iCol < 3; iCol++)
  {
    kRot[0][iCol] = akAxis[iCol].x;
    kRot[1][iCol] = akAxis[iCol].y;
    kRot[2][iCol] = akAxis[iCol].z;
  }
  
  FromRotationMatrix(kRot);
}
//-----------------------------------------------------------------------
void Quaternion::FromAxes (const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis)
{
  Matrix3 kRot;
  
  kRot[0][0] = xaxis.x;
  kRot[1][0] = xaxis.y;
  kRot[2][0] = xaxis.z;
  
  kRot[0][1] = yaxis.x;
  kRot[1][1] = yaxis.y;
  kRot[2][1] = yaxis.z;
  
  kRot[0][2] = zaxis.x;
  kRot[1][2] = zaxis.y;
  kRot[2][2] = zaxis.z;
  
  FromRotationMatrix(kRot);
  
}
//-----------------------------------------------------------------------
void Quaternion::ToAxes (Vector3* akAxis) const
{
  Matrix3 kRot;
  
  ToRotationMatrix(kRot);
  
  for (size_t iCol = 0; iCol < 3; iCol++)
  {
    akAxis[iCol].x = kRot[0][iCol];
    akAxis[iCol].y = kRot[1][iCol];
    akAxis[iCol].z = kRot[2][iCol];
  }
}
//-----------------------------------------------------------------------
Vector3 Quaternion::xAxis(void) const
{
  //Real fTx  = 2.0*x;
  Real fTy  = 2.0f*y;
  Real fTz  = 2.0f*z;
  Real fTwy = fTy*w;
  Real fTwz = fTz*w;
  Real fTxy = fTy*x;
  Real fTxz = fTz*x;
  Real fTyy = fTy*y;
  Real fTzz = fTz*z;
  
  return Vector3(1.0f-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
}
//-----------------------------------------------------------------------
Vector3 Quaternion::yAxis(void) const
{
  Real fTx  = 2.0f*x;
  Real fTy  = 2.0f*y;
  Real fTz  = 2.0f*z;
  Real fTwx = fTx*w;
  Real fTwz = fTz*w;
  Real fTxx = fTx*x;
  Real fTxy = fTy*x;
  Real fTyz = fTz*y;
  Real fTzz = fTz*z;
  
  return Vector3(fTxy-fTwz, 1.0f-(fTxx+fTzz), fTyz+fTwx);
}
//-----------------------------------------------------------------------
Vector3 Quaternion::zAxis(void) const
{
  Real fTx  = 2.0f*x;
  Real fTy  = 2.0f*y;
  Real fTz  = 2.0f*z;
  Real fTwx = fTx*w;
  Real fTwy = fTy*w;
  Real fTxx = fTx*x;
  Real fTxz = fTz*x;
  Real fTyy = fTy*y;
  Real fTyz = fTz*y;
  
  return Vector3(fTxz+fTwy, fTyz-fTwx, 1.0f-(fTxx+fTyy));
}
//-----------------------------------------------------------------------
void Quaternion::ToAxes (Vector3& xaxis, Vector3& yaxis, Vector3& zaxis) const
{
  Matrix3 kRot;
  
  ToRotationMatrix(kRot);
  
  xaxis.x = kRot[0][0];
  xaxis.y = kRot[1][0];
  xaxis.z = kRot[2][0];
  
  yaxis.x = kRot[0][1];
  yaxis.y = kRot[1][1];
  yaxis.z = kRot[2][1];
  
  zaxis.x = kRot[0][2];
  zaxis.y = kRot[1][2];
  zaxis.z = kRot[2][2];
}

//-----------------------------------------------------------------------
Quaternion Quaternion::operator+ (const Quaternion& rkQ) const
{
  return Quaternion(w+rkQ.w,x+rkQ.x,y+rkQ.y,z+rkQ.z);
}
//-----------------------------------------------------------------------
Quaternion Quaternion::operator- (const Quaternion& rkQ) const
{
  return Quaternion(w-rkQ.w,x-rkQ.x,y-rkQ.y,z-rkQ.z);
}
//-----------------------------------------------------------------------
Quaternion Quaternion::operator* (const Quaternion& rkQ) const
{
  // NOTE:  Multiplication is not generally commutative, so in most
  // cases p*q != q*p.
  
  return Quaternion
         (
           w * rkQ.w - x * rkQ.x - y * rkQ.y - z * rkQ.z,
           w * rkQ.x + x * rkQ.w + y * rkQ.z - z * rkQ.y,
           w * rkQ.y + y * rkQ.w + z * rkQ.x - x * rkQ.z,
           w * rkQ.z + z * rkQ.w + x * rkQ.y - y * rkQ.x
         );
}
//-----------------------------------------------------------------------
Quaternion Quaternion::operator* (Real fScalar) const
{
  return Quaternion(fScalar*w,fScalar*x,fScalar*y,fScalar*z);
}
//-----------------------------------------------------------------------
Quaternion operator* (Real fScalar, const Quaternion& rkQ)
{
  return Quaternion(fScalar*rkQ.w,fScalar*rkQ.x,fScalar*rkQ.y,
                    fScalar*rkQ.z);
}
//-----------------------------------------------------------------------
Quaternion Quaternion::operator- () const
{
  return Quaternion(-w,-x,-y,-z);
}
//-----------------------------------------------------------------------
Real Quaternion::Dot (const Quaternion& rkQ) const
{
  return w*rkQ.w+x*rkQ.x+y*rkQ.y+z*rkQ.z;
}
//-----------------------------------------------------------------------
Real Quaternion::Norm () const
{
  return w*w+x*x+y*y+z*z;
}
//-----------------------------------------------------------------------
Quaternion Quaternion::Inverse () const
{
  Real fNorm = w*w+x*x+y*y+z*z;
  if ( fNorm > 0.0 )
  {
    Real fInvNorm = 1.0f/fNorm;
    return Quaternion(w*fInvNorm,-x*fInvNorm,-y*fInvNorm,-z*fInvNorm);
  }
  else
  {
    // return an invalid result to flag the error
    return ZERO;
  }
}
//-----------------------------------------------------------------------
Quaternion Quaternion::UnitInverse () const
{
  // assert:  'this' is unit length
  return Quaternion(w,-x,-y,-z);
}
//-----------------------------------------------------------------------
Quaternion Quaternion::Exp () const
{
  // If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
  // exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
  // use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.
  
  Radian fAngle ( Math::Sqrt(x*x+y*y+z*z) );
  Real fSin = Math::Sin(fAngle);
  
  Quaternion kResult;
  kResult.w = Math::Cos(fAngle);
  
  if ( Math::Abs(fSin) >= ms_fEpsilon )
  {
    Real fCoeff = fSin/(fAngle.valueRadians());
    kResult.x = fCoeff*x;
    kResult.y = fCoeff*y;
    kResult.z = fCoeff*z;
  }
  else
  {
    kResult.x = x;
    kResult.y = y;
    kResult.z = z;
  }
  
  return kResult;
}
//-----------------------------------------------------------------------
Quaternion Quaternion::Log () const
{
  // If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
  // log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
  // sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.
  
  Quaternion kResult;
  kResult.w = 0.0;
  
  if ( Math::Abs(w) < 1.0 )
  {
    Radian fAngle ( Math::ACos(w) );
    Real fSin = Math::Sin(fAngle);
    if ( Math::Abs(fSin) >= ms_fEpsilon )
    {
      Real fCoeff = fAngle.valueRadians()/fSin;
      kResult.x = fCoeff*x;
      kResult.y = fCoeff*y;
      kResult.z = fCoeff*z;
      return kResult;
    }
  }
  
  kResult.x = x;
  kResult.y = y;
  kResult.z = z;
  
  return kResult;
}
//-----------------------------------------------------------------------
Vector3 Quaternion::operator* (const Vector3& v) const
{
  // nVidia SDK implementation
  Vector3 uv, uuv;
  Vector3 qvec(x, y, z);
  uv = qvec.crossProduct(v);
  uuv = qvec.crossProduct(uv);
  uv *= (2.0f * w);
  uuv *= 2.0f;
  
  return v + uv + uuv;
  
}
//-----------------------------------------------------------------------
bool Quaternion::equals(const Quaternion& rhs, const Radian& tolerance) const
{
  Real fCos = Dot(rhs);
  Radian angle = Math::ACos(fCos);
  
  return (Math::Abs(angle.valueRadians()) <= tolerance.valueRadians())
         || Math::RealEqual(angle.valueRadians(), Math::PI, tolerance.valueRadians());
         
         
}
//-----------------------------------------------------------------------
Quaternion Quaternion::Slerp (Real fT, const Quaternion& rkP,
                              const Quaternion& rkQ, bool shortestPath)
{
  Real fCos = rkP.Dot(rkQ);
  Quaternion rkT;
  
  // Do we need to invert rotation?
  if (fCos < 0.0f && shortestPath)
  {
    fCos = -fCos;
    rkT = -rkQ;
  }
  else
  {
    rkT = rkQ;
  }
  
  if (Math::Abs(fCos) < 1 - ms_fEpsilon)
  {
    // Standard case (slerp)
    Real fSin = Math::Sqrt(1 - Math::Sqr(fCos));
    Radian fAngle = Math::ATan2(fSin, fCos);
    Real fInvSin = 1.0f / fSin;
    Real fCoeff0 = Math::Sin((1.0f - fT) * fAngle) * fInvSin;
    Real fCoeff1 = Math::Sin(fT * fAngle) * fInvSin;
    return fCoeff0 * rkP + fCoeff1 * rkT;
  }
  else
  {
    // There are two situations:
    // 1. "rkP" and "rkQ" are very close (fCos ~= +1), so we can do a linear
    //    interpolation safely.
    // 2. "rkP" and "rkQ" are almost inverse of each other (fCos ~= -1), there
    //    are an infinite number of possibilities interpolation. but we haven't
    //    have method to fix this case, so just use linear interpolation here.
    Quaternion t = (1.0f - fT) * rkP + fT * rkT;
    // taking the complement requires renormalisation
    t.normalise();
    return t;
  }
}
//-----------------------------------------------------------------------
Quaternion Quaternion::SlerpExtraSpins (Real fT,
                                        const Quaternion& rkP, const Quaternion& rkQ, int iExtraSpins)
{
  Real fCos = rkP.Dot(rkQ);
  Radian fAngle ( Math::ACos(fCos) );
  
  if ( Math::Abs(fAngle.valueRadians()) < ms_fEpsilon )
    return rkP;
    
  Real fSin = Math::Sin(fAngle);
  Radian fPhase ( Math::PI*iExtraSpins*fT );
  Real fInvSin = 1.0f/fSin;
  Real fCoeff0 = Math::Sin((1.0f-fT)*fAngle - fPhase)*fInvSin;
  Real fCoeff1 = Math::Sin(fT*fAngle + fPhase)*fInvSin;
  return fCoeff0*rkP + fCoeff1*rkQ;
}
//-----------------------------------------------------------------------
void Quaternion::Intermediate (const Quaternion& rkQ0,
                               const Quaternion& rkQ1, const Quaternion& rkQ2,
                               Quaternion& rkA, Quaternion& rkB)
{
  // assert:  q0, q1, q2 are unit quaternions
  
  Quaternion kQ0inv = rkQ0.UnitInverse();
  Quaternion kQ1inv = rkQ1.UnitInverse();
  Quaternion rkP0 = kQ0inv*rkQ1;
  Quaternion rkP1 = kQ1inv*rkQ2;
  Quaternion kArg = 0.25*(rkP0.Log()-rkP1.Log());
  Quaternion kMinusArg = -kArg;
  
  rkA = rkQ1*kArg.Exp();
  rkB = rkQ1*kMinusArg.Exp();
}
//-----------------------------------------------------------------------
Quaternion Quaternion::Squad (Real fT,
                              const Quaternion& rkP, const Quaternion& rkA,
                              const Quaternion& rkB, const Quaternion& rkQ, bool shortestPath)
{
  Real fSlerpT = 2.0f*fT*(1.0f-fT);
  Quaternion kSlerpP = Slerp(fT, rkP, rkQ, shortestPath);
  Quaternion kSlerpQ = Slerp(fT, rkA, rkB);
  return Slerp(fSlerpT, kSlerpP ,kSlerpQ);
}
//-----------------------------------------------------------------------
Real Quaternion::normalise(void)
{
  Real len = Norm();
  Real factor = 1.0f / Math::Sqrt(len);
  *this = *this * factor;
  return len;
}
//-----------------------------------------------------------------------
Radian Quaternion::getRoll(bool reprojectAxis) const
{
  if (reprojectAxis)
  {
    // roll = atan2(localx.y, localx.x)
    // pick parts of xAxis() implementation that we need
//			Real fTx  = 2.0*x;
    Real fTy  = 2.0f*y;
    Real fTz  = 2.0f*z;
    Real fTwz = fTz*w;
    Real fTxy = fTy*x;
    Real fTyy = fTy*y;
    Real fTzz = fTz*z;
    
    // Vector3(1.0-(fTyy+fTzz), fTxy+fTwz, fTxz-fTwy);
    
    return Radian(Math::ATan2(fTxy+fTwz, 1.0f-(fTyy+fTzz)));
    
  }
  else
  {
    return Radian(Math::ATan2(2*(x*y + w*z), w*w + x*x - y*y - z*z));
  }
}
//-----------------------------------------------------------------------
Radian Quaternion::getPitch(bool reprojectAxis) const
{
  if (reprojectAxis)
  {
    // pitch = atan2(localy.z, localy.y)
    // pick parts of yAxis() implementation that we need
    Real fTx  = 2.0f*x;
//			Real fTy  = 2.0f*y;
    Real fTz  = 2.0f*z;
    Real fTwx = fTx*w;
    Real fTxx = fTx*x;
    Real fTyz = fTz*y;
    Real fTzz = fTz*z;
    
    // Vector3(fTxy-fTwz, 1.0-(fTxx+fTzz), fTyz+fTwx);
    return Radian(Math::ATan2(fTyz+fTwx, 1.0f-(fTxx+fTzz)));
  }
  else
  {
    // internal version
    return Radian(Math::ATan2(2*(y*z + w*x), w*w - x*x - y*y + z*z));
  }
}
//-----------------------------------------------------------------------
Radian Quaternion::getYaw(bool reprojectAxis) const
{
  if (reprojectAxis)
  {
    // yaw = atan2(localz.x, localz.z)
    // pick parts of zAxis() implementation that we need
    Real fTx  = 2.0f*x;
    Real fTy  = 2.0f*y;
    Real fTz  = 2.0f*z;
    Real fTwy = fTy*w;
    Real fTxx = fTx*x;
    Real fTxz = fTz*x;
    Real fTyy = fTy*y;
    
    // Vector3(fTxz+fTwy, fTyz-fTwx, 1.0-(fTxx+fTyy));
    
    return Radian(Math::ATan2(fTxz+fTwy, 1.0f-(fTxx+fTyy)));
    
  }
  else
  {
    // internal version
    return Radian(Math::ASin(-2*(x*z - w*y)));
  }
}
//-----------------------------------------------------------------------
Quaternion Quaternion::nlerp(Real fT, const Quaternion& rkP,
                             const Quaternion& rkQ, bool shortestPath)
{
  Quaternion result;
  Real fCos = rkP.Dot(rkQ);
  if (fCos < 0.0f && shortestPath)
  {
    result = rkP + fT * ((-rkQ) - rkP);
  }
  else
  {
    result = rkP + fT * (rkQ - rkP);
  }
  result.normalise();
  return result;
}

const Real Matrix3::EPSILON = 1e-06f;
const Matrix3 Matrix3::ZERO(0,0,0,0,0,0,0,0,0);
const Matrix3 Matrix3::IDENTITY(1,0,0,0,1,0,0,0,1);
const Real Matrix3::ms_fSvdEpsilon = 1e-04f;
const unsigned int Matrix3::ms_iSvdMaxIterations = 32;

//-----------------------------------------------------------------------
Vector3 Matrix3::GetColumn (size_t iCol) const
{
  assert( 0 <= iCol && iCol < 3 );
  return Vector3(m[0][iCol],m[1][iCol],
                 m[2][iCol]);
}
//-----------------------------------------------------------------------
void Matrix3::SetColumn(size_t iCol, const Vector3& vec)
{
  assert( 0 <= iCol && iCol < 3 );
  m[0][iCol] = vec.x;
  m[1][iCol] = vec.y;
  m[2][iCol] = vec.z;
  
}
//-----------------------------------------------------------------------
void Matrix3::FromAxes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis)
{
  SetColumn(0,xAxis);
  SetColumn(1,yAxis);
  SetColumn(2,zAxis);
  
}

//-----------------------------------------------------------------------
bool Matrix3::operator== (const Matrix3& rkMatrix) const
{
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
    {
      if ( m[iRow][iCol] != rkMatrix.m[iRow][iCol] )
        return false;
    }
  }
  
  return true;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::operator+ (const Matrix3& rkMatrix) const
{
  Matrix3 kSum;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
    {
      kSum.m[iRow][iCol] = m[iRow][iCol] +
                           rkMatrix.m[iRow][iCol];
    }
  }
  return kSum;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::operator- (const Matrix3& rkMatrix) const
{
  Matrix3 kDiff;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
    {
      kDiff.m[iRow][iCol] = m[iRow][iCol] -
                            rkMatrix.m[iRow][iCol];
    }
  }
  return kDiff;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::operator* (const Matrix3& rkMatrix) const
{
  Matrix3 kProd;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
    {
      kProd.m[iRow][iCol] =
        m[iRow][0]*rkMatrix.m[0][iCol] +
        m[iRow][1]*rkMatrix.m[1][iCol] +
        m[iRow][2]*rkMatrix.m[2][iCol];
    }
  }
  return kProd;
}
//-----------------------------------------------------------------------
Vector3 Matrix3::operator* (const Vector3& rkPoint) const
{
  Vector3 kProd;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    kProd[iRow] =
      m[iRow][0]*rkPoint[0] +
      m[iRow][1]*rkPoint[1] +
      m[iRow][2]*rkPoint[2];
  }
  return kProd;
}
//-----------------------------------------------------------------------
Vector3 operator* (const Vector3& rkPoint, const Matrix3& rkMatrix)
{
  Vector3 kProd;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    kProd[iRow] =
      rkPoint[0]*rkMatrix.m[0][iRow] +
      rkPoint[1]*rkMatrix.m[1][iRow] +
      rkPoint[2]*rkMatrix.m[2][iRow];
  }
  return kProd;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::operator- () const
{
  Matrix3 kNeg;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      kNeg[iRow][iCol] = -m[iRow][iCol];
  }
  return kNeg;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::operator* (Real fScalar) const
{
  Matrix3 kProd;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      kProd[iRow][iCol] = fScalar*m[iRow][iCol];
  }
  return kProd;
}
//-----------------------------------------------------------------------
Matrix3 operator* (Real fScalar, const Matrix3& rkMatrix)
{
  Matrix3 kProd;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      kProd[iRow][iCol] = fScalar*rkMatrix.m[iRow][iCol];
  }
  return kProd;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::Transpose () const
{
  Matrix3 kTranspose;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      kTranspose[iRow][iCol] = m[iCol][iRow];
  }
  return kTranspose;
}
//-----------------------------------------------------------------------
bool Matrix3::Inverse (Matrix3& rkInverse, Real fTolerance) const
{
  // Invert a 3x3 using cofactors.  This is about 8 times faster than
  // the Numerical Recipes code which uses Gaussian elimination.
  
  rkInverse[0][0] = m[1][1]*m[2][2] -
                    m[1][2]*m[2][1];
  rkInverse[0][1] = m[0][2]*m[2][1] -
                    m[0][1]*m[2][2];
  rkInverse[0][2] = m[0][1]*m[1][2] -
                    m[0][2]*m[1][1];
  rkInverse[1][0] = m[1][2]*m[2][0] -
                    m[1][0]*m[2][2];
  rkInverse[1][1] = m[0][0]*m[2][2] -
                    m[0][2]*m[2][0];
  rkInverse[1][2] = m[0][2]*m[1][0] -
                    m[0][0]*m[1][2];
  rkInverse[2][0] = m[1][0]*m[2][1] -
                    m[1][1]*m[2][0];
  rkInverse[2][1] = m[0][1]*m[2][0] -
                    m[0][0]*m[2][1];
  rkInverse[2][2] = m[0][0]*m[1][1] -
                    m[0][1]*m[1][0];
                    
  Real fDet =
    m[0][0]*rkInverse[0][0] +
    m[0][1]*rkInverse[1][0]+
    m[0][2]*rkInverse[2][0];
    
  if ( Math::Abs(fDet) <= fTolerance )
    return false;
    
  Real fInvDet = 1.0f/fDet;
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      rkInverse[iRow][iCol] *= fInvDet;
  }
  
  return true;
}
//-----------------------------------------------------------------------
Matrix3 Matrix3::Inverse (Real fTolerance) const
{
  Matrix3 kInverse = Matrix3::ZERO;
  Inverse(kInverse,fTolerance);
  return kInverse;
}
//-----------------------------------------------------------------------
Real Matrix3::Determinant () const
{
  Real fCofactor00 = m[1][1]*m[2][2] -
                     m[1][2]*m[2][1];
  Real fCofactor10 = m[1][2]*m[2][0] -
                     m[1][0]*m[2][2];
  Real fCofactor20 = m[1][0]*m[2][1] -
                     m[1][1]*m[2][0];
                     
  Real fDet =
    m[0][0]*fCofactor00 +
    m[0][1]*fCofactor10 +
    m[0][2]*fCofactor20;
    
  return fDet;
}
//-----------------------------------------------------------------------
void Matrix3::Bidiagonalize (Matrix3& kA, Matrix3& kL,
                             Matrix3& kR)
{
  Real afV[3], afW[3];
  Real fLength, fSign, fT1, fInvT1, fT2;
  bool bIdentity;
  
  // map first column to (*,0,0)
  fLength = Math::Sqrt(kA[0][0]*kA[0][0] + kA[1][0]*kA[1][0] +
                       kA[2][0]*kA[2][0]);
  if ( fLength > 0.0 )
  {
    fSign = (kA[0][0] > 0.0f ? 1.0f : -1.0f);
    fT1 = kA[0][0] + fSign*fLength;
    fInvT1 = 1.0f/fT1;
    afV[1] = kA[1][0]*fInvT1;
    afV[2] = kA[2][0]*fInvT1;
    
    fT2 = -2.0f/(1.0f+afV[1]*afV[1]+afV[2]*afV[2]);
    afW[0] = fT2*(kA[0][0]+kA[1][0]*afV[1]+kA[2][0]*afV[2]);
    afW[1] = fT2*(kA[0][1]+kA[1][1]*afV[1]+kA[2][1]*afV[2]);
    afW[2] = fT2*(kA[0][2]+kA[1][2]*afV[1]+kA[2][2]*afV[2]);
    kA[0][0] += afW[0];
    kA[0][1] += afW[1];
    kA[0][2] += afW[2];
    kA[1][1] += afV[1]*afW[1];
    kA[1][2] += afV[1]*afW[2];
    kA[2][1] += afV[2]*afW[1];
    kA[2][2] += afV[2]*afW[2];
    
    kL[0][0] = 1.0f+fT2;
    kL[0][1] = kL[1][0] = fT2*afV[1];
    kL[0][2] = kL[2][0] = fT2*afV[2];
    kL[1][1] = 1.0f+fT2*afV[1]*afV[1];
    kL[1][2] = kL[2][1] = fT2*afV[1]*afV[2];
    kL[2][2] = 1.0f+fT2*afV[2]*afV[2];
    bIdentity = false;
  }
  else
  {
    kL = Matrix3::IDENTITY;
    bIdentity = true;
  }
  
  // map first row to (*,*,0)
  fLength = Math::Sqrt(kA[0][1]*kA[0][1]+kA[0][2]*kA[0][2]);
  if ( fLength > 0.0 )
  {
    fSign = (kA[0][1] > 0.0f ? 1.0f : -1.0f);
    fT1 = kA[0][1] + fSign*fLength;
    afV[2] = kA[0][2]/fT1;
    
    fT2 = -2.0f/(1.0f+afV[2]*afV[2]);
    afW[0] = fT2*(kA[0][1]+kA[0][2]*afV[2]);
    afW[1] = fT2*(kA[1][1]+kA[1][2]*afV[2]);
    afW[2] = fT2*(kA[2][1]+kA[2][2]*afV[2]);
    kA[0][1] += afW[0];
    kA[1][1] += afW[1];
    kA[1][2] += afW[1]*afV[2];
    kA[2][1] += afW[2];
    kA[2][2] += afW[2]*afV[2];
    
    kR[0][0] = 1.0;
    kR[0][1] = kR[1][0] = 0.0;
    kR[0][2] = kR[2][0] = 0.0;
    kR[1][1] = 1.0f+fT2;
    kR[1][2] = kR[2][1] = fT2*afV[2];
    kR[2][2] = 1.0f+fT2*afV[2]*afV[2];
  }
  else
  {
    kR = Matrix3::IDENTITY;
  }
  
  // map second column to (*,*,0)
  fLength = Math::Sqrt(kA[1][1]*kA[1][1]+kA[2][1]*kA[2][1]);
  if ( fLength > 0.0 )
  {
    fSign = (kA[1][1] > 0.0f ? 1.0f : -1.0f);
    fT1 = kA[1][1] + fSign*fLength;
    afV[2] = kA[2][1]/fT1;
    
    fT2 = -2.0f/(1.0f+afV[2]*afV[2]);
    afW[1] = fT2*(kA[1][1]+kA[2][1]*afV[2]);
    afW[2] = fT2*(kA[1][2]+kA[2][2]*afV[2]);
    kA[1][1] += afW[1];
    kA[1][2] += afW[2];
    kA[2][2] += afV[2]*afW[2];
    
    Real fA = 1.0f+fT2;
    Real fB = fT2*afV[2];
    Real fC = 1.0f+fB*afV[2];
    
    if ( bIdentity )
    {
      kL[0][0] = 1.0;
      kL[0][1] = kL[1][0] = 0.0;
      kL[0][2] = kL[2][0] = 0.0;
      kL[1][1] = fA;
      kL[1][2] = kL[2][1] = fB;
      kL[2][2] = fC;
    }
    else
    {
      for (int iRow = 0; iRow < 3; iRow++)
      {
        Real fTmp0 = kL[iRow][1];
        Real fTmp1 = kL[iRow][2];
        kL[iRow][1] = fA*fTmp0+fB*fTmp1;
        kL[iRow][2] = fB*fTmp0+fC*fTmp1;
      }
    }
  }
}
//-----------------------------------------------------------------------
void Matrix3::GolubKahanStep (Matrix3& kA, Matrix3& kL,
                              Matrix3& kR)
{
  Real fT11 = kA[0][1]*kA[0][1]+kA[1][1]*kA[1][1];
  Real fT22 = kA[1][2]*kA[1][2]+kA[2][2]*kA[2][2];
  Real fT12 = kA[1][1]*kA[1][2];
  Real fTrace = fT11+fT22;
  Real fDiff = fT11-fT22;
  Real fDiscr = Math::Sqrt(fDiff*fDiff+4.0f*fT12*fT12);
  Real fRoot1 = 0.5f*(fTrace+fDiscr);
  Real fRoot2 = 0.5f*(fTrace-fDiscr);
  
  // adjust right
  Real fY = kA[0][0] - (Math::Abs(fRoot1-fT22) <=
                        Math::Abs(fRoot2-fT22) ? fRoot1 : fRoot2);
  Real fZ = kA[0][1];
  Real fInvLength = Math::InvSqrt(fY*fY+fZ*fZ);
  Real fSin = fZ*fInvLength;
  Real fCos = -fY*fInvLength;
  
  Real fTmp0 = kA[0][0];
  Real fTmp1 = kA[0][1];
  kA[0][0] = fCos*fTmp0-fSin*fTmp1;
  kA[0][1] = fSin*fTmp0+fCos*fTmp1;
  kA[1][0] = -fSin*kA[1][1];
  kA[1][1] *= fCos;
  
  size_t iRow;
  for (iRow = 0; iRow < 3; iRow++)
  {
    fTmp0 = kR[0][iRow];
    fTmp1 = kR[1][iRow];
    kR[0][iRow] = fCos*fTmp0-fSin*fTmp1;
    kR[1][iRow] = fSin*fTmp0+fCos*fTmp1;
  }
  
  // adjust left
  fY = kA[0][0];
  fZ = kA[1][0];
  fInvLength = Math::InvSqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;
  
  kA[0][0] = fCos*kA[0][0]-fSin*kA[1][0];
  fTmp0 = kA[0][1];
  fTmp1 = kA[1][1];
  kA[0][1] = fCos*fTmp0-fSin*fTmp1;
  kA[1][1] = fSin*fTmp0+fCos*fTmp1;
  kA[0][2] = -fSin*kA[1][2];
  kA[1][2] *= fCos;
  
  size_t iCol;
  for (iCol = 0; iCol < 3; iCol++)
  {
    fTmp0 = kL[iCol][0];
    fTmp1 = kL[iCol][1];
    kL[iCol][0] = fCos*fTmp0-fSin*fTmp1;
    kL[iCol][1] = fSin*fTmp0+fCos*fTmp1;
  }
  
  // adjust right
  fY = kA[0][1];
  fZ = kA[0][2];
  fInvLength = Math::InvSqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;
  
  kA[0][1] = fCos*kA[0][1]-fSin*kA[0][2];
  fTmp0 = kA[1][1];
  fTmp1 = kA[1][2];
  kA[1][1] = fCos*fTmp0-fSin*fTmp1;
  kA[1][2] = fSin*fTmp0+fCos*fTmp1;
  kA[2][1] = -fSin*kA[2][2];
  kA[2][2] *= fCos;
  
  for (iRow = 0; iRow < 3; iRow++)
  {
    fTmp0 = kR[1][iRow];
    fTmp1 = kR[2][iRow];
    kR[1][iRow] = fCos*fTmp0-fSin*fTmp1;
    kR[2][iRow] = fSin*fTmp0+fCos*fTmp1;
  }
  
  // adjust left
  fY = kA[1][1];
  fZ = kA[2][1];
  fInvLength = Math::InvSqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;
  
  kA[1][1] = fCos*kA[1][1]-fSin*kA[2][1];
  fTmp0 = kA[1][2];
  fTmp1 = kA[2][2];
  kA[1][2] = fCos*fTmp0-fSin*fTmp1;
  kA[2][2] = fSin*fTmp0+fCos*fTmp1;
  
  for (iCol = 0; iCol < 3; iCol++)
  {
    fTmp0 = kL[iCol][1];
    fTmp1 = kL[iCol][2];
    kL[iCol][1] = fCos*fTmp0-fSin*fTmp1;
    kL[iCol][2] = fSin*fTmp0+fCos*fTmp1;
  }
}
//-----------------------------------------------------------------------
void Matrix3::SingularValueDecomposition (Matrix3& kL, Vector3& kS,
    Matrix3& kR) const
{
  // temas: currently unused
  //const int iMax = 16;
  size_t iRow, iCol;
  
  Matrix3 kA = *this;
  Bidiagonalize(kA,kL,kR);
  
  for (unsigned int i = 0; i < ms_iSvdMaxIterations; i++)
  {
    Real fTmp, fTmp0, fTmp1;
    Real fSin0, fCos0, fTan0;
    Real fSin1, fCos1, fTan1;
    
    bool bTest1 = (Math::Abs(kA[0][1]) <=
                   ms_fSvdEpsilon*(Math::Abs(kA[0][0])+Math::Abs(kA[1][1])));
    bool bTest2 = (Math::Abs(kA[1][2]) <=
                   ms_fSvdEpsilon*(Math::Abs(kA[1][1])+Math::Abs(kA[2][2])));
    if ( bTest1 )
    {
      if ( bTest2 )
      {
        kS[0] = kA[0][0];
        kS[1] = kA[1][1];
        kS[2] = kA[2][2];
        break;
      }
      else
      {
        // 2x2 closed form factorization
        fTmp = (kA[1][1]*kA[1][1] - kA[2][2]*kA[2][2] +
                kA[1][2]*kA[1][2])/(kA[1][2]*kA[2][2]);
        fTan0 = 0.5f*(fTmp+Math::Sqrt(fTmp*fTmp + 4.0f));
        fCos0 = Math::InvSqrt(1.0f+fTan0*fTan0);
        fSin0 = fTan0*fCos0;
        
        for (iCol = 0; iCol < 3; iCol++)
        {
          fTmp0 = kL[iCol][1];
          fTmp1 = kL[iCol][2];
          kL[iCol][1] = fCos0*fTmp0-fSin0*fTmp1;
          kL[iCol][2] = fSin0*fTmp0+fCos0*fTmp1;
        }
        
        fTan1 = (kA[1][2]-kA[2][2]*fTan0)/kA[1][1];
        fCos1 = Math::InvSqrt(1.0f+fTan1*fTan1);
        fSin1 = -fTan1*fCos1;
        
        for (iRow = 0; iRow < 3; iRow++)
        {
          fTmp0 = kR[1][iRow];
          fTmp1 = kR[2][iRow];
          kR[1][iRow] = fCos1*fTmp0-fSin1*fTmp1;
          kR[2][iRow] = fSin1*fTmp0+fCos1*fTmp1;
        }
        
        kS[0] = kA[0][0];
        kS[1] = fCos0*fCos1*kA[1][1] -
                fSin1*(fCos0*kA[1][2]-fSin0*kA[2][2]);
        kS[2] = fSin0*fSin1*kA[1][1] +
                fCos1*(fSin0*kA[1][2]+fCos0*kA[2][2]);
        break;
      }
    }
    else
    {
      if ( bTest2 )
      {
        // 2x2 closed form factorization
        fTmp = (kA[0][0]*kA[0][0] + kA[1][1]*kA[1][1] -
                kA[0][1]*kA[0][1])/(kA[0][1]*kA[1][1]);
        fTan0 = 0.5f*(-fTmp+Math::Sqrt(fTmp*fTmp + 4.0f));
        fCos0 = Math::InvSqrt(1.0f+fTan0*fTan0);
        fSin0 = fTan0*fCos0;
        
        for (iCol = 0; iCol < 3; iCol++)
        {
          fTmp0 = kL[iCol][0];
          fTmp1 = kL[iCol][1];
          kL[iCol][0] = fCos0*fTmp0-fSin0*fTmp1;
          kL[iCol][1] = fSin0*fTmp0+fCos0*fTmp1;
        }
        
        fTan1 = (kA[0][1]-kA[1][1]*fTan0)/kA[0][0];
        fCos1 = Math::InvSqrt(1.0f+fTan1*fTan1);
        fSin1 = -fTan1*fCos1;
        
        for (iRow = 0; iRow < 3; iRow++)
        {
          fTmp0 = kR[0][iRow];
          fTmp1 = kR[1][iRow];
          kR[0][iRow] = fCos1*fTmp0-fSin1*fTmp1;
          kR[1][iRow] = fSin1*fTmp0+fCos1*fTmp1;
        }
        
        kS[0] = fCos0*fCos1*kA[0][0] -
                fSin1*(fCos0*kA[0][1]-fSin0*kA[1][1]);
        kS[1] = fSin0*fSin1*kA[0][0] +
                fCos1*(fSin0*kA[0][1]+fCos0*kA[1][1]);
        kS[2] = kA[2][2];
        break;
      }
      else
      {
        GolubKahanStep(kA,kL,kR);
      }
    }
  }
  
  // positize diagonal
  for (iRow = 0; iRow < 3; iRow++)
  {
    if ( kS[iRow] < 0.0 )
    {
      kS[iRow] = -kS[iRow];
      for (iCol = 0; iCol < 3; iCol++)
        kR[iRow][iCol] = -kR[iRow][iCol];
    }
  }
}
//-----------------------------------------------------------------------
void Matrix3::SingularValueComposition (const Matrix3& kL,
                                        const Vector3& kS, const Matrix3& kR)
{
  size_t iRow, iCol;
  Matrix3 kTmp;
  
  // product S*R
  for (iRow = 0; iRow < 3; iRow++)
  {
    for (iCol = 0; iCol < 3; iCol++)
      kTmp[iRow][iCol] = kS[iRow]*kR[iRow][iCol];
  }
  
  // product L*S*R
  for (iRow = 0; iRow < 3; iRow++)
  {
    for (iCol = 0; iCol < 3; iCol++)
    {
      m[iRow][iCol] = 0.0;
      for (int iMid = 0; iMid < 3; iMid++)
        m[iRow][iCol] += kL[iRow][iMid]*kTmp[iMid][iCol];
    }
  }
}
//-----------------------------------------------------------------------
void Matrix3::Orthonormalize ()
{
  // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
  // M = [m0|m1|m2], then orthonormal output matrix is Q = [q0|q1|q2],
  //
  //   q0 = m0/|m0|
  //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.
  
  // compute q0
  Real fInvLength = Math::InvSqrt(m[0][0]*m[0][0]
                                  + m[1][0]*m[1][0] +
                                  m[2][0]*m[2][0]);
                                  
  m[0][0] *= fInvLength;
  m[1][0] *= fInvLength;
  m[2][0] *= fInvLength;
  
  // compute q1
  Real fDot0 =
    m[0][0]*m[0][1] +
    m[1][0]*m[1][1] +
    m[2][0]*m[2][1];
    
  m[0][1] -= fDot0*m[0][0];
  m[1][1] -= fDot0*m[1][0];
  m[2][1] -= fDot0*m[2][0];
  
  fInvLength = Math::InvSqrt(m[0][1]*m[0][1] +
                             m[1][1]*m[1][1] +
                             m[2][1]*m[2][1]);
                             
  m[0][1] *= fInvLength;
  m[1][1] *= fInvLength;
  m[2][1] *= fInvLength;
  
  // compute q2
  Real fDot1 =
    m[0][1]*m[0][2] +
    m[1][1]*m[1][2] +
    m[2][1]*m[2][2];
    
  fDot0 =
    m[0][0]*m[0][2] +
    m[1][0]*m[1][2] +
    m[2][0]*m[2][2];
    
  m[0][2] -= fDot0*m[0][0] + fDot1*m[0][1];
  m[1][2] -= fDot0*m[1][0] + fDot1*m[1][1];
  m[2][2] -= fDot0*m[2][0] + fDot1*m[2][1];
  
  fInvLength = Math::InvSqrt(m[0][2]*m[0][2] +
                             m[1][2]*m[1][2] +
                             m[2][2]*m[2][2]);
                             
  m[0][2] *= fInvLength;
  m[1][2] *= fInvLength;
  m[2][2] *= fInvLength;
}
//-----------------------------------------------------------------------
void Matrix3::QDUDecomposition (Matrix3& kQ,
                                Vector3& kD, Vector3& kU) const
{
  // Factor M = QR = QDU where Q is orthogonal, D is diagonal,
  // and U is upper triangular with ones on its diagonal.  Algorithm uses
  // Gram-Schmidt orthogonalization (the QR algorithm).
  //
  // If M = [ m0 | m1 | m2 ] and Q = [ q0 | q1 | q2 ], then
  //
  //   q0 = m0/|m0|
  //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.  The matrix R has entries
  //
  //   r00 = q0*m0  r01 = q0*m1  r02 = q0*m2
  //   r10 = 0      r11 = q1*m1  r12 = q1*m2
  //   r20 = 0      r21 = 0      r22 = q2*m2
  //
  // so D = diag(r00,r11,r22) and U has entries u01 = r01/r00,
  // u02 = r02/r00, and u12 = r12/r11.
  
  // Q = rotation
  // D = scaling
  // U = shear
  
  // D stores the three diagonal entries r00, r11, r22
  // U stores the entries U[0] = u01, U[1] = u02, U[2] = u12
  
  // build orthogonal matrix Q
  Real fInvLength = Math::InvSqrt(m[0][0]*m[0][0]
                                  + m[1][0]*m[1][0] +
                                  m[2][0]*m[2][0]);
  kQ[0][0] = m[0][0]*fInvLength;
  kQ[1][0] = m[1][0]*fInvLength;
  kQ[2][0] = m[2][0]*fInvLength;
  
  Real fDot = kQ[0][0]*m[0][1] + kQ[1][0]*m[1][1] +
              kQ[2][0]*m[2][1];
  kQ[0][1] = m[0][1]-fDot*kQ[0][0];
  kQ[1][1] = m[1][1]-fDot*kQ[1][0];
  kQ[2][1] = m[2][1]-fDot*kQ[2][0];
  fInvLength = Math::InvSqrt(kQ[0][1]*kQ[0][1] + kQ[1][1]*kQ[1][1] +
                             kQ[2][1]*kQ[2][1]);
  kQ[0][1] *= fInvLength;
  kQ[1][1] *= fInvLength;
  kQ[2][1] *= fInvLength;
  
  fDot = kQ[0][0]*m[0][2] + kQ[1][0]*m[1][2] +
         kQ[2][0]*m[2][2];
  kQ[0][2] = m[0][2]-fDot*kQ[0][0];
  kQ[1][2] = m[1][2]-fDot*kQ[1][0];
  kQ[2][2] = m[2][2]-fDot*kQ[2][0];
  fDot = kQ[0][1]*m[0][2] + kQ[1][1]*m[1][2] +
         kQ[2][1]*m[2][2];
  kQ[0][2] -= fDot*kQ[0][1];
  kQ[1][2] -= fDot*kQ[1][1];
  kQ[2][2] -= fDot*kQ[2][1];
  fInvLength = Math::InvSqrt(kQ[0][2]*kQ[0][2] + kQ[1][2]*kQ[1][2] +
                             kQ[2][2]*kQ[2][2]);
  kQ[0][2] *= fInvLength;
  kQ[1][2] *= fInvLength;
  kQ[2][2] *= fInvLength;
  
  // guarantee that orthogonal matrix has determinant 1 (no reflections)
  Real fDet = kQ[0][0]*kQ[1][1]*kQ[2][2] + kQ[0][1]*kQ[1][2]*kQ[2][0] +
              kQ[0][2]*kQ[1][0]*kQ[2][1] - kQ[0][2]*kQ[1][1]*kQ[2][0] -
              kQ[0][1]*kQ[1][0]*kQ[2][2] - kQ[0][0]*kQ[1][2]*kQ[2][1];
              
  if ( fDet < 0.0 )
  {
    for (size_t iRow = 0; iRow < 3; iRow++)
      for (size_t iCol = 0; iCol < 3; iCol++)
        kQ[iRow][iCol] = -kQ[iRow][iCol];
  }
  
  // build "right" matrix R
  Matrix3 kR;
  kR[0][0] = kQ[0][0]*m[0][0] + kQ[1][0]*m[1][0] +
             kQ[2][0]*m[2][0];
  kR[0][1] = kQ[0][0]*m[0][1] + kQ[1][0]*m[1][1] +
             kQ[2][0]*m[2][1];
  kR[1][1] = kQ[0][1]*m[0][1] + kQ[1][1]*m[1][1] +
             kQ[2][1]*m[2][1];
  kR[0][2] = kQ[0][0]*m[0][2] + kQ[1][0]*m[1][2] +
             kQ[2][0]*m[2][2];
  kR[1][2] = kQ[0][1]*m[0][2] + kQ[1][1]*m[1][2] +
             kQ[2][1]*m[2][2];
  kR[2][2] = kQ[0][2]*m[0][2] + kQ[1][2]*m[1][2] +
             kQ[2][2]*m[2][2];
             
  // the scaling component
  kD[0] = kR[0][0];
  kD[1] = kR[1][1];
  kD[2] = kR[2][2];
  
  // the shear component
  Real fInvD0 = 1.0f/kD[0];
  kU[0] = kR[0][1]*fInvD0;
  kU[1] = kR[0][2]*fInvD0;
  kU[2] = kR[1][2]/kD[1];
}
//-----------------------------------------------------------------------
Real Matrix3::MaxCubicRoot (Real afCoeff[3])
{
  // Spectral norm is for A^T*A, so characteristic polynomial
  // P(x) = c[0]+c[1]*x+c[2]*x^2+x^3 has three positive real roots.
  // This yields the assertions c[0] < 0 and c[2]*c[2] >= 3*c[1].
  
  // quick out for uniform scale (triple root)
  const Real fOneThird = 1.0f/3.0f;
  const Real fEpsilon = 1e-06f;
  Real fDiscr = afCoeff[2]*afCoeff[2] - 3.0f*afCoeff[1];
  if ( fDiscr <= fEpsilon )
    return -fOneThird*afCoeff[2];
    
  // Compute an upper bound on roots of P(x).  This assumes that A^T*A
  // has been scaled by its largest entry.
  Real fX = 1.0;
  Real fPoly = afCoeff[0]+fX*(afCoeff[1]+fX*(afCoeff[2]+fX));
  if ( fPoly < 0.0 )
  {
    // uses a matrix norm to find an upper bound on maximum root
    fX = Math::Abs(afCoeff[0]);
    Real fTmp = 1.0f+Math::Abs(afCoeff[1]);
    if ( fTmp > fX )
      fX = fTmp;
    fTmp = 1.0f+Math::Abs(afCoeff[2]);
    if ( fTmp > fX )
      fX = fTmp;
  }
  
  // Newton's method to find root
  Real fTwoC2 = 2.0f*afCoeff[2];
  for (int i = 0; i < 16; i++)
  {
    fPoly = afCoeff[0]+fX*(afCoeff[1]+fX*(afCoeff[2]+fX));
    if ( Math::Abs(fPoly) <= fEpsilon )
      return fX;
      
    Real fDeriv = afCoeff[1]+fX*(fTwoC2+3.0f*fX);
    fX -= fPoly/fDeriv;
  }
  
  return fX;
}
//-----------------------------------------------------------------------
Real Matrix3::SpectralNorm () const
{
  Matrix3 kP;
  size_t iRow, iCol;
  Real fPmax = 0.0;
  for (iRow = 0; iRow < 3; iRow++)
  {
    for (iCol = 0; iCol < 3; iCol++)
    {
      kP[iRow][iCol] = 0.0;
      for (int iMid = 0; iMid < 3; iMid++)
      {
        kP[iRow][iCol] +=
          m[iMid][iRow]*m[iMid][iCol];
      }
      if ( kP[iRow][iCol] > fPmax )
        fPmax = kP[iRow][iCol];
    }
  }
  
  Real fInvPmax = 1.0f/fPmax;
  for (iRow = 0; iRow < 3; iRow++)
  {
    for (iCol = 0; iCol < 3; iCol++)
      kP[iRow][iCol] *= fInvPmax;
  }
  
  Real afCoeff[3];
  afCoeff[0] = -(kP[0][0]*(kP[1][1]*kP[2][2]-kP[1][2]*kP[2][1]) +
                 kP[0][1]*(kP[2][0]*kP[1][2]-kP[1][0]*kP[2][2]) +
                 kP[0][2]*(kP[1][0]*kP[2][1]-kP[2][0]*kP[1][1]));
  afCoeff[1] = kP[0][0]*kP[1][1]-kP[0][1]*kP[1][0] +
               kP[0][0]*kP[2][2]-kP[0][2]*kP[2][0] +
               kP[1][1]*kP[2][2]-kP[1][2]*kP[2][1];
  afCoeff[2] = -(kP[0][0]+kP[1][1]+kP[2][2]);
  
  Real fRoot = MaxCubicRoot(afCoeff);
  Real fNorm = Math::Sqrt(fPmax*fRoot);
  return fNorm;
}
//-----------------------------------------------------------------------
void Matrix3::ToAxisAngle (Vector3& rkAxis, Radian& rfRadians) const
{
  // Let (x,y,z) be the unit-length axis and let A be an angle of rotation.
  // The rotation matrix is R = I + sin(A)*P + (1-cos(A))*P^2 where
  // I is the identity and
  //
  //       +-        -+
  //   P = |  0 -z +y |
  //       | +z  0 -x |
  //       | -y +x  0 |
  //       +-        -+
  //
  // If A > 0, R represents a counterclockwise rotation about the axis in
  // the sense of looking from the tip of the axis vector towards the
  // origin.  Some algebra will show that
  //
  //   cos(A) = (trace(R)-1)/2  and  R - R^t = 2*sin(A)*P
  //
  // In the event that A = pi, R-R^t = 0 which prevents us from extracting
  // the axis through P.  Instead note that R = I+2*P^2 when A = pi, so
  // P^2 = (R-I)/2.  The diagonal entries of P^2 are x^2-1, y^2-1, and
  // z^2-1.  We can solve these for axis (x,y,z).  Because the angle is pi,
  // it does not matter which sign you choose on the square roots.
  
  Real fTrace = m[0][0] + m[1][1] + m[2][2];
  Real fCos = 0.5f*(fTrace-1.0f);
  rfRadians = Math::ACos(fCos);  // in [0,PI]
  
  if ( rfRadians > Radian(0.0) )
  {
    if ( rfRadians < Radian(Math::PI) )
    {
      rkAxis.x = m[2][1]-m[1][2];
      rkAxis.y = m[0][2]-m[2][0];
      rkAxis.z = m[1][0]-m[0][1];
      rkAxis.normalise();
    }
    else
    {
      // angle is PI
      float fHalfInverse;
      if ( m[0][0] >= m[1][1] )
      {
        // r00 >= r11
        if ( m[0][0] >= m[2][2] )
        {
          // r00 is maximum diagonal term
          rkAxis.x = 0.5f*Math::Sqrt(m[0][0] -
                                     m[1][1] - m[2][2] + 1.0f);
          fHalfInverse = 0.5f/rkAxis.x;
          rkAxis.y = fHalfInverse*m[0][1];
          rkAxis.z = fHalfInverse*m[0][2];
        }
        else
        {
          // r22 is maximum diagonal term
          rkAxis.z = 0.5f*Math::Sqrt(m[2][2] -
                                     m[0][0] - m[1][1] + 1.0f);
          fHalfInverse = 0.5f/rkAxis.z;
          rkAxis.x = fHalfInverse*m[0][2];
          rkAxis.y = fHalfInverse*m[1][2];
        }
      }
      else
      {
        // r11 > r00
        if ( m[1][1] >= m[2][2] )
        {
          // r11 is maximum diagonal term
          rkAxis.y = 0.5f*Math::Sqrt(m[1][1] -
                                     m[0][0] - m[2][2] + 1.0f);
          fHalfInverse  = 0.5f/rkAxis.y;
          rkAxis.x = fHalfInverse*m[0][1];
          rkAxis.z = fHalfInverse*m[1][2];
        }
        else
        {
          // r22 is maximum diagonal term
          rkAxis.z = 0.5f*Math::Sqrt(m[2][2] -
                                     m[0][0] - m[1][1] + 1.0f);
          fHalfInverse = 0.5f/rkAxis.z;
          rkAxis.x = fHalfInverse*m[0][2];
          rkAxis.y = fHalfInverse*m[1][2];
        }
      }
    }
  }
  else
  {
    // The angle is 0 and the matrix is the identity.  Any axis will
    // work, so just use the x-axis.
    rkAxis.x = 1.0;
    rkAxis.y = 0.0;
    rkAxis.z = 0.0;
  }
}
//-----------------------------------------------------------------------
void Matrix3::FromAxisAngle (const Vector3& rkAxis, const Radian& fRadians)
{
  Real fCos = Math::Cos(fRadians);
  Real fSin = Math::Sin(fRadians);
  Real fOneMinusCos = 1.0f-fCos;
  Real fX2 = rkAxis.x*rkAxis.x;
  Real fY2 = rkAxis.y*rkAxis.y;
  Real fZ2 = rkAxis.z*rkAxis.z;
  Real fXYM = rkAxis.x*rkAxis.y*fOneMinusCos;
  Real fXZM = rkAxis.x*rkAxis.z*fOneMinusCos;
  Real fYZM = rkAxis.y*rkAxis.z*fOneMinusCos;
  Real fXSin = rkAxis.x*fSin;
  Real fYSin = rkAxis.y*fSin;
  Real fZSin = rkAxis.z*fSin;
  
  m[0][0] = fX2*fOneMinusCos+fCos;
  m[0][1] = fXYM-fZSin;
  m[0][2] = fXZM+fYSin;
  m[1][0] = fXYM+fZSin;
  m[1][1] = fY2*fOneMinusCos+fCos;
  m[1][2] = fYZM-fXSin;
  m[2][0] = fXZM-fYSin;
  m[2][1] = fYZM+fXSin;
  m[2][2] = fZ2*fOneMinusCos+fCos;
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesXYZ (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz          -cy*sz           sy
  //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
  //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy
  
  rfPAngle = Radian(Math::ASin(m[0][2]));
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(-m[1][2],m[2][2]);
      rfRAngle = Math::ATan2(-m[0][1],m[0][0]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(m[1][0],m[1][1]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(m[1][0],m[1][1]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesXZY (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz          -sz              cz*sy
  //        sx*sy+cx*cy*sz  cx*cz          -cy*sx+cx*sy*sz
  //       -cx*sy+cy*sx*sz  cz*sx           cx*cy+sx*sy*sz
  
  rfPAngle = Math::ASin(-m[0][1]);
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(m[2][1],m[1][1]);
      rfRAngle = Math::ATan2(m[0][2],m[0][0]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(-m[2][0],m[2][2]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(-m[2][0],m[2][2]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesYXZ (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz+sx*sy*sz  cz*sx*sy-cy*sz  cx*sy
  //        cx*sz           cx*cz          -sx
  //       -cz*sy+cy*sx*sz  cy*cz*sx+sy*sz  cx*cy
  
  rfPAngle = Math::ASin(-m[1][2]);
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(m[0][2],m[2][2]);
      rfRAngle = Math::ATan2(m[1][0],m[1][1]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(-m[0][1],m[0][0]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(-m[0][1],m[0][0]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesYZX (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz           sx*sy-cx*cy*sz  cx*sy+cy*sx*sz
  //        sz              cx*cz          -cz*sx
  //       -cz*sy           cy*sx+cx*sy*sz  cx*cy-sx*sy*sz
  
  rfPAngle = Math::ASin(m[1][0]);
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(-m[2][0],m[0][0]);
      rfRAngle = Math::ATan2(-m[1][2],m[1][1]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(m[2][1],m[2][2]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(m[2][1],m[2][2]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesZXY (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
  //        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
  //       -cx*sy           sx              cx*cy
  
  rfPAngle = Math::ASin(m[2][1]);
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(-m[0][1],m[1][1]);
      rfRAngle = Math::ATan2(-m[2][0],m[2][2]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(m[0][2],m[0][0]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(m[0][2],m[0][0]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::ToEulerAnglesZYX (Radian& rfYAngle, Radian& rfPAngle,
                                Radian& rfRAngle) const
{
  // rot =  cy*cz           cz*sx*sy-cx*sz  cx*cz*sy+sx*sz
  //        cy*sz           cx*cz+sx*sy*sz -cz*sx+cx*sy*sz
  //       -sy              cy*sx           cx*cy
  
  rfPAngle = Math::ASin(-m[2][0]);
  if ( rfPAngle < Radian(Math::HALF_PI) )
  {
    if ( rfPAngle > Radian(-Math::HALF_PI) )
    {
      rfYAngle = Math::ATan2(m[1][0],m[0][0]);
      rfRAngle = Math::ATan2(m[2][1],m[2][2]);
      return true;
    }
    else
    {
      // WARNING.  Not a unique solution.
      Radian fRmY = Math::ATan2(-m[0][1],m[0][2]);
      rfRAngle = Radian(0.0);  // any angle works
      rfYAngle = rfRAngle - fRmY;
      return false;
    }
  }
  else
  {
    // WARNING.  Not a unique solution.
    Radian fRpY = Math::ATan2(-m[0][1],m[0][2]);
    rfRAngle = Radian(0.0);  // any angle works
    rfYAngle = fRpY - rfRAngle;
    return false;
  }
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesXYZ (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  *this = kXMat*(kYMat*kZMat);
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesXZY (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  *this = kXMat*(kZMat*kYMat);
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesYXZ (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  *this = kYMat*(kXMat*kZMat);
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesYZX (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  *this = kYMat*(kZMat*kXMat);
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesZXY (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  *this = kZMat*(kXMat*kYMat);
}
//-----------------------------------------------------------------------
void Matrix3::FromEulerAnglesZYX (const Radian& fYAngle, const Radian& fPAngle,
                                  const Radian& fRAngle)
{
  Real fCos, fSin;
  
  fCos = Math::Cos(fYAngle);
  fSin = Math::Sin(fYAngle);
  Matrix3 kZMat(fCos,-fSin,0.0,fSin,fCos,0.0,0.0,0.0,1.0);
  
  fCos = Math::Cos(fPAngle);
  fSin = Math::Sin(fPAngle);
  Matrix3 kYMat(fCos,0.0,fSin,0.0,1.0,0.0,-fSin,0.0,fCos);
  
  fCos = Math::Cos(fRAngle);
  fSin = Math::Sin(fRAngle);
  Matrix3 kXMat(1.0,0.0,0.0,0.0,fCos,-fSin,0.0,fSin,fCos);
  
  *this = kZMat*(kYMat*kXMat);
}
//-----------------------------------------------------------------------
void Matrix3::Tridiagonal (Real afDiag[3], Real afSubDiag[3])
{
  // Householder reduction T = Q^t M Q
  //   Input:
  //     mat, symmetric 3x3 matrix M
  //   Output:
  //     mat, orthogonal matrix Q
  //     diag, diagonal entries of T
  //     subd, subdiagonal entries of T (T is symmetric)
  
  Real fA = m[0][0];
  Real fB = m[0][1];
  Real fC = m[0][2];
  Real fD = m[1][1];
  Real fE = m[1][2];
  Real fF = m[2][2];
  
  afDiag[0] = fA;
  afSubDiag[2] = 0.0;
  if ( Math::Abs(fC) >= EPSILON )
  {
    Real fLength = Math::Sqrt(fB*fB+fC*fC);
    Real fInvLength = 1.0f/fLength;
    fB *= fInvLength;
    fC *= fInvLength;
    Real fQ = 2.0f*fB*fE+fC*(fF-fD);
    afDiag[1] = fD+fC*fQ;
    afDiag[2] = fF-fC*fQ;
    afSubDiag[0] = fLength;
    afSubDiag[1] = fE-fB*fQ;
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = fB;
    m[1][2] = fC;
    m[2][0] = 0.0;
    m[2][1] = fC;
    m[2][2] = -fB;
  }
  else
  {
    afDiag[1] = fD;
    afDiag[2] = fF;
    afSubDiag[0] = fB;
    afSubDiag[1] = fE;
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
  }
}
//-----------------------------------------------------------------------
bool Matrix3::QLAlgorithm (Real afDiag[3], Real afSubDiag[3])
{
  // QL iteration with implicit shifting to reduce matrix from tridiagonal
  // to diagonal
  
  for (int i0 = 0; i0 < 3; i0++)
  {
    const unsigned int iMaxIter = 32;
    unsigned int iIter;
    for (iIter = 0; iIter < iMaxIter; iIter++)
    {
      int i1;
      for (i1 = i0; i1 <= 1; i1++)
      {
        Real fSum = Math::Abs(afDiag[i1]) +
                    Math::Abs(afDiag[i1+1]);
        if ( Math::Abs(afSubDiag[i1]) + fSum == fSum )
          break;
      }
      if ( i1 == i0 )
        break;
        
      Real fTmp0 = (afDiag[i0+1]-afDiag[i0])/(2.0f*afSubDiag[i0]);
      Real fTmp1 = Math::Sqrt(fTmp0*fTmp0+1.0f);
      if ( fTmp0 < 0.0 )
        fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0-fTmp1);
      else
        fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0+fTmp1);
      Real fSin = 1.0;
      Real fCos = 1.0;
      Real fTmp2 = 0.0;
      for (int i2 = i1-1; i2 >= i0; i2--)
      {
        Real fTmp3 = fSin*afSubDiag[i2];
        Real fTmp4 = fCos*afSubDiag[i2];
        if ( Math::Abs(fTmp3) >= Math::Abs(fTmp0) )
        {
          fCos = fTmp0/fTmp3;
          fTmp1 = Math::Sqrt(fCos*fCos+1.0f);
          afSubDiag[i2+1] = fTmp3*fTmp1;
          fSin = 1.0f/fTmp1;
          fCos *= fSin;
        }
        else
        {
          fSin = fTmp3/fTmp0;
          fTmp1 = Math::Sqrt(fSin*fSin+1.0f);
          afSubDiag[i2+1] = fTmp0*fTmp1;
          fCos = 1.0f/fTmp1;
          fSin *= fCos;
        }
        fTmp0 = afDiag[i2+1]-fTmp2;
        fTmp1 = (afDiag[i2]-fTmp0)*fSin+2.0f*fTmp4*fCos;
        fTmp2 = fSin*fTmp1;
        afDiag[i2+1] = fTmp0+fTmp2;
        fTmp0 = fCos*fTmp1-fTmp4;
        
        for (int iRow = 0; iRow < 3; iRow++)
        {
          fTmp3 = m[iRow][i2+1];
          m[iRow][i2+1] = fSin*m[iRow][i2] +
                          fCos*fTmp3;
          m[iRow][i2] = fCos*m[iRow][i2] -
                        fSin*fTmp3;
        }
      }
      afDiag[i0] -= fTmp2;
      afSubDiag[i0] = fTmp0;
      afSubDiag[i1] = 0.0;
    }
    
    if ( iIter == iMaxIter )
    {
      // should not get here under normal circumstances
      return false;
    }
  }
  
  return true;
}
//-----------------------------------------------------------------------
void Matrix3::EigenSolveSymmetric (Real afEigenvalue[3],
                                   Vector3 akEigenvector[3]) const
{
  Matrix3 kMatrix = *this;
  Real afSubDiag[3];
  kMatrix.Tridiagonal(afEigenvalue,afSubDiag);
  kMatrix.QLAlgorithm(afEigenvalue,afSubDiag);
  
  for (size_t i = 0; i < 3; i++)
  {
    akEigenvector[i][0] = kMatrix[0][i];
    akEigenvector[i][1] = kMatrix[1][i];
    akEigenvector[i][2] = kMatrix[2][i];
  }
  
  // make eigenvectors form a right--handed system
  Vector3 kCross = akEigenvector[1].crossProduct(akEigenvector[2]);
  Real fDet = akEigenvector[0].dotProduct(kCross);
  if ( fDet < 0.0 )
  {
    akEigenvector[2][0] = - akEigenvector[2][0];
    akEigenvector[2][1] = - akEigenvector[2][1];
    akEigenvector[2][2] = - akEigenvector[2][2];
  }
}
//-----------------------------------------------------------------------
void Matrix3::TensorProduct (const Vector3& rkU, const Vector3& rkV,
                             Matrix3& rkProduct)
{
  for (size_t iRow = 0; iRow < 3; iRow++)
  {
    for (size_t iCol = 0; iCol < 3; iCol++)
      rkProduct[iRow][iCol] = rkU[iRow]*rkV[iCol];
  }
}
//-----------------------------------------------------------------------

const Matrix4 Matrix4::ZERO(
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0 );
  
const Matrix4 Matrix4::IDENTITY(
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1 );
  
const Matrix4 Matrix4::CLIPSPACE2DTOIMAGESPACE(
  0.5,    0,  0, 0.5,
  0, -0.5,  0, 0.5,
  0,    0,  1,   0,
  0,    0,  0,   1);
  
//-----------------------------------------------------------------------
inline static Real
MINOR(const Matrix4& m, const size_t r0, const size_t r1, const size_t r2,
      const size_t c0, const size_t c1, const size_t c2)
{
  return m[r0][c0] * (m[r1][c1] * m[r2][c2] - m[r2][c1] * m[r1][c2]) -
         m[r0][c1] * (m[r1][c0] * m[r2][c2] - m[r2][c0] * m[r1][c2]) +
         m[r0][c2] * (m[r1][c0] * m[r2][c1] - m[r2][c0] * m[r1][c1]);
}
//-----------------------------------------------------------------------
Matrix4 Matrix4::adjoint() const
{
  return Matrix4( MINOR(*this, 1, 2, 3, 1, 2, 3),
                  -MINOR(*this, 0, 2, 3, 1, 2, 3),
                  MINOR(*this, 0, 1, 3, 1, 2, 3),
                  -MINOR(*this, 0, 1, 2, 1, 2, 3),
                  
                  -MINOR(*this, 1, 2, 3, 0, 2, 3),
                  MINOR(*this, 0, 2, 3, 0, 2, 3),
                  -MINOR(*this, 0, 1, 3, 0, 2, 3),
                  MINOR(*this, 0, 1, 2, 0, 2, 3),
                  
                  MINOR(*this, 1, 2, 3, 0, 1, 3),
                  -MINOR(*this, 0, 2, 3, 0, 1, 3),
                  MINOR(*this, 0, 1, 3, 0, 1, 3),
                  -MINOR(*this, 0, 1, 2, 0, 1, 3),
                  
                  -MINOR(*this, 1, 2, 3, 0, 1, 2),
                  MINOR(*this, 0, 2, 3, 0, 1, 2),
                  -MINOR(*this, 0, 1, 3, 0, 1, 2),
                  MINOR(*this, 0, 1, 2, 0, 1, 2));
}
//-----------------------------------------------------------------------
Real Matrix4::determinant() const
{
  return m[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
         m[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
         m[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
         m[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
}
//-----------------------------------------------------------------------
Matrix4 Matrix4::inverse() const
{
  Real m00 = m[0][0], m01 = m[0][1], m02 = m[0][2], m03 = m[0][3];
  Real m10 = m[1][0], m11 = m[1][1], m12 = m[1][2], m13 = m[1][3];
  Real m20 = m[2][0], m21 = m[2][1], m22 = m[2][2], m23 = m[2][3];
  Real m30 = m[3][0], m31 = m[3][1], m32 = m[3][2], m33 = m[3][3];
  
  Real v0 = m20 * m31 - m21 * m30;
  Real v1 = m20 * m32 - m22 * m30;
  Real v2 = m20 * m33 - m23 * m30;
  Real v3 = m21 * m32 - m22 * m31;
  Real v4 = m21 * m33 - m23 * m31;
  Real v5 = m22 * m33 - m23 * m32;
  
  Real t00 = + (v5 * m11 - v4 * m12 + v3 * m13);
  Real t10 = - (v5 * m10 - v2 * m12 + v1 * m13);
  Real t20 = + (v4 * m10 - v2 * m11 + v0 * m13);
  Real t30 = - (v3 * m10 - v1 * m11 + v0 * m12);
  
  Real invDet = 1 / (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);
  
  Real d00 = t00 * invDet;
  Real d10 = t10 * invDet;
  Real d20 = t20 * invDet;
  Real d30 = t30 * invDet;
  
  Real d01 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  Real d11 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  Real d21 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  Real d31 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;
  
  v0 = m10 * m31 - m11 * m30;
  v1 = m10 * m32 - m12 * m30;
  v2 = m10 * m33 - m13 * m30;
  v3 = m11 * m32 - m12 * m31;
  v4 = m11 * m33 - m13 * m31;
  v5 = m12 * m33 - m13 * m32;
  
  Real d02 = + (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  Real d12 = - (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  Real d22 = + (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  Real d32 = - (v3 * m00 - v1 * m01 + v0 * m02) * invDet;
  
  v0 = m21 * m10 - m20 * m11;
  v1 = m22 * m10 - m20 * m12;
  v2 = m23 * m10 - m20 * m13;
  v3 = m22 * m11 - m21 * m12;
  v4 = m23 * m11 - m21 * m13;
  v5 = m23 * m12 - m22 * m13;
  
  Real d03 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
  Real d13 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
  Real d23 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
  Real d33 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;
  
  return Matrix4(
           d00, d01, d02, d03,
           d10, d11, d12, d13,
           d20, d21, d22, d23,
           d30, d31, d32, d33);
}
//-----------------------------------------------------------------------
Matrix4 Matrix4::inverseAffine(void) const
{
  assert(isAffine());
  
  Real m10 = m[1][0], m11 = m[1][1], m12 = m[1][2];
  Real m20 = m[2][0], m21 = m[2][1], m22 = m[2][2];
  
  Real t00 = m22 * m11 - m21 * m12;
  Real t10 = m20 * m12 - m22 * m10;
  Real t20 = m21 * m10 - m20 * m11;
  
  Real m00 = m[0][0], m01 = m[0][1], m02 = m[0][2];
  
  Real invDet = 1 / (m00 * t00 + m01 * t10 + m02 * t20);
  
  t00 *= invDet;
  t10 *= invDet;
  t20 *= invDet;
  
  m00 *= invDet;
  m01 *= invDet;
  m02 *= invDet;
  
  Real r00 = t00;
  Real r01 = m02 * m21 - m01 * m22;
  Real r02 = m01 * m12 - m02 * m11;
  
  Real r10 = t10;
  Real r11 = m00 * m22 - m02 * m20;
  Real r12 = m02 * m10 - m00 * m12;
  
  Real r20 = t20;
  Real r21 = m01 * m20 - m00 * m21;
  Real r22 = m00 * m11 - m01 * m10;
  
  Real m03 = m[0][3], m13 = m[1][3], m23 = m[2][3];
  
  Real r03 = - (r00 * m03 + r01 * m13 + r02 * m23);
  Real r13 = - (r10 * m03 + r11 * m13 + r12 * m23);
  Real r23 = - (r20 * m03 + r21 * m13 + r22 * m23);
  
  return Matrix4(
           r00, r01, r02, r03,
           r10, r11, r12, r13,
           r20, r21, r22, r23,
           0,   0,   0,   1);
}
//-----------------------------------------------------------------------
void Matrix4::makeTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation)
{
  // Ordering:
  //    1. Scale
  //    2. Rotate
  //    3. Translate
  
  Matrix3 rot3x3;
  orientation.ToRotationMatrix(rot3x3);
  
  // Set up final matrix with scale, rotation and translation
  m[0][0] = scale.x * rot3x3[0][0];
  m[0][1] = scale.y * rot3x3[0][1];
  m[0][2] = scale.z * rot3x3[0][2];
  m[0][3] = position.x;
  m[1][0] = scale.x * rot3x3[1][0];
  m[1][1] = scale.y * rot3x3[1][1];
  m[1][2] = scale.z * rot3x3[1][2];
  m[1][3] = position.y;
  m[2][0] = scale.x * rot3x3[2][0];
  m[2][1] = scale.y * rot3x3[2][1];
  m[2][2] = scale.z * rot3x3[2][2];
  m[2][3] = position.z;
  
  // No projection term
  m[3][0] = 0;
  m[3][1] = 0;
  m[3][2] = 0;
  m[3][3] = 1;
}
//-----------------------------------------------------------------------
void Matrix4::makeInverseTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation)
{
  // Invert the parameters
  Vector3 invTranslate = -position;
  Vector3 invScale(1 / scale.x, 1 / scale.y, 1 / scale.z);
  Quaternion invRot = orientation.Inverse();
  
  // Because we're inverting, order is translation, rotation, scale
  // So make translation relative to scale & rotation
  invTranslate = invRot * invTranslate; // rotate
  invTranslate *= invScale; // scale
  
  // Next, make a 3x3 rotation matrix
  Matrix3 rot3x3;
  invRot.ToRotationMatrix(rot3x3);
  
  // Set up final matrix with scale, rotation and translation
  m[0][0] = invScale.x * rot3x3[0][0];
  m[0][1] = invScale.x * rot3x3[0][1];
  m[0][2] = invScale.x * rot3x3[0][2];
  m[0][3] = invTranslate.x;
  m[1][0] = invScale.y * rot3x3[1][0];
  m[1][1] = invScale.y * rot3x3[1][1];
  m[1][2] = invScale.y * rot3x3[1][2];
  m[1][3] = invTranslate.y;
  m[2][0] = invScale.z * rot3x3[2][0];
  m[2][1] = invScale.z * rot3x3[2][1];
  m[2][2] = invScale.z * rot3x3[2][2];
  m[2][3] = invTranslate.z;
  
  // No projection term
  m[3][0] = 0;
  m[3][1] = 0;
  m[3][2] = 0;
  m[3][3] = 1;
}
//-----------------------------------------------------------------------
void Matrix4::decomposition(Vector3& position, Vector3& scale, Quaternion& orientation) const
{
  assert(isAffine());
  
  Matrix3 m3x3;
  extract3x3Matrix(m3x3);
  
  Matrix3 matQ;
  Vector3 vecU;
  m3x3.QDUDecomposition( matQ, scale, vecU );
  
  orientation = Quaternion( matQ );
  position = Vector3( m[0][3], m[1][3], m[2][3] );
}

} // namespace OgreMath
