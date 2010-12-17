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


#ifndef OGREMATHLIB_MATH_LIB_H
#define OGREMATHLIB_MATH_LIB_H

#include <string>
#include <ostream>
#include <sstream>
#include "assert.h"

#define OGREMATHLIB_EXPORT

namespace OgreMathLib
{

typedef float Real;
class Vector3;
class Quaternion;
class Matrix3;
class Matrix4;
class Degree;
class Radian;
class Math;

/** Wrapper class which indicates a given angle value is in Radians.
   @remarks
       Radian values are interchangeable with Degree values, and conversions
       will be done automatically between them.
*/
class OGREMATHLIB_EXPORT Radian
{
  Real mRad;
  
public:
  explicit Radian ( Real r=0 ) : mRad(r) {}
  Radian ( const Degree& d );
  Radian& operator = ( const Real& f ) {
    mRad = f;
    return *this;
  }
  Radian& operator = ( const Radian& r ) {
    mRad = r.mRad;
    return *this;
  }
  Radian& operator = ( const Degree& d );
  
  Real valueDegrees() const; // see bottom of this file
  Real valueRadians() const {
    return mRad;
  }
  Real valueAngleUnits() const;
  
  const Radian& operator + () const {
    return *this;
  }
  Radian operator + ( const Radian& r ) const {
    return Radian ( mRad + r.mRad );
  }
  Radian operator + ( const Degree& d ) const;
  Radian& operator += ( const Radian& r ) {
    mRad += r.mRad;
    return *this;
  }
  Radian& operator += ( const Degree& d );
  Radian operator - () const {
    return Radian(-mRad);
  }
  Radian operator - ( const Radian& r ) const {
    return Radian ( mRad - r.mRad );
  }
  Radian operator - ( const Degree& d ) const;
  Radian& operator -= ( const Radian& r ) {
    mRad -= r.mRad;
    return *this;
  }
  Radian& operator -= ( const Degree& d );
  Radian operator * ( Real f ) const {
    return Radian ( mRad * f );
  }
  Radian operator * ( const Radian& f ) const {
    return Radian ( mRad * f.mRad );
  }
  Radian& operator *= ( Real f ) {
    mRad *= f;
    return *this;
  }
  Radian operator / ( Real f ) const {
    return Radian ( mRad / f );
  }
  Radian& operator /= ( Real f ) {
    mRad /= f;
    return *this;
  }
  
  bool operator <  ( const Radian& r ) const {
    return mRad <  r.mRad;
  }
  bool operator <= ( const Radian& r ) const {
    return mRad <= r.mRad;
  }
  bool operator == ( const Radian& r ) const {
    return mRad == r.mRad;
  }
  bool operator != ( const Radian& r ) const {
    return mRad != r.mRad;
  }
  bool operator >= ( const Radian& r ) const {
    return mRad >= r.mRad;
  }
  bool operator >  ( const Radian& r ) const {
    return mRad >  r.mRad;
  }
  
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Radian& v )
  {
    o << "Radian(" << v.valueRadians() << ")";
    return o;
  }
  
  inline std::string to_s() const
  {
    std::ostringstream o;
    o << mRad;
    return o.str();
  }
  
};

/** Wrapper class which indicates a given angle value is in Degrees.
@remarks
    Degree values are interchangeable with Radian values, and conversions
    will be done automatically between them.
*/
class OGREMATHLIB_EXPORT Degree
{
  Real mDeg; // if you get an error here - make sure to define/typedef 'Real' first
  
public:
  explicit Degree ( Real d=0 ) : mDeg(d) {}
  Degree ( const Radian& r ) : mDeg(r.valueDegrees()) {}
  Degree& operator = ( const Real& f ) {
    mDeg = f;
    return *this;
  }
  Degree& operator = ( const Degree& d ) {
    mDeg = d.mDeg;
    return *this;
  }
  Degree& operator = ( const Radian& r ) {
    mDeg = r.valueDegrees();
    return *this;
  }
  
  Real valueDegrees() const {
    return mDeg;
  }
  Real valueRadians() const; // see bottom of this file
  Real valueAngleUnits() const;
  
  const Degree& operator + () const {
    return *this;
  }
  Degree operator + ( const Degree& d ) const {
    return Degree ( mDeg + d.mDeg );
  }
  Degree operator + ( const Radian& r ) const {
    return Degree ( mDeg + r.valueDegrees() );
  }
  Degree& operator += ( const Degree& d ) {
    mDeg += d.mDeg;
    return *this;
  }
  Degree& operator += ( const Radian& r ) {
    mDeg += r.valueDegrees();
    return *this;
  }
  Degree operator - () const {
    return Degree(-mDeg);
  }
  Degree operator - ( const Degree& d ) const {
    return Degree ( mDeg - d.mDeg );
  }
  Degree operator - ( const Radian& r ) const {
    return Degree ( mDeg - r.valueDegrees() );
  }
  Degree& operator -= ( const Degree& d ) {
    mDeg -= d.mDeg;
    return *this;
  }
  Degree& operator -= ( const Radian& r ) {
    mDeg -= r.valueDegrees();
    return *this;
  }
  Degree operator * ( Real f ) const {
    return Degree ( mDeg * f );
  }
  Degree operator * ( const Degree& f ) const {
    return Degree ( mDeg * f.mDeg );
  }
  Degree& operator *= ( Real f ) {
    mDeg *= f;
    return *this;
  }
  Degree operator / ( Real f ) const {
    return Degree ( mDeg / f );
  }
  Degree& operator /= ( Real f ) {
    mDeg /= f;
    return *this;
  }
  
  bool operator <  ( const Degree& d ) const {
    return mDeg <  d.mDeg;
  }
  bool operator <= ( const Degree& d ) const {
    return mDeg <= d.mDeg;
  }
  bool operator == ( const Degree& d ) const {
    return mDeg == d.mDeg;
  }
  bool operator != ( const Degree& d ) const {
    return mDeg != d.mDeg;
  }
  bool operator >= ( const Degree& d ) const {
    return mDeg >= d.mDeg;
  }
  bool operator >  ( const Degree& d ) const {
    return mDeg >  d.mDeg;
  }
  
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Degree& v )
  {
    o << "Degree(" << v.valueDegrees() << ")";
    return o;
  }
  
  inline std::string to_s() const
  {
    std::ostringstream o;
    o << mDeg;
    return o.str();
  }
  
};

/** Class to provide access to common mathematical functions.
       @remarks
           Most of the maths functions are aliased versions of the C runtime
           library functions. They are aliased here to provide future
           optimisation opportunities, either from faster RTLs or custom
           math approximations.
       @note
           <br>This is based on MgcMath.h from
           <a href="http://www.geometrictools.com/">Wild Magic</a>.
   */
class OGREMATHLIB_EXPORT Math
{
public:
  /** The angular units used by the API. This functionality is now deprecated in favor
    of discreet angular unit types ( see Degree and Radian above ). The only place
  this functionality is actually still used is when parsing files. Search for
  usage of the Angle class for those instances
  */
  enum AngleUnit
  {
    AU_DEGREE,
    AU_RADIAN
  };
  
protected:


  /// Size of the trig tables as determined by constructor.
  static int mTrigTableSize;
  
  /// Radian -> index factor value ( mTrigTableSize / 2 * PI )
  static Real mTrigTableFactor;
  static Real* mSinTable;
  static Real* mTanTable;
  
  /** Private function to build trig tables.
  */
  void buildTrigTables();
  
  static Real SinTable (Real fValue);
  static Real TanTable (Real fValue);
public:
  /** Default constructor.
      @param
          trigTableSize Optional parameter to set the size of the
          tables used to implement Sin, Cos, Tan
  */
  Math(unsigned int trigTableSize = 4096);
  
  /** Default destructor.
  */
  ~Math();
  
  static inline int IAbs (int iValue) {
    return ( iValue >= 0 ? iValue : -iValue );
  }
  static inline int ICeil (float fValue) {
    return int(ceil(fValue));
  }
  static inline int IFloor (float fValue) {
    return int(floor(fValue));
  }
  static int ISign (int iValue);
  
  static inline Real Abs (Real fValue) {
    return Real(fabs(fValue));
  }
  static inline Degree Abs (const Degree& dValue) {
    return Degree(fabs(dValue.valueDegrees()));
  }
  static inline Radian Abs (const Radian& rValue) {
    return Radian(fabs(rValue.valueRadians()));
  }
  static Radian ACos (Real fValue);
  static Radian ASin (Real fValue);
  static inline Radian ATan (Real fValue) {
    return Radian(atan(fValue));
  }
  static inline Radian ATan2 (Real fY, Real fX) {
    return Radian(atan2(fY,fX));
  }
  static inline Real Ceil (Real fValue) {
    return Real(ceil(fValue));
  }
  static inline bool isNaN(Real f)
  {
    // std::isnan() is C99, not supported by all compilers
    // However NaN always fails this next test, no other number does.
    return f != f;
  }
  
  /** Cosine function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Cos (const Radian& fValue, bool useTables = false) {
    return (!useTables) ? Real(cos(fValue.valueRadians())) : SinTable(fValue.valueRadians() + HALF_PI);
  }
  /** Cosine function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Cos (Real fValue, bool useTables = false) {
    return (!useTables) ? Real(cos(fValue)) : SinTable(fValue + HALF_PI);
  }
  
  static inline Real Exp (Real fValue) {
    return Real(exp(fValue));
  }
  
  static inline Real Floor (Real fValue) {
    return Real(floor(fValue));
  }
  
  static inline Real Log (Real fValue) {
    return Real(log(fValue));
  }
  
  /// Stored value of log(2) for frequent use
  static const Real LOG2;
  
  static inline Real Log2 (Real fValue) {
    return Real(log(fValue)/LOG2);
  }
  
  static inline Real LogN (Real base, Real fValue) {
    return Real(log(fValue)/log(base));
  }
  
  static inline Real Pow (Real fBase, Real fExponent) {
    return Real(pow(fBase,fExponent));
  }
  
  static Real Sign (Real fValue);
  static inline Radian Sign ( const Radian& rValue )
  {
    return Radian(Sign(rValue.valueRadians()));
  }
  static inline Degree Sign ( const Degree& dValue )
  {
    return Degree(Sign(dValue.valueDegrees()));
  }
  
  /** Sine function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Sin (const Radian& fValue, bool useTables = false) {
    return (!useTables) ? Real(sin(fValue.valueRadians())) : SinTable(fValue.valueRadians());
  }
  /** Sine function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Sin (Real fValue, bool useTables = false) {
    return (!useTables) ? Real(sin(fValue)) : SinTable(fValue);
  }
  
  static inline Real Sqr (Real fValue) {
    return fValue*fValue;
  }
  
  static inline Real Sqrt (Real fValue) {
    return Real(sqrt(fValue));
  }
  
  static inline Radian Sqrt (const Radian& fValue) {
    return Radian(sqrt(fValue.valueRadians()));
  }
  
  static inline Degree Sqrt (const Degree& fValue) {
    return Degree(sqrt(fValue.valueDegrees()));
  }
  
  /** Inverse square root i.e. 1 / Sqrt(x), good for vector
      normalisation.
  */
  static Real InvSqrt(Real fValue);
  
  static Real UnitRandom ();  // in [0,1]
  
  static Real RangeRandom (Real fLow, Real fHigh);  // in [fLow,fHigh]
  
  static Real SymmetricRandom ();  // in [-1,1]
  
  /** Tangent function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Tan (const Radian& fValue, bool useTables = false) {
    return (!useTables) ? Real(tan(fValue.valueRadians())) : TanTable(fValue.valueRadians());
  }
  /** Tangent function.
      @param
          fValue Angle in radians
      @param
          useTables If true, uses lookup tables rather than
          calculation - faster but less accurate.
  */
  static inline Real Tan (Real fValue, bool useTables = false) {
    return (!useTables) ? Real(tan(fValue)) : TanTable(fValue);
  }
  
  static inline Real DegreesToRadians(Real degrees) {
    return degrees * fDeg2Rad;
  }
  static inline Real RadiansToDegrees(Real radians) {
    return radians * fRad2Deg;
  }
  
  /** Compare 2 reals, using tolerance for inaccuracies.
  */
  static bool RealEqual(Real a, Real b,
                        Real tolerance = std::numeric_limits<Real>::epsilon());
                        
  /** Generates a value based on the Gaussian (normal) distribution function
  	with the given offset and scale parameters.
  */
  static Real gaussianDistribution(Real x, Real offset = 0.0f, Real scale = 1.0f);
  
  /** Clamp a value within an inclusive range. */
  template <typename T>
  static T Clamp(T val, T minval, T maxval)
  {
    assert (minval < maxval && "Invalid clamp range");
    return std::max(std::min(val, maxval), minval);
  }
  
  
  
  static const Real POS_INFINITY;
  static const Real NEG_INFINITY;
  static const Real PI;
  static const Real TWO_PI;
  static const Real HALF_PI;
  static const Real fDeg2Rad;
  static const Real fRad2Deg;
  
};


inline Radian operator * ( Real a, const Radian& b )
{
  return Radian ( a * b.valueRadians() );
}

inline Radian operator / ( Real a, const Radian& b )
{
  return Radian ( a / b.valueRadians() );
}

inline Degree operator * ( Real a, const Degree& b )
{
  return Degree ( a * b.valueDegrees() );
}

inline Degree operator / ( Real a, const Degree& b )
{
  return Degree ( a / b.valueDegrees() );
}

/** Standard 3-dimensional vector.
      @remarks
          A direction in 3D space represented as distances along the 3
          orthogonal axes (x, y, z). Note that positions, directions and
          scaling factors can be represented by a vector, depending on how
          you interpret the values.
  */
class Vector3
{
public:
  Real x, y, z;
  
public:
  inline Vector3()
  {
  }
  
  inline Vector3( const Real fX, const Real fY, const Real fZ )
      : x( fX ), y( fY ), z( fZ )
  {
  }
  
  inline explicit Vector3( const Real afCoordinate[3] )
      : x( afCoordinate[0] ),
      y( afCoordinate[1] ),
      z( afCoordinate[2] )
  {
  }
  
  inline explicit Vector3( const int afCoordinate[3] )
  {
    x = (Real)afCoordinate[0];
    y = (Real)afCoordinate[1];
    z = (Real)afCoordinate[2];
  }
  
  inline explicit Vector3( Real* const r )
      : x( r[0] ), y( r[1] ), z( r[2] )
  {
  }
  
  inline explicit Vector3( const Real scaler )
      : x( scaler )
      , y( scaler )
      , z( scaler )
  {
  }
  
  
  /** Exchange the contents of this vector with another.
  */
  inline void swap(Vector3& other)
  {
    std::swap(x, other.x);
    std::swap(y, other.y);
    std::swap(z, other.z);
  }
  
  inline Real operator [] ( const size_t i ) const
  {
    assert( i < 3 );
    
    return *(&x+i);
  }
  
  inline Real& operator [] ( const size_t i )
  {
    assert( i < 3 );
    
    return *(&x+i);
  }
  /// Pointer accessor for direct copying
  inline Real* ptr()
  {
    return &x;
  }
  /// Pointer accessor for direct copying
  inline const Real* ptr() const
  {
    return &x;
  }
  
  /** Assigns the value of the other vector.
      @param
          rkVector The other vector
  */
  inline Vector3& operator = ( const Vector3& rkVector )
  {
    x = rkVector.x;
    y = rkVector.y;
    z = rkVector.z;
    
    return *this;
  }
  
  inline Vector3& operator = ( const Real fScaler )
  {
    x = fScaler;
    y = fScaler;
    z = fScaler;
    
    return *this;
  }
  
  inline bool operator == ( const Vector3& rkVector ) const
  {
    return ( x == rkVector.x && y == rkVector.y && z == rkVector.z );
  }
  
  inline bool operator != ( const Vector3& rkVector ) const
  {
    return ( x != rkVector.x || y != rkVector.y || z != rkVector.z );
  }
  
  // arithmetic operations
  inline Vector3 operator + ( const Vector3& rkVector ) const
  {
    return Vector3(
             x + rkVector.x,
             y + rkVector.y,
             z + rkVector.z);
  }
  
  inline Vector3 operator - ( const Vector3& rkVector ) const
  {
    return Vector3(
             x - rkVector.x,
             y - rkVector.y,
             z - rkVector.z);
  }
  
  inline Vector3 operator * ( const Real fScalar ) const
  {
    return Vector3(
             x * fScalar,
             y * fScalar,
             z * fScalar);
  }
  
  inline Vector3 operator * ( const Vector3& rhs) const
  {
    return Vector3(
             x * rhs.x,
             y * rhs.y,
             z * rhs.z);
  }
  
  inline Vector3 operator / ( const Real fScalar ) const
  {
    assert( fScalar != 0.0 );
    
    Real fInv = 1.0f / fScalar;
    
    return Vector3(
             x * fInv,
             y * fInv,
             z * fInv);
  }
  
  inline Vector3 operator / ( const Vector3& rhs) const
  {
    return Vector3(
             x / rhs.x,
             y / rhs.y,
             z / rhs.z);
  }
  
  inline const Vector3& operator + () const
  {
    return *this;
  }
  
  inline Vector3 operator - () const
  {
    return Vector3(-x, -y, -z);
  }
  
  // overloaded operators to help Vector3
  inline friend Vector3 operator * ( const Real fScalar, const Vector3& rkVector )
  {
    return Vector3(
             fScalar * rkVector.x,
             fScalar * rkVector.y,
             fScalar * rkVector.z);
  }
  
  inline friend Vector3 operator / ( const Real fScalar, const Vector3& rkVector )
  {
    return Vector3(
             fScalar / rkVector.x,
             fScalar / rkVector.y,
             fScalar / rkVector.z);
  }
  
  inline friend Vector3 operator + (const Vector3& lhs, const Real rhs)
  {
    return Vector3(
             lhs.x + rhs,
             lhs.y + rhs,
             lhs.z + rhs);
  }
  
  inline friend Vector3 operator + (const Real lhs, const Vector3& rhs)
  {
    return Vector3(
             lhs + rhs.x,
             lhs + rhs.y,
             lhs + rhs.z);
  }
  
  inline friend Vector3 operator - (const Vector3& lhs, const Real rhs)
  {
    return Vector3(
             lhs.x - rhs,
             lhs.y - rhs,
             lhs.z - rhs);
  }
  
  inline friend Vector3 operator - (const Real lhs, const Vector3& rhs)
  {
    return Vector3(
             lhs - rhs.x,
             lhs - rhs.y,
             lhs - rhs.z);
  }
  
  // arithmetic updates
  inline Vector3& operator += ( const Vector3& rkVector )
  {
    x += rkVector.x;
    y += rkVector.y;
    z += rkVector.z;
    
    return *this;
  }
  
  inline Vector3& operator += ( const Real fScalar )
  {
    x += fScalar;
    y += fScalar;
    z += fScalar;
    return *this;
  }
  
  inline Vector3& operator -= ( const Vector3& rkVector )
  {
    x -= rkVector.x;
    y -= rkVector.y;
    z -= rkVector.z;
    
    return *this;
  }
  
  inline Vector3& operator -= ( const Real fScalar )
  {
    x -= fScalar;
    y -= fScalar;
    z -= fScalar;
    return *this;
  }
  
  inline Vector3& operator *= ( const Real fScalar )
  {
    x *= fScalar;
    y *= fScalar;
    z *= fScalar;
    return *this;
  }
  
  inline Vector3& operator *= ( const Vector3& rkVector )
  {
    x *= rkVector.x;
    y *= rkVector.y;
    z *= rkVector.z;
    
    return *this;
  }
  
  inline Vector3& operator /= ( const Real fScalar )
  {
    assert( fScalar != 0.0 );
    
    Real fInv = 1.0f / fScalar;
    
    x *= fInv;
    y *= fInv;
    z *= fInv;
    
    return *this;
  }
  
  inline Vector3& operator /= ( const Vector3& rkVector )
  {
    x /= rkVector.x;
    y /= rkVector.y;
    z /= rkVector.z;
    
    return *this;
  }
  
  
  /** Returns the length (magnitude) of the vector.
      @warning
          This operation requires a square root and is expensive in
          terms of CPU operations. If you don't need to know the exact
          length (e.g. for just comparing lengths) use squaredLength()
          instead.
  */
  inline Real length () const
  {
    return Math::Sqrt( x * x + y * y + z * z );
  }
  
  /** Returns the square of the length(magnitude) of the vector.
      @remarks
          This  method is for efficiency - calculating the actual
          length of a vector requires a square root, which is expensive
          in terms of the operations required. This method returns the
          square of the length of the vector, i.e. the same as the
          length but before the square root is taken. Use this if you
          want to find the longest / shortest vector without incurring
          the square root.
  */
  inline Real squaredLength () const
  {
    return x * x + y * y + z * z;
  }
  
  /** Returns the distance to another vector.
      @warning
          This operation requires a square root and is expensive in
          terms of CPU operations. If you don't need to know the exact
          distance (e.g. for just comparing distances) use squaredDistance()
          instead.
  */
  inline Real distance(const Vector3& rhs) const
  {
    return (*this - rhs).length();
  }
  
  /** Returns the square of the distance to another vector.
      @remarks
          This method is for efficiency - calculating the actual
          distance to another vector requires a square root, which is
          expensive in terms of the operations required. This method
          returns the square of the distance to another vector, i.e.
          the same as the distance but before the square root is taken.
          Use this if you want to find the longest / shortest distance
          without incurring the square root.
  */
  inline Real squaredDistance(const Vector3& rhs) const
  {
    return (*this - rhs).squaredLength();
  }
  
  /** Calculates the dot (scalar) product of this vector with another.
      @remarks
          The dot product can be used to calculate the angle between 2
          vectors. If both are unit vectors, the dot product is the
          cosine of the angle; otherwise the dot product must be
          divided by the product of the lengths of both vectors to get
          the cosine of the angle. This result can further be used to
          calculate the distance of a point from a plane.
      @param
          vec Vector with which to calculate the dot product (together
          with this one).
      @returns
          A float representing the dot product value.
  */
  inline Real dotProduct(const Vector3& vec) const
  {
    return x * vec.x + y * vec.y + z * vec.z;
  }
  
  /** Calculates the absolute dot (scalar) product of this vector with another.
      @remarks
          This function work similar dotProduct, except it use absolute value
          of each component of the vector to computing.
      @param
          vec Vector with which to calculate the absolute dot product (together
          with this one).
      @returns
          A Real representing the absolute dot product value.
  */
  inline Real absDotProduct(const Vector3& vec) const
  {
    return Math::Abs(x * vec.x) + Math::Abs(y * vec.y) + Math::Abs(z * vec.z);
  }
  
  /** Normalises the vector.
      @remarks
          This method normalises the vector such that it's
          length / magnitude is 1. The result is called a unit vector.
      @note
          This function will not crash for zero-sized vectors, but there
          will be no changes made to their components.
      @returns The previous length of the vector.
  */
  inline Real normalise()
  {
    Real fLength = Math::Sqrt( x * x + y * y + z * z );
    
    // Will also work for zero-sized vectors, but will change nothing
    if ( fLength > 1e-08 )
    {
      Real fInvLength = 1.0f / fLength;
      x *= fInvLength;
      y *= fInvLength;
      z *= fInvLength;
    }
    
    return fLength;
  }
  
  /** Calculates the cross-product of 2 vectors, i.e. the vector that
      lies perpendicular to them both.
      @remarks
          The cross-product is normally used to calculate the normal
          vector of a plane, by calculating the cross-product of 2
          non-equivalent vectors which lie on the plane (e.g. 2 edges
          of a triangle).
      @param
          vec Vector which, together with this one, will be used to
          calculate the cross-product.
      @returns
          A vector which is the result of the cross-product. This
          vector will <b>NOT</b> be normalised, to maximise efficiency
          - call Vector3::normalise on the result if you wish this to
          be done. As for which side the resultant vector will be on, the
          returned vector will be on the side from which the arc from 'this'
          to rkVector is anticlockwise, e.g. UNIT_Y.crossProduct(UNIT_Z)
          = UNIT_X, whilst UNIT_Z.crossProduct(UNIT_Y) = -UNIT_X.
  This is because OGRE uses a right-handed coordinate system.
      @par
          For a clearer explanation, look a the left and the bottom edges
          of your monitor's screen. Assume that the first vector is the
          left edge and the second vector is the bottom edge, both of
          them starting from the lower-left corner of the screen. The
          resulting vector is going to be perpendicular to both of them
          and will go <i>inside</i> the screen, towards the cathode tube
          (assuming you're using a CRT monitor, of course).
  */
  inline Vector3 crossProduct( const Vector3& rkVector ) const
  {
    return Vector3(
             y * rkVector.z - z * rkVector.y,
             z * rkVector.x - x * rkVector.z,
             x * rkVector.y - y * rkVector.x);
  }
  
  /** Returns a vector at a point half way between this and the passed
      in vector.
  */
  inline Vector3 midPoint( const Vector3& vec ) const
  {
    return Vector3(
             ( x + vec.x ) * 0.5f,
             ( y + vec.y ) * 0.5f,
             ( z + vec.z ) * 0.5f );
  }
  
  /** Returns true if the vector's scalar components are all greater
      that the ones of the vector it is compared against.
  */
  inline bool operator < ( const Vector3& rhs ) const
  {
    if ( x < rhs.x && y < rhs.y && z < rhs.z )
      return true;
    return false;
  }
  
  /** Returns true if the vector's scalar components are all smaller
      that the ones of the vector it is compared against.
  */
  inline bool operator > ( const Vector3& rhs ) const
  {
    if ( x > rhs.x && y > rhs.y && z > rhs.z )
      return true;
    return false;
  }
  
  /** Sets this vector's components to the minimum of its own and the
      ones of the passed in vector.
      @remarks
          'Minimum' in this case means the combination of the lowest
          value of x, y and z from both vectors. Lowest is taken just
          numerically, not magnitude, so -1 < 0.
  */
  inline void makeFloor( const Vector3& cmp )
  {
    if ( cmp.x < x ) x = cmp.x;
    if ( cmp.y < y ) y = cmp.y;
    if ( cmp.z < z ) z = cmp.z;
  }
  
  /** Sets this vector's components to the maximum of its own and the
      ones of the passed in vector.
      @remarks
          'Maximum' in this case means the combination of the highest
          value of x, y and z from both vectors. Highest is taken just
          numerically, not magnitude, so 1 > -3.
  */
  inline void makeCeil( const Vector3& cmp )
  {
    if ( cmp.x > x ) x = cmp.x;
    if ( cmp.y > y ) y = cmp.y;
    if ( cmp.z > z ) z = cmp.z;
  }
  
  /** Generates a vector perpendicular to this vector (eg an 'up' vector).
      @remarks
          This method will return a vector which is perpendicular to this
          vector. There are an infinite number of possibilities but this
          method will guarantee to generate one of them. If you need more
          control you should use the Quaternion class.
  */
  inline Vector3 perpendicular(void) const
  {
    static const Real fSquareZero = (Real)(1e-06 * 1e-06);
    
    Vector3 perp = this->crossProduct( Vector3::UNIT_X );
    
    // Check length
    if ( perp.squaredLength() < fSquareZero )
    {
      /* This vector is the Y axis multiplied by a scalar, so we have
         to use another axis.
      */
      perp = this->crossProduct( Vector3::UNIT_Y );
    }
    perp.normalise();
    
    return perp;
  }
  /** Generates a new random vector which deviates from this vector by a
      given angle in a random direction.
      @remarks
          This method assumes that the random number generator has already
          been seeded appropriately.
      @param
          angle The angle at which to deviate
      @param
          up Any vector perpendicular to this one (which could generated
          by cross-product of this vector and any other non-colinear
          vector). If you choose not to provide this the function will
          derive one on it's own, however if you provide one yourself the
          function will be faster (this allows you to reuse up vectors if
          you call this method more than once)
      @returns
          A random vector which deviates from this vector by angle. This
          vector will not be normalised, normalise it if you wish
          afterwards.
  */
  inline Vector3 randomDeviant(
    const Radian& angle,
    const Vector3& up = Vector3::ZERO ) const;
    
  /** Gets the angle between 2 vectors.
  @remarks
  	Vectors do not have to be unit-length but must represent directions.
  */
  inline Radian angleBetween(const Vector3& dest);
  
  /** Gets the shortest arc quaternion to rotate this vector to the destination
      vector.
  @remarks
      If you call this with a dest vector that is close to the inverse
      of this vector, we will rotate 180 degrees around the 'fallbackAxis'
  (if specified, or a generated axis if not) since in this case
  ANY axis of rotation is valid.
  */
  Quaternion getRotationTo(const Vector3& dest,
                           const Vector3& fallbackAxis = Vector3::ZERO) const;
                           
  /** Returns true if this vector is zero length. */
  inline bool isZeroLength(void) const
  {
    Real sqlen = (x * x) + (y * y) + (z * z);
    return (sqlen < (1e-06 * 1e-06));
    
  }
  
  /** As normalise, except that this vector is unaffected and the
      normalised vector is returned as a copy. */
  inline Vector3 normalisedCopy(void) const
  {
    Vector3 ret = *this;
    ret.normalise();
    return ret;
  }
  
  /** Calculates a reflection vector to the plane with the given normal .
  @remarks NB assumes 'this' is pointing AWAY FROM the plane, invert if it is not.
  */
  inline Vector3 reflect(const Vector3& normal) const
  {
    return Vector3( *this - ( 2 * this->dotProduct(normal) * normal ) );
  }
  
  /** Returns whether this vector is within a positional tolerance
  	of another vector.
  @param rhs The vector to compare with
  @param tolerance The amount that each element of the vector may vary by
  	and still be considered equal
  */
  inline bool positionEquals(const Vector3& rhs, Real tolerance = 1e-03) const
  {
    return Math::RealEqual(x, rhs.x, tolerance) &&
           Math::RealEqual(y, rhs.y, tolerance) &&
           Math::RealEqual(z, rhs.z, tolerance);
           
  }
  
  /** Returns whether this vector is within a positional tolerance
  	of another vector, also take scale of the vectors into account.
  @param rhs The vector to compare with
  @param tolerance The amount (related to the scale of vectors) that distance
          of the vector may vary by and still be considered close
  */
  inline bool positionCloses(const Vector3& rhs, Real tolerance = 1e-03f) const
  {
    return squaredDistance(rhs) <=
           (squaredLength() + rhs.squaredLength()) * tolerance;
  }
  
  /** Returns whether this vector is within a directional tolerance
  	of another vector.
  @param rhs The vector to compare with
  @param tolerance The maximum angle by which the vectors may vary and
  	still be considered equal
  @note Both vectors should be normalised.
  */
  inline bool directionEquals(const Vector3& rhs,
                              const Radian& tolerance) const
  {
    Real dot = dotProduct(rhs);
    Radian angle = Math::ACos(dot);
    
    return Math::Abs(angle.valueRadians()) <= tolerance.valueRadians();
    
  }
  
  /// Check whether this vector contains valid values
  inline bool isNaN() const
  {
    return Math::isNaN(x) || Math::isNaN(y) || Math::isNaN(z);
  }
  
  /// Extract the primary (dominant) axis from this direction vector
  inline Vector3 primaryAxis() const
  {
    Real absx = Math::Abs(x);
    Real absy = Math::Abs(y);
    Real absz = Math::Abs(z);
    if (absx > absy)
      if (absx > absz)
        return x > 0 ? Vector3::UNIT_X : Vector3::NEGATIVE_UNIT_X;
      else
        return z > 0 ? Vector3::UNIT_Z : Vector3::NEGATIVE_UNIT_Z;
    else // absx <= absy
      if (absy > absz)
        return y > 0 ? Vector3::UNIT_Y : Vector3::NEGATIVE_UNIT_Y;
      else
        return z > 0 ? Vector3::UNIT_Z : Vector3::NEGATIVE_UNIT_Z;
        
        
  }
  
  // special points
  static const Vector3 ZERO;
  static const Vector3 UNIT_X;
  static const Vector3 UNIT_Y;
  static const Vector3 UNIT_Z;
  static const Vector3 NEGATIVE_UNIT_X;
  static const Vector3 NEGATIVE_UNIT_Y;
  static const Vector3 NEGATIVE_UNIT_Z;
  static const Vector3 UNIT_SCALE;
  
  /** Function for writing to a stream.
  */
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Vector3& v )
  {
    o << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
    return o;
  }
  
  inline std::string to_s() const
  {
    std::ostringstream o;
    o << x << ", " << y << "," << z;
    return o.str();
  }
  
};

/** Implementation of a Quaternion, i.e. a rotation around an axis.
    */
class OGREMATHLIB_EXPORT Quaternion
{
public:
  /// Default constructor, initializes to identity rotation (aka 0�)
  inline Quaternion ()
      : w(1), x(0), y(0), z(0)
  {
  }
  /// Construct from an explicit list of values
  inline Quaternion (
    Real fW,
    Real fX, Real fY, Real fZ)
      : w(fW), x(fX), y(fY), z(fZ)
  {
  }
  /// Construct a quaternion from a rotation matrix
  inline Quaternion(const Matrix3& rot)
  {
    this->FromRotationMatrix(rot);
  }
  /// Construct a quaternion from an angle/axis
  inline Quaternion(const Radian& rfAngle, const Vector3& rkAxis)
  {
    this->FromAngleAxis(rfAngle, rkAxis);
  }
  /// Construct a quaternion from 3 orthonormal local axes
  inline Quaternion(const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis)
  {
    this->FromAxes(xaxis, yaxis, zaxis);
  }
  /// Construct a quaternion from 3 orthonormal local axes
  inline Quaternion(const Vector3* akAxis)
  {
    this->FromAxes(akAxis);
  }
  /// Construct a quaternion from 4 manual w/x/y/z values
  inline Quaternion(Real* valptr)
  {
    memcpy(&w, valptr, sizeof(Real)*4);
  }
  
  /** Exchange the contents of this quaternion with another.
  */
  inline void swap(Quaternion& other)
  {
    std::swap(w, other.w);
    std::swap(x, other.x);
    std::swap(y, other.y);
    std::swap(z, other.z);
  }
  
  /// Array accessor operator
  inline Real operator [] ( const size_t i ) const
  {
    assert( i < 4 );
    
    return *(&w+i);
  }
  
  /// Array accessor operator
  inline Real& operator [] ( const size_t i )
  {
    assert( i < 4 );
    
    return *(&w+i);
  }
  
  /// Pointer accessor for direct copying
  inline Real* ptr()
  {
    return &w;
  }
  
  /// Pointer accessor for direct copying
  inline const Real* ptr() const
  {
    return &w;
  }
  
  void FromRotationMatrix (const Matrix3& kRot);
  void ToRotationMatrix (Matrix3& kRot) const;
  void FromAngleAxis (const Radian& rfAngle, const Vector3& rkAxis);
  void ToAngleAxis (Radian& rfAngle, Vector3& rkAxis) const;
  inline void ToAngleAxis (Degree& dAngle, Vector3& rkAxis) const {
    Radian rAngle;
    ToAngleAxis ( rAngle, rkAxis );
    dAngle = rAngle;
  }
  void FromAxes (const Vector3* akAxis);
  void FromAxes (const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis);
  void ToAxes (Vector3* akAxis) const;
  void ToAxes (Vector3& xAxis, Vector3& yAxis, Vector3& zAxis) const;
  /// Get the local x-axis
  Vector3 xAxis(void) const;
  /// Get the local y-axis
  Vector3 yAxis(void) const;
  /// Get the local z-axis
  Vector3 zAxis(void) const;
  
  inline Quaternion& operator= (const Quaternion& rkQ)
  {
    w = rkQ.w;
    x = rkQ.x;
    y = rkQ.y;
    z = rkQ.z;
    return *this;
  }
  Quaternion operator+ (const Quaternion& rkQ) const;
  Quaternion operator- (const Quaternion& rkQ) const;
  Quaternion operator* (const Quaternion& rkQ) const;
  Quaternion operator* (Real fScalar) const;
  OGREMATHLIB_EXPORT friend Quaternion operator* (Real fScalar,
      const Quaternion& rkQ);
  Quaternion operator- () const;
  inline bool operator== (const Quaternion& rhs) const
  {
    return (rhs.x == x) && (rhs.y == y) &&
           (rhs.z == z) && (rhs.w == w);
  }
  inline bool operator!= (const Quaternion& rhs) const
  {
    return !operator==(rhs);
  }
  // functions of a quaternion
  Real Dot (const Quaternion& rkQ) const;  // dot product
  Real Norm () const;  // squared-length
  /// Normalises this quaternion, and returns the previous length
  Real normalise(void);
  Quaternion Inverse () const;  // apply to non-zero quaternion
  Quaternion UnitInverse () const;  // apply to unit-length quaternion
  Quaternion Exp () const;
  Quaternion Log () const;
  
  // rotation of a vector by a quaternion
  Vector3 operator* (const Vector3& rkVector) const;
  
  /** Calculate the local roll element of this quaternion.
  @param reprojectAxis By default the method returns the 'intuitive' result
  that is, if you projected the local Y of the quaternion onto the X and
  Y axes, the angle between them is returned. If set to false though, the
  result is the actual yaw that will be used to implement the quaternion,
  which is the shortest possible path to get to the same orientation and
  may involve less axial rotation.
  */
  Radian getRoll(bool reprojectAxis = true) const;
  /** Calculate the local pitch element of this quaternion
  @param reprojectAxis By default the method returns the 'intuitive' result
  that is, if you projected the local Z of the quaternion onto the X and
  Y axes, the angle between them is returned. If set to true though, the
  result is the actual yaw that will be used to implement the quaternion,
  which is the shortest possible path to get to the same orientation and
  may involve less axial rotation.
  */
  Radian getPitch(bool reprojectAxis = true) const;
  /** Calculate the local yaw element of this quaternion
  @param reprojectAxis By default the method returns the 'intuitive' result
  that is, if you projected the local Z of the quaternion onto the X and
  Z axes, the angle between them is returned. If set to true though, the
  result is the actual yaw that will be used to implement the quaternion,
  which is the shortest possible path to get to the same orientation and
  may involve less axial rotation.
  */
  Radian getYaw(bool reprojectAxis = true) const;
  /// Equality with tolerance (tolerance is max angle difference)
  bool equals(const Quaternion& rhs, const Radian& tolerance) const;
  
  // spherical linear interpolation
  static Quaternion Slerp (Real fT, const Quaternion& rkP,
                           const Quaternion& rkQ, bool shortestPath = false);
                           
  static Quaternion SlerpExtraSpins (Real fT,
                                     const Quaternion& rkP, const Quaternion& rkQ,
                                     int iExtraSpins);
                                     
  // setup for spherical quadratic interpolation
  static void Intermediate (const Quaternion& rkQ0,
                            const Quaternion& rkQ1, const Quaternion& rkQ2,
                            Quaternion& rka, Quaternion& rkB);
                            
  // spherical quadratic interpolation
  static Quaternion Squad (Real fT, const Quaternion& rkP,
                           const Quaternion& rkA, const Quaternion& rkB,
                           const Quaternion& rkQ, bool shortestPath = false);
                           
  // normalised linear interpolation - faster but less accurate (non-constant rotation velocity)
  static Quaternion nlerp(Real fT, const Quaternion& rkP,
                          const Quaternion& rkQ, bool shortestPath = false);
                          
  // cutoff for sine near zero
  static const Real ms_fEpsilon;
  
  // special values
  static const Quaternion ZERO;
  static const Quaternion IDENTITY;
  
  Real w, x, y, z;
  
  /// Check whether this quaternion contains valid values
  inline bool isNaN() const
  {
    return Math::isNaN(x) || Math::isNaN(y) || Math::isNaN(z) || Math::isNaN(w);
  }
  
  /** Function for writing to a stream. Outputs "Quaternion(w, x, y, z)" with w,x,y,z
      being the member values of the quaternion.
  */
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Quaternion& q )
  {
    o << "Quaternion(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
    return o;
  }
  
  inline std::string to_s() const
  {
    std::ostringstream o;
    o << w << "," << x << ", " << y << "," << z;
    return o.str();
  }
  
};


// NB All code adapted from Wild Magic 0.2 Matrix math (free source code)
// http://www.geometrictools.com/

// NOTE.  The (x,y,z) coordinate system is assumed to be right-handed.
// Coordinate axis rotation matrices are of the form
//   RX =    1       0       0
//           0     cos(t) -sin(t)
//           0     sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the yz-plane
//   RY =  cos(t)    0     sin(t)
//           0       1       0
//        -sin(t)    0     cos(t)
// where t > 0 indicates a counterclockwise rotation in the zx-plane
//   RZ =  cos(t) -sin(t)    0
//         sin(t)  cos(t)    0
//           0       0       1
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

/** A 3x3 matrix which can represent rotations around axes.
      @note
          <b>All the code is adapted from the Wild Magic 0.2 Matrix
          library (http://www.geometrictools.com/).</b>
      @par
          The coordinate system is assumed to be <b>right-handed</b>.
  */
class OGREMATHLIB_EXPORT Matrix3
{
public:
  /** Default constructor.
      @note
          It does <b>NOT</b> initialize the matrix for efficiency.
  */
  inline Matrix3 () {}
  inline explicit Matrix3 (const Real arr[3][3])
  {
    memcpy(m,arr,9*sizeof(Real));
  }
  inline Matrix3 (const Matrix3& rkMatrix)
  {
    memcpy(m,rkMatrix.m,9*sizeof(Real));
  }
  Matrix3 (Real fEntry00, Real fEntry01, Real fEntry02,
           Real fEntry10, Real fEntry11, Real fEntry12,
           Real fEntry20, Real fEntry21, Real fEntry22)
  {
    m[0][0] = fEntry00;
    m[0][1] = fEntry01;
    m[0][2] = fEntry02;
    m[1][0] = fEntry10;
    m[1][1] = fEntry11;
    m[1][2] = fEntry12;
    m[2][0] = fEntry20;
    m[2][1] = fEntry21;
    m[2][2] = fEntry22;
  }
  
  /** Exchange the contents of this matrix with another.
  */
  inline void swap(Matrix3& other)
  {
    std::swap(m[0][0], other.m[0][0]);
    std::swap(m[0][1], other.m[0][1]);
    std::swap(m[0][2], other.m[0][2]);
    std::swap(m[1][0], other.m[1][0]);
    std::swap(m[1][1], other.m[1][1]);
    std::swap(m[1][2], other.m[1][2]);
    std::swap(m[2][0], other.m[2][0]);
    std::swap(m[2][1], other.m[2][1]);
    std::swap(m[2][2], other.m[2][2]);
  }
  
  // member access, allows use of construct mat[r][c]
  inline Real* operator[] (size_t iRow) const
  {
    return (Real*)m[iRow];
  }
  /*inline operator Real* ()
  {
  return (Real*)m[0];
  }*/
  Vector3 GetColumn (size_t iCol) const;
  void SetColumn(size_t iCol, const Vector3& vec);
  void FromAxes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis);
  
  // assignment and comparison
  inline Matrix3& operator= (const Matrix3& rkMatrix)
  {
    memcpy(m,rkMatrix.m,9*sizeof(Real));
    return *this;
  }
  bool operator== (const Matrix3& rkMatrix) const;
  inline bool operator!= (const Matrix3& rkMatrix) const
  {
    return !operator==(rkMatrix);
  }
  
  // arithmetic operations
  Matrix3 operator+ (const Matrix3& rkMatrix) const;
  Matrix3 operator- (const Matrix3& rkMatrix) const;
  Matrix3 operator* (const Matrix3& rkMatrix) const;
  Matrix3 operator- () const;
  
  // matrix * vector [3x3 * 3x1 = 3x1]
  Vector3 operator* (const Vector3& rkVector) const;
  
  // vector * matrix [1x3 * 3x3 = 1x3]
  OGREMATHLIB_EXPORT friend Vector3 operator* (const Vector3& rkVector,
      const Matrix3& rkMatrix);
      
  // matrix * scalar
  Matrix3 operator* (Real fScalar) const;
  
  // scalar * matrix
  OGREMATHLIB_EXPORT friend Matrix3 operator* (Real fScalar, const Matrix3& rkMatrix);
  
  // utilities
  Matrix3 Transpose () const;
  bool Inverse (Matrix3& rkInverse, Real fTolerance = 1e-06) const;
  Matrix3 Inverse (Real fTolerance = 1e-06) const;
  Real Determinant () const;
  
  // singular value decomposition
  void SingularValueDecomposition (Matrix3& rkL, Vector3& rkS,
                                   Matrix3& rkR) const;
  void SingularValueComposition (const Matrix3& rkL,
                                 const Vector3& rkS, const Matrix3& rkR);
                                 
  // Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
  void Orthonormalize ();
  
  // orthogonal Q, diagonal D, upper triangular U stored as (u01,u02,u12)
  void QDUDecomposition (Matrix3& rkQ, Vector3& rkD,
                         Vector3& rkU) const;
                         
  Real SpectralNorm () const;
  
  // matrix must be orthonormal
  void ToAxisAngle (Vector3& rkAxis, Radian& rfAngle) const;
  inline void ToAxisAngle (Vector3& rkAxis, Degree& rfAngle) const {
    Radian r;
    ToAxisAngle ( rkAxis, r );
    rfAngle = r;
  }
  void FromAxisAngle (const Vector3& rkAxis, const Radian& fRadians);
  
  // The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
  // where yaw is rotation about the Up vector, pitch is rotation about the
  // Right axis, and roll is rotation about the Direction axis.
  bool ToEulerAnglesXYZ (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  bool ToEulerAnglesXZY (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  bool ToEulerAnglesYXZ (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  bool ToEulerAnglesYZX (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  bool ToEulerAnglesZXY (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  bool ToEulerAnglesZYX (Radian& rfYAngle, Radian& rfPAngle,
                         Radian& rfRAngle) const;
  void FromEulerAnglesXYZ (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  void FromEulerAnglesXZY (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  void FromEulerAnglesYXZ (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  void FromEulerAnglesYZX (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  void FromEulerAnglesZXY (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  void FromEulerAnglesZYX (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle);
  // eigensolver, matrix must be symmetric
  void EigenSolveSymmetric (Real afEigenvalue[3],
                            Vector3 akEigenvector[3]) const;
                            
  static void TensorProduct (const Vector3& rkU, const Vector3& rkV,
                             Matrix3& rkProduct);
                             
  /** Determines if this matrix involves a scaling. */
  inline bool hasScale() const
  {
    // check magnitude of column vectors (==local axes)
    Real t = m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
    t = m[0][1] * m[0][1] + m[1][1] * m[1][1] + m[2][1] * m[2][1];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
    t = m[0][2] * m[0][2] + m[1][2] * m[1][2] + m[2][2] * m[2][2];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
      
    return false;
  }
  
  /** Function for writing to a stream.
  */
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Matrix3& mat )
  {
    o << "Matrix3(" << mat[0][0] << ", " << mat[0][1] << ", " << mat[0][2] << ", "
    << mat[1][0] << ", " << mat[1][1] << ", " << mat[1][2] << ", "
    << mat[2][0] << ", " << mat[2][1] << ", " << mat[2][2] << ")";
    return o;
  }
  
  static const Real EPSILON;
  static const Matrix3 ZERO;
  static const Matrix3 IDENTITY;
  
protected:
  // support for eigensolver
  void Tridiagonal (Real afDiag[3], Real afSubDiag[3]);
  bool QLAlgorithm (Real afDiag[3], Real afSubDiag[3]);
  
  // support for singular value decomposition
  static const Real ms_fSvdEpsilon;
  static const unsigned int ms_iSvdMaxIterations;
  static void Bidiagonalize (Matrix3& kA, Matrix3& kL,
                             Matrix3& kR);
  static void GolubKahanStep (Matrix3& kA, Matrix3& kL,
                              Matrix3& kR);
                              
  // support for spectral norm
  static Real MaxCubicRoot (Real afCoeff[3]);
  
  Real m[3][3];
  
  // for faster access
  friend class Matrix4;
};

/** Class encapsulating a standard 4x4 homogeneous matrix.
      @remarks
          OGRE uses column vectors when applying matrix multiplications,
          This means a vector is represented as a single column, 4-row
          matrix. This has the effect that the transformations implemented
          by the matrices happens right-to-left e.g. if vector V is to be
          transformed by M1 then M2 then M3, the calculation would be
          M3 * M2 * M1 * V. The order that matrices are concatenated is
          vital since matrix multiplication is not commutative, i.e. you
          can get a different result if you concatenate in the wrong order.
      @par
          The use of column vectors and right-to-left ordering is the
          standard in most mathematical texts, and is the same as used in
          OpenGL. It is, however, the opposite of Direct3D, which has
          inexplicably chosen to differ from the accepted standard and uses
          row vectors and left-to-right matrix multiplication.
      @par
          OGRE deals with the differences between D3D and OpenGL etc.
          internally when operating through different render systems. OGRE
          users only need to conform to standard maths conventions, i.e.
          right-to-left matrix multiplication, (OGRE transposes matrices it
          passes to D3D to compensate).
      @par
          The generic form M * V which shows the layout of the matrix
          entries is shown below:
          <pre>
              [ m[0][0]  m[0][1]  m[0][2]  m[0][3] ]   {x}
              | m[1][0]  m[1][1]  m[1][2]  m[1][3] | * {y}
              | m[2][0]  m[2][1]  m[2][2]  m[2][3] |   {z}
              [ m[3][0]  m[3][1]  m[3][2]  m[3][3] ]   {1}
          </pre>
  */
class OGREMATHLIB_EXPORT Matrix4
{
protected:
  /// The matrix entries, indexed by [row][col].
  union {
    Real m[4][4];
    Real _m[16];
  };
public:
  /** Default constructor.
      @note
          It does <b>NOT</b> initialize the matrix for efficiency.
  */
  inline Matrix4()
  {
  }
  
  inline Matrix4(
    Real m00, Real m01, Real m02, Real m03,
    Real m10, Real m11, Real m12, Real m13,
    Real m20, Real m21, Real m22, Real m23,
    Real m30, Real m31, Real m32, Real m33 )
  {
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;
    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;
    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;
    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
  }
  
  /** Creates a standard 4x4 transformation matrix with a zero translation part from a rotation/scaling 3x3 matrix.
   */
  
  inline Matrix4(const Matrix3& m3x3)
  {
    operator=(IDENTITY);
    operator=(m3x3);
  }
  
  /** Creates a standard 4x4 transformation matrix with a zero translation part from a rotation/scaling Quaternion.
   */
  
  inline Matrix4(const Quaternion& rot)
  {
    Matrix3 m3x3;
    rot.ToRotationMatrix(m3x3);
    operator=(IDENTITY);
    operator=(m3x3);
  }
  
  
  /** Exchange the contents of this matrix with another.
  */
  inline void swap(Matrix4& other)
  {
    std::swap(m[0][0], other.m[0][0]);
    std::swap(m[0][1], other.m[0][1]);
    std::swap(m[0][2], other.m[0][2]);
    std::swap(m[0][3], other.m[0][3]);
    std::swap(m[1][0], other.m[1][0]);
    std::swap(m[1][1], other.m[1][1]);
    std::swap(m[1][2], other.m[1][2]);
    std::swap(m[1][3], other.m[1][3]);
    std::swap(m[2][0], other.m[2][0]);
    std::swap(m[2][1], other.m[2][1]);
    std::swap(m[2][2], other.m[2][2]);
    std::swap(m[2][3], other.m[2][3]);
    std::swap(m[3][0], other.m[3][0]);
    std::swap(m[3][1], other.m[3][1]);
    std::swap(m[3][2], other.m[3][2]);
    std::swap(m[3][3], other.m[3][3]);
  }
  
  inline Real* operator [] ( size_t iRow )
  {
    assert( iRow < 4 );
    return m[iRow];
  }
  
  inline const Real *operator [] ( size_t iRow ) const
  {
    assert( iRow < 4 );
    return m[iRow];
  }
  
  inline Matrix4 concatenate(const Matrix4 &m2) const
  {
    Matrix4 r;
    r.m[0][0] = m[0][0] * m2.m[0][0] + m[0][1] * m2.m[1][0] + m[0][2] * m2.m[2][0] + m[0][3] * m2.m[3][0];
    r.m[0][1] = m[0][0] * m2.m[0][1] + m[0][1] * m2.m[1][1] + m[0][2] * m2.m[2][1] + m[0][3] * m2.m[3][1];
    r.m[0][2] = m[0][0] * m2.m[0][2] + m[0][1] * m2.m[1][2] + m[0][2] * m2.m[2][2] + m[0][3] * m2.m[3][2];
    r.m[0][3] = m[0][0] * m2.m[0][3] + m[0][1] * m2.m[1][3] + m[0][2] * m2.m[2][3] + m[0][3] * m2.m[3][3];
    
    r.m[1][0] = m[1][0] * m2.m[0][0] + m[1][1] * m2.m[1][0] + m[1][2] * m2.m[2][0] + m[1][3] * m2.m[3][0];
    r.m[1][1] = m[1][0] * m2.m[0][1] + m[1][1] * m2.m[1][1] + m[1][2] * m2.m[2][1] + m[1][3] * m2.m[3][1];
    r.m[1][2] = m[1][0] * m2.m[0][2] + m[1][1] * m2.m[1][2] + m[1][2] * m2.m[2][2] + m[1][3] * m2.m[3][2];
    r.m[1][3] = m[1][0] * m2.m[0][3] + m[1][1] * m2.m[1][3] + m[1][2] * m2.m[2][3] + m[1][3] * m2.m[3][3];
    
    r.m[2][0] = m[2][0] * m2.m[0][0] + m[2][1] * m2.m[1][0] + m[2][2] * m2.m[2][0] + m[2][3] * m2.m[3][0];
    r.m[2][1] = m[2][0] * m2.m[0][1] + m[2][1] * m2.m[1][1] + m[2][2] * m2.m[2][1] + m[2][3] * m2.m[3][1];
    r.m[2][2] = m[2][0] * m2.m[0][2] + m[2][1] * m2.m[1][2] + m[2][2] * m2.m[2][2] + m[2][3] * m2.m[3][2];
    r.m[2][3] = m[2][0] * m2.m[0][3] + m[2][1] * m2.m[1][3] + m[2][2] * m2.m[2][3] + m[2][3] * m2.m[3][3];
    
    r.m[3][0] = m[3][0] * m2.m[0][0] + m[3][1] * m2.m[1][0] + m[3][2] * m2.m[2][0] + m[3][3] * m2.m[3][0];
    r.m[3][1] = m[3][0] * m2.m[0][1] + m[3][1] * m2.m[1][1] + m[3][2] * m2.m[2][1] + m[3][3] * m2.m[3][1];
    r.m[3][2] = m[3][0] * m2.m[0][2] + m[3][1] * m2.m[1][2] + m[3][2] * m2.m[2][2] + m[3][3] * m2.m[3][2];
    r.m[3][3] = m[3][0] * m2.m[0][3] + m[3][1] * m2.m[1][3] + m[3][2] * m2.m[2][3] + m[3][3] * m2.m[3][3];
    
    return r;
  }
  
  /** Matrix concatenation using '*'.
  */
  inline Matrix4 operator * ( const Matrix4 &m2 ) const
  {
    return concatenate( m2 );
  }
  
  /** Vector transformation using '*'.
      @remarks
          Transforms the given 3-D vector by the matrix, projecting the
          result back into <i>w</i> = 1.
      @note
          This means that the initial <i>w</i> is considered to be 1.0,
          and then all the tree elements of the resulting 3-D vector are
          divided by the resulting <i>w</i>.
  */
  inline Vector3 operator * ( const Vector3 &v ) const
  {
    Vector3 r;
    
    Real fInvW = 1.0f / ( m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] );
    
    r.x = ( m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] ) * fInvW;
    r.y = ( m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] ) * fInvW;
    r.z = ( m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] ) * fInvW;
    
    return r;
  }
  
  
  
  /** Matrix addition.
  */
  inline Matrix4 operator + ( const Matrix4 &m2 ) const
  {
    Matrix4 r;
    
    r.m[0][0] = m[0][0] + m2.m[0][0];
    r.m[0][1] = m[0][1] + m2.m[0][1];
    r.m[0][2] = m[0][2] + m2.m[0][2];
    r.m[0][3] = m[0][3] + m2.m[0][3];
    
    r.m[1][0] = m[1][0] + m2.m[1][0];
    r.m[1][1] = m[1][1] + m2.m[1][1];
    r.m[1][2] = m[1][2] + m2.m[1][2];
    r.m[1][3] = m[1][3] + m2.m[1][3];
    
    r.m[2][0] = m[2][0] + m2.m[2][0];
    r.m[2][1] = m[2][1] + m2.m[2][1];
    r.m[2][2] = m[2][2] + m2.m[2][2];
    r.m[2][3] = m[2][3] + m2.m[2][3];
    
    r.m[3][0] = m[3][0] + m2.m[3][0];
    r.m[3][1] = m[3][1] + m2.m[3][1];
    r.m[3][2] = m[3][2] + m2.m[3][2];
    r.m[3][3] = m[3][3] + m2.m[3][3];
    
    return r;
  }
  
  /** Matrix subtraction.
  */
  inline Matrix4 operator - ( const Matrix4 &m2 ) const
  {
    Matrix4 r;
    r.m[0][0] = m[0][0] - m2.m[0][0];
    r.m[0][1] = m[0][1] - m2.m[0][1];
    r.m[0][2] = m[0][2] - m2.m[0][2];
    r.m[0][3] = m[0][3] - m2.m[0][3];
    
    r.m[1][0] = m[1][0] - m2.m[1][0];
    r.m[1][1] = m[1][1] - m2.m[1][1];
    r.m[1][2] = m[1][2] - m2.m[1][2];
    r.m[1][3] = m[1][3] - m2.m[1][3];
    
    r.m[2][0] = m[2][0] - m2.m[2][0];
    r.m[2][1] = m[2][1] - m2.m[2][1];
    r.m[2][2] = m[2][2] - m2.m[2][2];
    r.m[2][3] = m[2][3] - m2.m[2][3];
    
    r.m[3][0] = m[3][0] - m2.m[3][0];
    r.m[3][1] = m[3][1] - m2.m[3][1];
    r.m[3][2] = m[3][2] - m2.m[3][2];
    r.m[3][3] = m[3][3] - m2.m[3][3];
    
    return r;
  }
  
  /** Tests 2 matrices for equality.
  */
  inline bool operator == ( const Matrix4& m2 ) const
  {
    if (
      m[0][0] != m2.m[0][0] || m[0][1] != m2.m[0][1] || m[0][2] != m2.m[0][2] || m[0][3] != m2.m[0][3] ||
      m[1][0] != m2.m[1][0] || m[1][1] != m2.m[1][1] || m[1][2] != m2.m[1][2] || m[1][3] != m2.m[1][3] ||
      m[2][0] != m2.m[2][0] || m[2][1] != m2.m[2][1] || m[2][2] != m2.m[2][2] || m[2][3] != m2.m[2][3] ||
      m[3][0] != m2.m[3][0] || m[3][1] != m2.m[3][1] || m[3][2] != m2.m[3][2] || m[3][3] != m2.m[3][3] )
      return false;
    return true;
  }
  
  /** Tests 2 matrices for inequality.
  */
  inline bool operator != ( const Matrix4& m2 ) const
  {
    if (
      m[0][0] != m2.m[0][0] || m[0][1] != m2.m[0][1] || m[0][2] != m2.m[0][2] || m[0][3] != m2.m[0][3] ||
      m[1][0] != m2.m[1][0] || m[1][1] != m2.m[1][1] || m[1][2] != m2.m[1][2] || m[1][3] != m2.m[1][3] ||
      m[2][0] != m2.m[2][0] || m[2][1] != m2.m[2][1] || m[2][2] != m2.m[2][2] || m[2][3] != m2.m[2][3] ||
      m[3][0] != m2.m[3][0] || m[3][1] != m2.m[3][1] || m[3][2] != m2.m[3][2] || m[3][3] != m2.m[3][3] )
      return true;
    return false;
  }
  
  /** Assignment from 3x3 matrix.
  */
  inline void operator = ( const Matrix3& mat3 )
  {
    m[0][0] = mat3.m[0][0];
    m[0][1] = mat3.m[0][1];
    m[0][2] = mat3.m[0][2];
    m[1][0] = mat3.m[1][0];
    m[1][1] = mat3.m[1][1];
    m[1][2] = mat3.m[1][2];
    m[2][0] = mat3.m[2][0];
    m[2][1] = mat3.m[2][1];
    m[2][2] = mat3.m[2][2];
  }
  
  inline Matrix4 transpose(void) const
  {
    return Matrix4(m[0][0], m[1][0], m[2][0], m[3][0],
                   m[0][1], m[1][1], m[2][1], m[3][1],
                   m[0][2], m[1][2], m[2][2], m[3][2],
                   m[0][3], m[1][3], m[2][3], m[3][3]);
  }
  
  /*
  -----------------------------------------------------------------------
  Translation Transformation
  -----------------------------------------------------------------------
  */
  /** Sets the translation transformation part of the matrix.
  */
  inline void setTrans( const Vector3& v )
  {
    m[0][3] = v.x;
    m[1][3] = v.y;
    m[2][3] = v.z;
  }
  
  /** Extracts the translation transformation part of the matrix.
   */
  inline Vector3 getTrans() const
  {
    return Vector3(m[0][3], m[1][3], m[2][3]);
  }
  
  
  /** Builds a translation matrix
  */
  inline void makeTrans( const Vector3& v )
  {
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[0][3] = v.x;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[1][3] = v.y;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
    m[2][3] = v.z;
    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;
  }
  
  inline void makeTrans( Real tx, Real ty, Real tz )
  {
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[0][3] = tx;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[1][3] = ty;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
    m[2][3] = tz;
    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;
  }
  
  /** Gets a translation matrix.
  */
  inline static Matrix4 getTrans( const Vector3& v )
  {
    Matrix4 r;
    
    r.m[0][0] = 1.0;
    r.m[0][1] = 0.0;
    r.m[0][2] = 0.0;
    r.m[0][3] = v.x;
    r.m[1][0] = 0.0;
    r.m[1][1] = 1.0;
    r.m[1][2] = 0.0;
    r.m[1][3] = v.y;
    r.m[2][0] = 0.0;
    r.m[2][1] = 0.0;
    r.m[2][2] = 1.0;
    r.m[2][3] = v.z;
    r.m[3][0] = 0.0;
    r.m[3][1] = 0.0;
    r.m[3][2] = 0.0;
    r.m[3][3] = 1.0;
    
    return r;
  }
  
  /** Gets a translation matrix - variation for not using a vector.
  */
  inline static Matrix4 getTrans( Real t_x, Real t_y, Real t_z )
  {
    Matrix4 r;
    
    r.m[0][0] = 1.0;
    r.m[0][1] = 0.0;
    r.m[0][2] = 0.0;
    r.m[0][3] = t_x;
    r.m[1][0] = 0.0;
    r.m[1][1] = 1.0;
    r.m[1][2] = 0.0;
    r.m[1][3] = t_y;
    r.m[2][0] = 0.0;
    r.m[2][1] = 0.0;
    r.m[2][2] = 1.0;
    r.m[2][3] = t_z;
    r.m[3][0] = 0.0;
    r.m[3][1] = 0.0;
    r.m[3][2] = 0.0;
    r.m[3][3] = 1.0;
    
    return r;
  }
  
  /*
  -----------------------------------------------------------------------
  Scale Transformation
  -----------------------------------------------------------------------
  */
  /** Sets the scale part of the matrix.
  */
  inline void setScale( const Vector3& v )
  {
    m[0][0] = v.x;
    m[1][1] = v.y;
    m[2][2] = v.z;
  }
  
  /** Gets a scale matrix.
  */
  inline static Matrix4 getScale( const Vector3& v )
  {
    Matrix4 r;
    r.m[0][0] = v.x;
    r.m[0][1] = 0.0;
    r.m[0][2] = 0.0;
    r.m[0][3] = 0.0;
    r.m[1][0] = 0.0;
    r.m[1][1] = v.y;
    r.m[1][2] = 0.0;
    r.m[1][3] = 0.0;
    r.m[2][0] = 0.0;
    r.m[2][1] = 0.0;
    r.m[2][2] = v.z;
    r.m[2][3] = 0.0;
    r.m[3][0] = 0.0;
    r.m[3][1] = 0.0;
    r.m[3][2] = 0.0;
    r.m[3][3] = 1.0;
    
    return r;
  }
  
  /** Gets a scale matrix - variation for not using a vector.
  */
  inline static Matrix4 getScale( Real s_x, Real s_y, Real s_z )
  {
    Matrix4 r;
    r.m[0][0] = s_x;
    r.m[0][1] = 0.0;
    r.m[0][2] = 0.0;
    r.m[0][3] = 0.0;
    r.m[1][0] = 0.0;
    r.m[1][1] = s_y;
    r.m[1][2] = 0.0;
    r.m[1][3] = 0.0;
    r.m[2][0] = 0.0;
    r.m[2][1] = 0.0;
    r.m[2][2] = s_z;
    r.m[2][3] = 0.0;
    r.m[3][0] = 0.0;
    r.m[3][1] = 0.0;
    r.m[3][2] = 0.0;
    r.m[3][3] = 1.0;
    
    return r;
  }
  
  /** Extracts the rotation / scaling part of the Matrix as a 3x3 matrix.
  @param m3x3 Destination Matrix3
  */
  inline void extract3x3Matrix(Matrix3& m3x3) const
  {
    m3x3.m[0][0] = m[0][0];
    m3x3.m[0][1] = m[0][1];
    m3x3.m[0][2] = m[0][2];
    m3x3.m[1][0] = m[1][0];
    m3x3.m[1][1] = m[1][1];
    m3x3.m[1][2] = m[1][2];
    m3x3.m[2][0] = m[2][0];
    m3x3.m[2][1] = m[2][1];
    m3x3.m[2][2] = m[2][2];
    
  }
  
  /** Determines if this matrix involves a scaling. */
  inline bool hasScale() const
  {
    // check magnitude of column vectors (==local axes)
    Real t = m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
    t = m[0][1] * m[0][1] + m[1][1] * m[1][1] + m[2][1] * m[2][1];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
    t = m[0][2] * m[0][2] + m[1][2] * m[1][2] + m[2][2] * m[2][2];
    if (!Math::RealEqual(t, 1.0, (Real)1e-04))
      return true;
      
    return false;
  }
  
  /** Determines if this matrix involves a negative scaling. */
  inline bool hasNegativeScale() const
  {
    return determinant() < 0;
  }
  
  /** Extracts the rotation / scaling part as a quaternion from the Matrix.
       */
  inline Quaternion extractQuaternion() const
  {
    Matrix3 m3x3;
    extract3x3Matrix(m3x3);
    return Quaternion(m3x3);
  }
  
  static const Matrix4 ZERO;
  static const Matrix4 IDENTITY;
  /** Useful little matrix which takes 2D clipspace {-1, 1} to {0,1}
      and inverts the Y. */
  static const Matrix4 CLIPSPACE2DTOIMAGESPACE;
  
  inline Matrix4 operator*(Real scalar) const
  {
    return Matrix4(
             scalar*m[0][0], scalar*m[0][1], scalar*m[0][2], scalar*m[0][3],
             scalar*m[1][0], scalar*m[1][1], scalar*m[1][2], scalar*m[1][3],
             scalar*m[2][0], scalar*m[2][1], scalar*m[2][2], scalar*m[2][3],
             scalar*m[3][0], scalar*m[3][1], scalar*m[3][2], scalar*m[3][3]);
  }
  
  /** Function for writing to a stream.
  */
  inline OGREMATHLIB_EXPORT friend std::ostream& operator <<
  ( std::ostream& o, const Matrix4& mat )
  {
    o << "Matrix4(";
    for (size_t i = 0; i < 4; ++i)
    {
      o << " row" << (unsigned)i << "{";
      for (size_t j = 0; j < 4; ++j)
      {
        o << mat[i][j] << " ";
      }
      o << "}";
    }
    o << ")";
    return o;
  }
  
  Matrix4 adjoint() const;
  Real determinant() const;
  Matrix4 inverse() const;
  
  /** Building a Matrix4 from orientation / scale / position.
  @remarks
      Transform is performed in the order scale, rotate, translation, i.e. translation is independent
      of orientation axes, scale does not affect size of translation, rotation and scaling are always
      centered on the origin.
  */
  void makeTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation);
  
  /** Building an inverse Matrix4 from orientation / scale / position.
  @remarks
      As makeTransform except it build the inverse given the same data as makeTransform, so
      performing -translation, -rotate, 1/scale in that order.
  */
  void makeInverseTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation);
  
  /** Decompose a Matrix4 to orientation / scale / position.
  */
  void decomposition(Vector3& position, Vector3& scale, Quaternion& orientation) const;
  
  /** Check whether or not the matrix is affine matrix.
      @remarks
          An affine matrix is a 4x4 matrix with row 3 equal to (0, 0, 0, 1),
          e.g. no projective coefficients.
  */
  inline bool isAffine(void) const
  {
    return m[3][0] == 0 && m[3][1] == 0 && m[3][2] == 0 && m[3][3] == 1;
  }
  
  /** Returns the inverse of the affine matrix.
      @note
          The matrix must be an affine matrix. @see Matrix4::isAffine.
  */
  Matrix4 inverseAffine(void) const;
  
  /** Concatenate two affine matrices.
      @note
          The matrices must be affine matrix. @see Matrix4::isAffine.
  */
  inline Matrix4 concatenateAffine(const Matrix4 &m2) const
  {
    assert(isAffine() && m2.isAffine());
    
    return Matrix4(
             m[0][0] * m2.m[0][0] + m[0][1] * m2.m[1][0] + m[0][2] * m2.m[2][0],
             m[0][0] * m2.m[0][1] + m[0][1] * m2.m[1][1] + m[0][2] * m2.m[2][1],
             m[0][0] * m2.m[0][2] + m[0][1] * m2.m[1][2] + m[0][2] * m2.m[2][2],
             m[0][0] * m2.m[0][3] + m[0][1] * m2.m[1][3] + m[0][2] * m2.m[2][3] + m[0][3],
             
             m[1][0] * m2.m[0][0] + m[1][1] * m2.m[1][0] + m[1][2] * m2.m[2][0],
             m[1][0] * m2.m[0][1] + m[1][1] * m2.m[1][1] + m[1][2] * m2.m[2][1],
             m[1][0] * m2.m[0][2] + m[1][1] * m2.m[1][2] + m[1][2] * m2.m[2][2],
             m[1][0] * m2.m[0][3] + m[1][1] * m2.m[1][3] + m[1][2] * m2.m[2][3] + m[1][3],
             
             m[2][0] * m2.m[0][0] + m[2][1] * m2.m[1][0] + m[2][2] * m2.m[2][0],
             m[2][0] * m2.m[0][1] + m[2][1] * m2.m[1][1] + m[2][2] * m2.m[2][1],
             m[2][0] * m2.m[0][2] + m[2][1] * m2.m[1][2] + m[2][2] * m2.m[2][2],
             m[2][0] * m2.m[0][3] + m[2][1] * m2.m[1][3] + m[2][2] * m2.m[2][3] + m[2][3],
             
             0, 0, 0, 1);
  }
  
  /** 3-D Vector transformation specially for an affine matrix.
      @remarks
          Transforms the given 3-D vector by the matrix, projecting the
          result back into <i>w</i> = 1.
      @note
          The matrix must be an affine matrix. @see Matrix4::isAffine.
  */
  inline Vector3 transformAffine(const Vector3& v) const
  {
    assert(isAffine());
    
    return Vector3(
             m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3],
             m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3],
             m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3]);
  }
  
};

} // namespace OgreMathLib

#endif