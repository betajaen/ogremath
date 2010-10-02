/*
-----------------------------------------------------------------------------

OgreMath.h

This source file was part of OGRE
(Object-oriented Graphics Rendering Engine)
For the latest info, see http://www.ogre3d.org/

Copyright (c) 2000-2009 Torus Knot Software Ltd
Re-arranged by Robin Southern 2010 http://github.com/betajaen/ogremath/

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
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

Some Code is is based on material from:

Geometric Tools, LLC
Copyright (c) 1998-2010
Distributed under the Boost Software License, Version 1.0.
http://www.boost.org/LICENSE_1_0.txt
http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

-----------------------------------------------------------------------------
*/

#ifndef ogremath_h
#define ogremath_h

#include <ostream>
#include <vector>
#include <list>
#include <map>
#include "assert.h"

#define OGREMATH_ALLOC_T(TYPE, SIZE, CATEGORY) (TYPE*) malloc(sizeof(TYPE) * SIZE)
#define OGREMATH_FREE(MEM, CATEGORY) free(MEM)

namespace OgreMath
{

 typedef float  Real;

 namespace Private
 {
  template<typename=void> class MathT;
  template<typename=void> class DegreeT;
  template<typename=void> class RadianT;
  template<typename=void> class Vector2T;
  template<typename=void> class Vector3T;
  template<typename=void> class Vector4T;
  template<typename=void> class QuaternionT;
  template<typename=void> class RayT;
  template<typename=void> class SphereT;
  template<typename=void> class PlaneT;
  template<typename=void> class AxisAlignedBoxT;
  template<typename=void> class Matrix3T;
  template<typename=void> class Matrix4T;
 };

 template <typename T> struct vector 
 { 
  typedef typename std::vector<typename T> type;
 }; 

 template <typename T> struct list 
 { 
  typedef typename std::list<typename T> type;
 };

 namespace Private
 {

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
  template<typename o_O> class MathT
  {
  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

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
   // angle units used by the api
   static AngleUnit msAngleUnit;

   /// Size of the trig tables as determined by constructor.
   static int mTrigTableSize;

   /// Radian -> index factor value ( mTrigTableSize / 2 * PI )
   static Real mTrigTableFactor;
   static Real* mSinTable;
   static Real* mTanTable;

   /** Private function to build trig tables.
   */
   void buildTrigTables()
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
   static Real SinTable (Real fValue)
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

   static Real TanTable (Real fValue)
   {
    // Convert range to index values, wrap if required
    int idx = int(fValue *= mTrigTableFactor) % mTrigTableSize;
    return mTanTable[idx];
   }

  public:
   /** Default constructor.
   @param
   trigTableSize Optional parameter to set the size of the
   tables used to implement Sin, Cos, Tan
   */
   MathT(unsigned int trigTableSize = 4096)
   {
    msAngleUnit = AU_DEGREE;

    mTrigTableSize = trigTableSize;
    mTrigTableFactor = mTrigTableSize / Math::TWO_PI;

    mSinTable = OGREMATH_ALLOC_T(Real, mTrigTableSize, MEMCATEGORY_GENERAL);
    mTanTable = OGREMATH_ALLOC_T(Real, mTrigTableSize, MEMCATEGORY_GENERAL);

    buildTrigTables();
   }

   /** Default destructor.
   */
   ~MathT()
   {
    OGREMATH_FREE(mSinTable, MEMCATEGORY_GENERAL);
    OGREMATH_FREE(mTanTable, MEMCATEGORY_GENERAL);
   }

   static inline int IAbs (int iValue) { return ( iValue >= 0 ? iValue : -iValue ); }
   static inline int ICeil (float fValue) { return int(ceil(fValue)); }
   static inline int IFloor (float fValue) { return int(floor(fValue)); }
   static int ISign (int iValue)
   {
    return ( iValue > 0 ? +1 : ( iValue < 0 ? -1 : 0 ) );
   }


   static inline Real Abs (Real fValue) { return Real(fabs(fValue)); }
   static inline Degree Abs (const Degree& dValue) { return Degree(fabs(dValue.valueDegrees())); }
   static inline Radian Abs (const Radian& rValue) { return Radian(fabs(rValue.valueRadians())); }
   static Radian ACos (Real fValue)
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

   static Radian ASin (Real fValue)
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
   static inline Radian ATan (Real fValue) { return Radian(atan(fValue)); }
   static inline Radian ATan2 (Real fY, Real fX) { return Radian(atan2(fY,fX)); }
   static inline Real Ceil (Real fValue) { return Real(ceil(fValue)); }
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

   static inline Real Exp (Real fValue) { return Real(exp(fValue)); }

   static inline Real Floor (Real fValue) { return Real(floor(fValue)); }

   static inline Real Log (Real fValue) { return Real(log(fValue)); }

   /// Stored value of log(2) for frequent use
   static const Real LOG2;

   static inline Real Log2 (Real fValue) { return Real(log(fValue)/LOG2); }

   static inline Real LogN (Real base, Real fValue) { return Real(log(fValue)/log(base)); }

   static inline Real Pow (Real fBase, Real fExponent) { return Real(pow(fBase,fExponent)); }

   static Real Sign (Real fValue)
   {
    if ( fValue > 0.0 )
     return 1.0;
    if ( fValue < 0.0 )
     return -1.0;
    return 0.0;
   }

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

   static inline Real Sqr (Real fValue) { return fValue*fValue; }

   static inline Real Sqrt (Real fValue) { return Real(sqrt(fValue)); }

   static inline Radian Sqrt (const Radian& fValue) { return Radian(sqrt(fValue.valueRadians())); }

   static inline Degree Sqrt (const Degree& fValue) { return Degree(sqrt(fValue.valueDegrees())); }

   /** Inverse square root i.e. 1 / Sqrt(x), good for vector
   normalisation.
   */
   static Real InvSqrt(Real fValue)
   {
    return Real(asm_rsq(fValue));
   }

   static Real UnitRandom ()   // in [0,1]
   {
    return asm_rand() / asm_rand_max();
   }

   static Real RangeRandom (Real fLow, Real fHigh)   // in [fLow,fHigh]
   {
    return (fHigh-fLow)*UnitRandom() + fLow;
   }

   static Real SymmetricRandom ()   // in [-1,1]
   {
    return 2.0f * UnitRandom() - 1.0f;
   }

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

   static inline Real DegreesToRadians(Real degrees) { return degrees * fDeg2Rad; }
   static inline Real RadiansToDegrees(Real radians) { return radians * fRad2Deg; }

   /** These functions used to set the assumed angle units (radians or degrees) 
   expected when using the Angle type.
   @par
   You can set this directly after creating a new Root, and also before/after resource creation,
   depending on whether you want the change to affect resource files.
   */
   static void setAngleUnit(AngleUnit unit)
   {
    msAngleUnit = unit;
   }

   /** Get the unit being used for angles. */
   static AngleUnit getAngleUnit(void)
   {
    return msAngleUnit;
   }

   /** Convert from the current AngleUnit to radians. */
   static Real AngleUnitsToRadians(Real units)
   {
    if (msAngleUnit == AU_DEGREE)
     return angleunits * fDeg2Rad;
    else
     return angleunits;
   }
   /** Convert from radians to the current AngleUnit . */
   static Real RadiansToAngleUnits(Real radians)
   {
    if (msAngleUnit == AU_DEGREE)
     return radians * fRad2Deg;
    else
     return radians;
   }

   /** Convert from the current AngleUnit to degrees. */
   static Real AngleUnitsToDegrees(Real units)
   {
    if (msAngleUnit == AU_RADIAN)
     return angleunits * fRad2Deg;
    else
     return angleunits;
   }

   /** Convert from degrees to the current AngleUnit. */
   static Real DegreesToAngleUnits(Real degrees)
   {
    if (msAngleUnit == AU_RADIAN)
     return degrees * fDeg2Rad;
    else
     return degrees;
   }

   /** Checks whether a given point is inside a triangle, in a
   2-dimensional (Cartesian) space.
   @remarks
   The vertices of the triangle must be given in either
   trigonometrical (anticlockwise) or inverse trigonometrical
   (clockwise) order.
   @param
   p The point.
   @param
   a The triangle's first vertex.
   @param
   b The triangle's second vertex.
   @param
   c The triangle's third vertex.
   @returns
   If the point resides in the triangle, <b>true</b> is
   returned.
   @par
   If the point is outside the triangle, <b>false</b> is
   returned.
   */
   static bool pointInTri2D(const Vector2& p, const Vector2& a, 
    const Vector2& b, const Vector2& c)
   {
    // Winding must be consistent from all edges for point to be inside
    Vector2 v1, v2;
    Real dot[3];
    bool zeroDot[3];

    v1 = b - a;
    v2 = p - a;

    // Note we don't care about normalisation here since sign is all we need
    // It means we don't have to worry about magnitude of cross products either
    dot[0] = v1.crossProduct(v2);
    zeroDot[0] = Math::RealEqual(dot[0], 0.0f, 1e-3);


    v1 = c - b;
    v2 = p - b;

    dot[1] = v1.crossProduct(v2);
    zeroDot[1] = Math::RealEqual(dot[1], 0.0f, 1e-3);

    // Compare signs (ignore colinear / coincident points)
    if(!zeroDot[0] && !zeroDot[1] 
    && Math::Sign(dot[0]) != Math::Sign(dot[1]))
    {
     return false;
    }

    v1 = a - c;
    v2 = p - c;

    dot[2] = v1.crossProduct(v2);
    zeroDot[2] = Math::RealEqual(dot[2], 0.0f, 1e-3);
    // Compare signs (ignore colinear / coincident points)
    if((!zeroDot[0] && !zeroDot[2] 
    && Math::Sign(dot[0]) != Math::Sign(dot[2])) ||
     (!zeroDot[1] && !zeroDot[2] 
    && Math::Sign(dot[1]) != Math::Sign(dot[2])))
    {
     return false;
    }


    return true;
   }

   /** Checks whether a given 3D point is inside a triangle.
   @remarks
   The vertices of the triangle must be given in either
   trigonometrical (anticlockwise) or inverse trigonometrical
   (clockwise) order, and the point must be guaranteed to be in the
   same plane as the triangle
   @param
   p The point.
   @param
   a The triangle's first vertex.
   @param
   b The triangle's second vertex.
   @param
   c The triangle's third vertex.
   @param 
   normal The triangle plane's normal (passed in rather than calculated
   on demand since the caller may already have it)
   @returns
   If the point resides in the triangle, <b>true</b> is
   returned.
   @par
   If the point is outside the triangle, <b>false</b> is
   returned.
   */
   static bool pointInTri3D(const Vector3& p, const Vector3& a, 
    const Vector3& b, const Vector3& c, const Vector3& normal)
   {
    // Winding must be consistent from all edges for point to be inside
    Vector3 v1, v2;
    Real dot[3];
    bool zeroDot[3];

    v1 = b - a;
    v2 = p - a;

    // Note we don't care about normalisation here since sign is all we need
    // It means we don't have to worry about magnitude of cross products either
    dot[0] = v1.crossProduct(v2).dotProduct(normal);
    zeroDot[0] = Math::RealEqual(dot[0], 0.0f, 1e-3);


    v1 = c - b;
    v2 = p - b;

    dot[1] = v1.crossProduct(v2).dotProduct(normal);
    zeroDot[1] = Math::RealEqual(dot[1], 0.0f, 1e-3);

    // Compare signs (ignore colinear / coincident points)
    if(!zeroDot[0] && !zeroDot[1] 
    && Math::Sign(dot[0]) != Math::Sign(dot[1]))
    {
     return false;
    }

    v1 = a - c;
    v2 = p - c;

    dot[2] = v1.crossProduct(v2).dotProduct(normal);
    zeroDot[2] = Math::RealEqual(dot[2], 0.0f, 1e-3);
    // Compare signs (ignore colinear / coincident points)
    if((!zeroDot[0] && !zeroDot[2] 
    && Math::Sign(dot[0]) != Math::Sign(dot[2])) ||
     (!zeroDot[1] && !zeroDot[2] 
    && Math::Sign(dot[1]) != Math::Sign(dot[2])))
    {
     return false;
    }


    return true;
   }


   /** Ray / plane intersection, returns boolean result and distance. */
   static std::pair<bool, Real> intersects(const Ray& ray, const Plane& plane)
   {
    Real denom = plane.normal.dotProduct(ray.getDirection());
    if (Math::Abs(denom) < std::numeric_limits<Real>::epsilon())
    {
     // Parallel
     return std::pair<bool, Real>(false, 0);
    }
    else
    {
     Real nom = plane.normal.dotProduct(ray.getOrigin()) + plane.d;
     Real t = -(nom/denom);
     return std::pair<bool, Real>(t >= 0, t);
    }

   }

   /** Ray / sphere intersection, returns boolean result and distance. */
   static std::pair<bool, Real> intersects(const Ray& ray, const Sphere& sphere, 
    bool discardInside = true)
   {
    const Vector3& raydir = ray.getDirection();
    // Adjust ray origin relative to sphere center
    const Vector3& rayorig = ray.getOrigin() - sphere.getCenter();
    Real radius = sphere.getRadius();

    // Check origin inside first
    if (rayorig.squaredLength() <= radius*radius && discardInside)
    {
     return std::pair<bool, Real>(true, 0);
    }

    // Mmm, quadratics
    // Build coeffs which can be used with std quadratic solver
    // ie t = (-b +/- sqrt(b*b + 4ac)) / 2a
    Real a = raydir.dotProduct(raydir);
    Real b = 2 * rayorig.dotProduct(raydir);
    Real c = rayorig.dotProduct(rayorig) - radius*radius;

    // Calc determinant
    Real d = (b*b) - (4 * a * c);
    if (d < 0)
    {
     // No intersection
     return std::pair<bool, Real>(false, 0);
    }
    else
    {
     // BTW, if d=0 there is one intersection, if d > 0 there are 2
     // But we only want the closest one, so that's ok, just use the 
     // '-' version of the solver
     Real t = ( -b - Math::Sqrt(d) ) / (2 * a);
     if (t < 0)
      t = ( -b + Math::Sqrt(d) ) / (2 * a);
     return std::pair<bool, Real>(true, t);
    }
   }

   /** Ray / box intersection, returns boolean result and distance. */
   static std::pair<bool, Real> intersects(const Ray& ray, const AxisAlignedBox& box)
   {
    if (box.isNull()) return std::pair<bool, Real>(false, 0);
    if (box.isInfinite()) return std::pair<bool, Real>(true, 0);

    Real lowt = 0.0f;
    Real t;
    bool hit = false;
    Vector3 hitpoint;
    const Vector3& min = box.getMinimum();
    const Vector3& max = box.getMaximum();
    const Vector3& rayorig = ray.getOrigin();
    const Vector3& raydir = ray.getDirection();

    // Check origin inside first
    if ( rayorig > min && rayorig < max )
    {
     return std::pair<bool, Real>(true, 0);
    }

    // Check each face in turn, only check closest 3
    // Min x
    if (rayorig.x <= min.x && raydir.x > 0)
    {
     t = (min.x - rayorig.x) / raydir.x;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
       hitpoint.z >= min.z && hitpoint.z <= max.z &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }
    // Max x
    if (rayorig.x >= max.x && raydir.x < 0)
    {
     t = (max.x - rayorig.x) / raydir.x;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
       hitpoint.z >= min.z && hitpoint.z <= max.z &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }
    // Min y
    if (rayorig.y <= min.y && raydir.y > 0)
    {
     t = (min.y - rayorig.y) / raydir.y;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
       hitpoint.z >= min.z && hitpoint.z <= max.z &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }
    // Max y
    if (rayorig.y >= max.y && raydir.y < 0)
    {
     t = (max.y - rayorig.y) / raydir.y;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
       hitpoint.z >= min.z && hitpoint.z <= max.z &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }
    // Min z
    if (rayorig.z <= min.z && raydir.z > 0)
    {
     t = (min.z - rayorig.z) / raydir.z;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
       hitpoint.y >= min.y && hitpoint.y <= max.y &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }
    // Max z
    if (rayorig.z >= max.z && raydir.z < 0)
    {
     t = (max.z - rayorig.z) / raydir.z;
     if (t >= 0)
     {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
       hitpoint.y >= min.y && hitpoint.y <= max.y &&
       (!hit || t < lowt))
      {
       hit = true;
       lowt = t;
      }
     }
    }

    return std::pair<bool, Real>(hit, lowt);

   }

   /** Ray / box intersection, returns boolean result and two intersection distance.
   @param
   ray The ray.
   @param
   box The box.
   @param
   d1 A real pointer to retrieve the near intersection distance
   from the ray origin, maybe <b>null</b> which means don't care
   about the near intersection distance.
   @param
   d2 A real pointer to retrieve the far intersection distance
   from the ray origin, maybe <b>null</b> which means don't care
   about the far intersection distance.
   @returns
   If the ray is intersects the box, <b>true</b> is returned, and
   the near intersection distance is return by <i>d1</i>, the
   far intersection distance is return by <i>d2</i>. Guarantee
   <b>0</b> <= <i>d1</i> <= <i>d2</i>.
   @par
   If the ray isn't intersects the box, <b>false</b> is returned, and
   <i>d1</i> and <i>d2</i> is unmodified.
   */
   static bool intersects(const Ray& ray, const AxisAlignedBox& box,
    Real* d1, Real* d2)
   {
    if (box.isNull())
     return false;

    if (box.isInfinite())
    {
     if (d1) *d1 = 0;
     if (d2) *d2 = Math::POS_INFINITY;
     return true;
    }

    const Vector3& min = box.getMinimum();
    const Vector3& max = box.getMaximum();
    const Vector3& rayorig = ray.getOrigin();
    const Vector3& raydir = ray.getDirection();

    Vector3 absDir;
    absDir[0] = Math::Abs(raydir[0]);
    absDir[1] = Math::Abs(raydir[1]);
    absDir[2] = Math::Abs(raydir[2]);

    // Sort the axis, ensure check minimise floating error axis first
    int imax = 0, imid = 1, imin = 2;
    if (absDir[0] < absDir[2])
    {
     imax = 2;
     imin = 0;
    }
    if (absDir[1] < absDir[imin])
    {
     imid = imin;
     imin = 1;
    }
    else if (absDir[1] > absDir[imax])
    {
     imid = imax;
     imax = 1;
    }

    Real start = 0, end = Math::POS_INFINITY;

#define _CALC_AXIS(i)                                       \
 do {                                                    \
 Real denom = 1 / raydir[i];                         \
 Real newstart = (min[i] - rayorig[i]) * denom;      \
 Real newend = (max[i] - rayorig[i]) * denom;        \
 if (newstart > newend) std::swap(newstart, newend); \
 if (newstart > end || newend < start) return false; \
 if (newstart > start) start = newstart;             \
 if (newend < end) end = newend;                     \
 } while(0)

    // Check each axis in turn

    _CALC_AXIS(imax);

    if (absDir[imid] < std::numeric_limits<Real>::epsilon())
    {
     // Parallel with middle and minimise axis, check bounds only
     if (rayorig[imid] < min[imid] || rayorig[imid] > max[imid] ||
      rayorig[imin] < min[imin] || rayorig[imin] > max[imin])
      return false;
    }
    else
    {
     _CALC_AXIS(imid);

     if (absDir[imin] < std::numeric_limits<Real>::epsilon())
     {
      // Parallel with minimise axis, check bounds only
      if (rayorig[imin] < min[imin] || rayorig[imin] > max[imin])
       return false;
     }
     else
     {
      _CALC_AXIS(imin);
     }
    }
#undef _CALC_AXIS

    if (d1) *d1 = start;
    if (d2) *d2 = end;

    return true;

   }

   /** Ray / triangle intersection, returns boolean result and distance.
   @param
   ray The ray.
   @param
   a The triangle's first vertex.
   @param
   b The triangle's second vertex.
   @param
   c The triangle's third vertex.
   @param 
   normal The triangle plane's normal (passed in rather than calculated
   on demand since the caller may already have it), doesn't need
   normalised since we don't care.
   @param
   positiveSide Intersect with "positive side" of the triangle
   @param
   negativeSide Intersect with "negative side" of the triangle
   @returns
   If the ray is intersects the triangle, a pair of <b>true</b> and the
   distance between intersection point and ray origin returned.
   @par
   If the ray isn't intersects the triangle, a pair of <b>false</b> and
   <b>0</b> returned.
   */
   static std::pair<bool, Real> intersects(const Ray& ray, const Vector3& a,
    const Vector3& b, const Vector3& c, const Vector3& normal,
    bool positiveSide = true, bool negativeSide = true)
   {
    //
    // Calculate intersection with plane.
    //
    Real t;
    {
     Real denom = normal.dotProduct(ray.getDirection());

     // Check intersect side
     if (denom > + std::numeric_limits<Real>::epsilon())
     {
      if (!negativeSide)
       return std::pair<bool, Real>(false, 0);
     }
     else if (denom < - std::numeric_limits<Real>::epsilon())
     {
      if (!positiveSide)
       return std::pair<bool, Real>(false, 0);
     }
     else
     {
      // Parallel or triangle area is close to zero when
      // the plane normal not normalised.
      return std::pair<bool, Real>(false, 0);
     }

     t = normal.dotProduct(a - ray.getOrigin()) / denom;

     if (t < 0)
     {
      // Intersection is behind origin
      return std::pair<bool, Real>(false, 0);
     }
    }

    //
    // Calculate the largest area projection plane in X, Y or Z.
    //
    size_t i0, i1;
    {
     Real n0 = Math::Abs(normal[0]);
     Real n1 = Math::Abs(normal[1]);
     Real n2 = Math::Abs(normal[2]);

     i0 = 1; i1 = 2;
     if (n1 > n2)
     {
      if (n1 > n0) i0 = 0;
     }
     else
     {
      if (n2 > n0) i1 = 0;
     }
    }

    //
    // Check the intersection point is inside the triangle.
    //
    {
     Real u1 = b[i0] - a[i0];
     Real v1 = b[i1] - a[i1];
     Real u2 = c[i0] - a[i0];
     Real v2 = c[i1] - a[i1];
     Real u0 = t * ray.getDirection()[i0] + ray.getOrigin()[i0] - a[i0];
     Real v0 = t * ray.getDirection()[i1] + ray.getOrigin()[i1] - a[i1];

     Real alpha = u0 * v2 - u2 * v0;
     Real beta  = u1 * v0 - u0 * v1;
     Real area  = u1 * v2 - u2 * v1;

     // epsilon to avoid float precision error
     const Real EPSILON = 1e-6f;

     Real tolerance = - EPSILON * area;

     if (area > 0)
     {
      if (alpha < tolerance || beta < tolerance || alpha+beta > area-tolerance)
       return std::pair<bool, Real>(false, 0);
     }
     else
     {
      if (alpha > tolerance || beta > tolerance || alpha+beta < area-tolerance)
       return std::pair<bool, Real>(false, 0);
     }
    }

    return std::pair<bool, Real>(true, t);

   }

   /** Ray / triangle intersection, returns boolean result and distance.
   @param
   ray The ray.
   @param
   a The triangle's first vertex.
   @param
   b The triangle's second vertex.
   @param
   c The triangle's third vertex.
   @param
   positiveSide Intersect with "positive side" of the triangle
   @param
   negativeSide Intersect with "negative side" of the triangle
   @returns
   If the ray is intersects the triangle, a pair of <b>true</b> and the
   distance between intersection point and ray origin returned.
   @par
   If the ray isn't intersects the triangle, a pair of <b>false</b> and
   <b>0</b> returned.
   */
   static std::pair<bool, Real> intersects(const Ray& ray, const Vector3& a,
    const Vector3& b, const Vector3& c,
    bool positiveSide = true, bool negativeSide = true)
   {
    Vector3 normal = calculateBasicFaceNormalWithoutNormalize(a, b, c);
    return intersects(ray, a, b, c, normal, positiveSide, negativeSide);
   }

   /** Sphere / box intersection test. */
   static bool intersects(const Sphere& sphere, const AxisAlignedBox& box)
   {
    if (box.isNull()) return false;
    if (box.isInfinite()) return true;

    // Use splitting planes
    const Vector3& center = sphere.getCenter();
    Real radius = sphere.getRadius();
    const Vector3& min = box.getMinimum();
    const Vector3& max = box.getMaximum();

    // Arvo's algorithm
    Real s, d = 0;
    for (int i = 0; i < 3; ++i)
    {
     if (center.ptr()[i] < min.ptr()[i])
     {
      s = center.ptr()[i] - min.ptr()[i];
      d += s * s; 
     }
     else if(center.ptr()[i] > max.ptr()[i])
     {
      s = center.ptr()[i] - max.ptr()[i];
      d += s * s; 
     }
    }
    return d <= radius * radius;

   }

   /** Plane / box intersection test. */
   static bool intersects(const Plane& plane, const AxisAlignedBox& box)
   {
    return (plane.getSide(box) == Plane::BOTH_SIDE);
   }

   /** Ray / convex plane list intersection test. 
   @param ray The ray to test with
   @param plaeList List of planes which form a convex volume
   @param normalIsOutside Does the normal point outside the volume
   */
   static std::pair<bool, Real> intersects(
    const Ray& ray, const typename vector<Plane>::type& planeList, 
    bool normalIsOutside)
   {
    list<Plane>::type planesList;
    for (vector<Plane>::type::const_iterator i = planes.begin(); i != planes.end(); ++i)
    {
     planesList.push_back(*i);
    }
    return intersects(ray, planesList, normalIsOutside);
   }

   /** Ray / convex plane list intersection test. 
   @param ray The ray to test with
   @param plaeList List of planes which form a convex volume
   @param normalIsOutside Does the normal point outside the volume
   */
   static std::pair<bool, Real> intersects(
    const Ray& ray, const typename list<Plane>::type& planeList, 
    bool normalIsOutside)
   {
    list<Plane>::type::const_iterator planeit, planeitend;
    planeitend = planes.end();
    bool allInside = true;
    std::pair<bool, Real> ret;
    std::pair<bool, Real> end;
    ret.first = false;
    ret.second = 0.0f;
    end.first = false;
    end.second = 0;


    // derive side
    // NB we don't pass directly since that would require Plane::Side in 
    // interface, which results in recursive includes since Math is so fundamental
    Plane::Side outside = normalIsOutside ? Plane::POSITIVE_SIDE : Plane::NEGATIVE_SIDE;

    for (planeit = planes.begin(); planeit != planeitend; ++planeit)
    {
     const Plane& plane = *planeit;
     // is origin outside?
     if (plane.getSide(ray.getOrigin()) == outside)
     {
      allInside = false;
      // Test single plane
      std::pair<bool, Real> planeRes = 
       ray.intersects(plane);
      if (planeRes.first)
      {
       // Ok, we intersected
       ret.first = true;
       // Use the most distant result since convex volume
       ret.second = std::max(ret.second, planeRes.second);
      }
      else
      {
       ret.first =false;
       ret.second=0.0f;
       return ret;
      }
     }
     else
     {
      std::pair<bool, Real> planeRes = 
       ray.intersects(plane);
      if (planeRes.first)
      {
       if( !end.first )
       {
        end.first = true;
        end.second = planeRes.second;
       }
       else
       {
        end.second = std::min( planeRes.second, end.second );
       }

      }

     }
    }

    if (allInside)
    {
     // Intersecting at 0 distance since inside the volume!
     ret.first = true;
     ret.second = 0.0f;
     return ret;
    }

    if( end.first )
    {
     if( end.second < ret.second )
     {
      ret.first = false;
      return ret;
     }
    }
    return ret;

   }

   /** Sphere / plane intersection test. 
   @remarks NB just do a plane.getDistance(sphere.getCenter()) for more detail!
   */
   static bool intersects(const Sphere& sphere, const Plane& plane)
   {
    return (
     Math::Abs(plane.getDistance(sphere.getCenter()))
     <= sphere.getRadius() );
   }

   /** Compare 2 reals, using tolerance for inaccuracies.
   */
   static bool RealEqual(Real a, Real b,
    Real tolerance = std::numeric_limits<Real>::epsilon())
   {
    if (fabs(b-a) <= tolerance)
     return true;
    else
     return false;
   }

   /** Calculates the tangent space vector for a given set of positions / texture coords. */
   static Vector3 calculateTangentSpaceVector(
    const Vector3& position1, const Vector3& position2, const Vector3& position3,
    Real u1, Real v1, Real u2, Real v2, Real u3, Real v3);

   /** Build a reflection matrix for the passed in plane. */
   static Matrix4 buildReflectionMatrix(const Plane& p)
   {
    return Matrix4(
     -2 * p.normal.x * p.normal.x + 1,   -2 * p.normal.x * p.normal.y,       -2 * p.normal.x * p.normal.z,       -2 * p.normal.x * p.d, 
     -2 * p.normal.y * p.normal.x,       -2 * p.normal.y * p.normal.y + 1,   -2 * p.normal.y * p.normal.z,       -2 * p.normal.y * p.d, 
     -2 * p.normal.z * p.normal.x,       -2 * p.normal.z * p.normal.y,       -2 * p.normal.z * p.normal.z + 1,   -2 * p.normal.z * p.d, 
     0,                                  0,                                  0,                                  1);
   }
   /** Calculate a face normal, including the w component which is the offset from the origin. */
   static Vector4 calculateFaceNormal(const Vector3& v1, const Vector3& v2, const Vector3& v3)
   {
    Vector3 normal = calculateBasicFaceNormal(v1, v2, v3);
    // Now set up the w (distance of tri from origin
    return Vector4(normal.x, normal.y, normal.z, -(normal.dotProduct(v1)));
   }

   /** Calculate a face normal, no w-information. */
   static Vector3 calculateBasicFaceNormal(const Vector3& v1, const Vector3& v2, const Vector3& v3)
   {
    Vector3 normal = (v2 - v1).crossProduct(v3 - v1);
    normal.normalise();
    return normal;
   }
   /** Calculate a face normal without normalize, including the w component which is the offset from the origin. */
   static Vector4 calculateFaceNormalWithoutNormalize(const Vector3& v1, const Vector3& v2, const Vector3& v3)
   {
    Vector3 normal = calculateBasicFaceNormalWithoutNormalize(v1, v2, v3);
    // Now set up the w (distance of tri from origin)
    return Vector4(normal.x, normal.y, normal.z, -(normal.dotProduct(v1)));
   }

   /** Calculate a face normal without normalize, no w-information. */
   static Vector3 calculateBasicFaceNormalWithoutNormalize(const Vector3& v1, const Vector3& v2, const Vector3& v3)
   {
    Vector3 normal = (v2 - v1).crossProduct(v3 - v1);
    return normal;
   }

   /** Generates a value based on the Gaussian (normal) distribution function
   with the given offset and scale parameters.
   */
   static Real gaussianDistribution(Real x, Real offset = 0.0f, Real scale = 1.0f)
   {
    Real nom = Math::Exp(
     -Math::Sqr(x - offset) / (2 * Math::Sqr(scale)));
    Real denom = scale * Math::Sqrt(2 * Math::PI);

    return nom / denom;

   }

   /** Clamp a value within an inclusive range. */
   template <typename T>
   static T Clamp(T val, T minval, T maxval)
   {
    assert (minval < maxval && "Invalid clamp range");
    return std::max(std::min(val, maxval), minval);
   }

   static Matrix4 makeViewMatrix(const Vector3& position, const Quaternion& orientation, 
    const Matrix4* reflectMatrix = 0)
   {
    Matrix4 viewMatrix;

    // View matrix is:
    //
    //  [ Lx  Uy  Dz  Tx  ]
    //  [ Lx  Uy  Dz  Ty  ]
    //  [ Lx  Uy  Dz  Tz  ]
    //  [ 0   0   0   1   ]
    //
    // Where T = -(Transposed(Rot) * Pos)

    // This is most efficiently done using 3x3 Matrices
    Matrix3 rot;
    orientation.ToRotationMatrix(rot);

    // Make the translation relative to new axes
    Matrix3 rotT = rot.Transpose();
    Vector3 trans = -rotT * position;

    // Make final matrix
    viewMatrix = Matrix4::IDENTITY;
    viewMatrix = rotT; // fills upper 3x3
    viewMatrix[0][3] = trans.x;
    viewMatrix[1][3] = trans.y;
    viewMatrix[2][3] = trans.z;

    // Deal with reflections
    if (reflectMatrix)
    {
     viewMatrix = viewMatrix * (*reflectMatrix);
    }

    return viewMatrix;


   }

   /** Get a bounding radius value from a bounding box. */
   static Real boundingRadiusFromAABB(const AxisAlignedBox& aabb)
   {
    Vector3 max = aabb.getMaximum();
    Vector3 min = aabb.getMinimum();

    Vector3 magnitude = max;
    magnitude.makeCeil(-max);
    magnitude.makeCeil(min);
    magnitude.makeCeil(-min);

    return magnitude.length();
   }



   static const Real POS_INFINITY;
   static const Real NEG_INFINITY;
   static const Real PI;
   static const Real TWO_PI;
   static const Real HALF_PI;
   static const Real fDeg2Rad;
   static const Real fRad2Deg;

  }; // class Math



  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Wrapper class which indicates a given angle value is in Radians.
  @remarks
  Radian values are interchangeable with Degree values, and conversions
  will be done automatically between them.
  */
  template<typename o_O> class RadianT
  {
   Real mRad;

  public:

   
   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   
  // these functions could not be defined within the class definition of class
  // Radian because they required class Degree to be defined
  inline Radian ( const Degree& d ) : mRad(d.valueRadians()) {
  }
  inline Radian& operator = ( const Degree& d ) {
   mRad = d.valueRadians(); return *this;
  }
  inline Radian operator + ( const Degree& d ) const {
   return Radian ( mRad + d.valueRadians() );
  }
  inline Radian& operator += ( const Degree& d ) {
   mRad += d.valueRadians();
   return *this;
  }
  inline Radian operator - ( const Degree& d ) const {
   return Radian ( mRad - d.valueRadians() );
  }
  inline Radian& operator -= ( const Degree& d ) {
   mRad -= d.valueRadians();
   return *this;
  }
   explicit Radian ( Real r=0 ) : mRad(r) {}

   Radian& operator = ( const Real& f ) { mRad = f; return *this; }
   Radian& operator = ( const Radian& r ) { mRad = r.mRad; return *this; }

   Real valueDegrees() const; // see bottom of this file
   Real valueRadians() const { return mRad; }
   Real valueAngleUnits() const;

   const Radian& operator + () const { return *this; }
   Radian operator + ( const Radian& r ) const { return Radian ( mRad + r.mRad ); }
   Radian& operator += ( const Radian& r ) { mRad += r.mRad; return *this; }
   Radian operator - () const { return Radian(-mRad); }
   Radian operator - ( const Radian& r ) const { return Radian ( mRad - r.mRad ); }
   Radian& operator -= ( const Radian& r ) { mRad -= r.mRad; return *this; }
   Radian operator * ( Real f ) const { return Radian ( mRad * f ); }
   Radian operator * ( const Radian& f ) const { return Radian ( mRad * f.mRad ); }
   Radian& operator *= ( Real f ) { mRad *= f; return *this; }
   Radian operator / ( Real f ) const { return Radian ( mRad / f ); }
   Radian& operator /= ( Real f ) { mRad /= f; return *this; }

   bool operator <  ( const Radian& r ) const { return mRad <  r.mRad; }
   bool operator <= ( const Radian& r ) const { return mRad <= r.mRad; }
   bool operator == ( const Radian& r ) const { return mRad == r.mRad; }
   bool operator != ( const Radian& r ) const { return mRad != r.mRad; }
   bool operator >= ( const Radian& r ) const { return mRad >= r.mRad; }
   bool operator >  ( const Radian& r ) const { return mRad >  r.mRad; }

   inline friend std::ostream& operator <<
    ( std::ostream& o, const Radian& v )
   {
    o << "Radian(" << v.valueRadians() << ")";
    return o;
   }
  };

  /** Wrapper class which indicates a given angle value is in Degrees.
  @remarks
  Degree values are interchangeable with Radian values, and conversions
  will be done automatically between them.
  */
  template<typename o_O> class DegreeT
  {
   Real mDeg; // if you get an error here - make sure to define/typedef 'Real' first

  public:
   
   
   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;
   
   explicit Degree ( Real d=0 ) : mDeg(d) {}
   Degree ( const Radian& r ) : mDeg(r.valueDegrees()) {}
   Degree& operator = ( const Real& f ) { mDeg = f; return *this; }
   Degree& operator = ( const Degree& d ) { mDeg = d.mDeg; return *this; }
   Degree& operator = ( const Radian& r ) { mDeg = r.valueDegrees(); return *this; }

   Real valueDegrees() const { return mDeg; }
   Real valueRadians() const; // see bottom of this file
   Real valueAngleUnits() const;

   const Degree& operator + () const { return *this; }
   Degree operator + ( const Degree& d ) const { return Degree ( mDeg + d.mDeg ); }
   Degree operator + ( const Radian& r ) const { return Degree ( mDeg + r.valueDegrees() ); }
   Degree& operator += ( const Degree& d ) { mDeg += d.mDeg; return *this; }
   Degree& operator += ( const Radian& r ) { mDeg += r.valueDegrees(); return *this; }
   Degree operator - () const { return Degree(-mDeg); }
   Degree operator - ( const Degree& d ) const { return Degree ( mDeg - d.mDeg ); }
   Degree operator - ( const Radian& r ) const { return Degree ( mDeg - r.valueDegrees() ); }
   Degree& operator -= ( const Degree& d ) { mDeg -= d.mDeg; return *this; }
   Degree& operator -= ( const Radian& r ) { mDeg -= r.valueDegrees(); return *this; }
   Degree operator * ( Real f ) const { return Degree ( mDeg * f ); }
   Degree operator * ( const Degree& f ) const { return Degree ( mDeg * f.mDeg ); }
   Degree& operator *= ( Real f ) { mDeg *= f; return *this; }
   Degree operator / ( Real f ) const { return Degree ( mDeg / f ); }
   Degree& operator /= ( Real f ) { mDeg /= f; return *this; }

   bool operator <  ( const Degree& d ) const { return mDeg <  d.mDeg; }
   bool operator <= ( const Degree& d ) const { return mDeg <= d.mDeg; }
   bool operator == ( const Degree& d ) const { return mDeg == d.mDeg; }
   bool operator != ( const Degree& d ) const { return mDeg != d.mDeg; }
   bool operator >= ( const Degree& d ) const { return mDeg >= d.mDeg; }
   bool operator >  ( const Degree& d ) const { return mDeg >  d.mDeg; }

   inline friend std::ostream& operator <<
    ( std::ostream& o, const Degree& v )
   {
    o << "Degree(" << v.valueDegrees() << ")";
    return o;
   }
  };

  /** Wrapper class which identifies a value as the currently default angle 
  type, as defined by Math::setAngleUnit.
  @remarks
  Angle values will be automatically converted between radians and degrees,
  as appropriate.
  */
  class Angle
  {
   Real mAngle;
   
   typedef Private::MathT<>            Math;
   typedef Private::DegreeT<>          Degree;
   typedef Private::RadianT<>          Radian;
   typedef Private::Vector2T<>         Vector2;
   typedef Private::Vector3T<>         Vector3;
   typedef Private::Vector4T<>         Vector4;
   typedef Private::QuaternionT<>      Quaternion;
   typedef Private::RayT<>             Ray;
   typedef Private::SphereT<>          Sphere;
   typedef Private::PlaneT<>           Plane;
   typedef Private::AxisAlignedBoxT<>  AxisAlignedBox;
   typedef Private::Matrix3T<>         Matrix3;
   typedef Private::Matrix4T<>         Matrix4;

  public:
   explicit Angle ( Real angle ) : mAngle(angle) {}
   operator Radian() const;
   operator Degree() const;
  };







  // -------------------------------------------------------------------------------





  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Standard 2-dimensional vector.
  @remarks
  A direction in 2D space represented as distances along the 2
  orthogonal axes (x, y). Note that positions, directions and
  scaling factors can be represented by a vector, depending on how
  you interpret the values.
  */
  template<typename o_O> class Vector2T
  {

  public:

   Real x, y;

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   inline Vector2()
   {
   }

   inline Vector2(const Real fX, const Real fY )
    : x( fX ), y( fY )
   {
   }

   inline explicit Vector2( const Real scaler )
    : x( scaler), y( scaler )
   {
   }

   inline explicit Vector2( const Real afCoordinate[2] )
    : x( afCoordinate[0] ),
    y( afCoordinate[1] )
   {
   }

   inline explicit Vector2( const int afCoordinate[2] )
   {
    x = (Real)afCoordinate[0];
    y = (Real)afCoordinate[1];
   }

   inline explicit Vector2( Real* const r )
    : x( r[0] ), y( r[1] )
   {
   }

   /** Exchange the contents of this vector with another. 
   */
   inline void swap(Vector2& other)
   {
    std::swap(x, other.x);
    std::swap(y, other.y);
   }

   inline Real operator [] ( const size_t i ) const
   {
    assert( i < 2 );

    return *(&x+i);
   }

   inline Real& operator [] ( const size_t i )
   {
    assert( i < 2 );

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
   inline Vector2& operator = ( const Vector2& rkVector )
   {
    x = rkVector.x;
    y = rkVector.y;

    return *this;
   }

   inline Vector2& operator = ( const Real fScalar)
   {
    x = fScalar;
    y = fScalar;

    return *this;
   }

   inline bool operator == ( const Vector2& rkVector ) const
   {
    return ( x == rkVector.x && y == rkVector.y );
   }

   inline bool operator != ( const Vector2& rkVector ) const
   {
    return ( x != rkVector.x || y != rkVector.y  );
   }

   // arithmetic operations
   inline Vector2 operator + ( const Vector2& rkVector ) const
   {
    return Vector2(
     x + rkVector.x,
     y + rkVector.y);
   }

   inline Vector2 operator - ( const Vector2& rkVector ) const
   {
    return Vector2(
     x - rkVector.x,
     y - rkVector.y);
   }

   inline Vector2 operator * ( const Real fScalar ) const
   {
    return Vector2(
     x * fScalar,
     y * fScalar);
   }

   inline Vector2 operator * ( const Vector2& rhs) const
   {
    return Vector2(
     x * rhs.x,
     y * rhs.y);
   }

   inline Vector2 operator / ( const Real fScalar ) const
   {
    assert( fScalar != 0.0 );

    Real fInv = 1.0f / fScalar;

    return Vector2(
     x * fInv,
     y * fInv);
   }

   inline Vector2 operator / ( const Vector2& rhs) const
   {
    return Vector2(
     x / rhs.x,
     y / rhs.y);
   }

   inline const Vector2& operator + () const
   {
    return *this;
   }

   inline Vector2 operator - () const
   {
    return Vector2(-x, -y);
   }

   // overloaded operators to help Vector2
   inline friend Vector2 operator * ( const Real fScalar, const Vector2& rkVector )
   {
    return Vector2(
     fScalar * rkVector.x,
     fScalar * rkVector.y);
   }

   inline friend Vector2 operator / ( const Real fScalar, const Vector2& rkVector )
   {
    return Vector2(
     fScalar / rkVector.x,
     fScalar / rkVector.y);
   }

   inline friend Vector2 operator + (const Vector2& lhs, const Real rhs)
   {
    return Vector2(
     lhs.x + rhs,
     lhs.y + rhs);
   }

   inline friend Vector2 operator + (const Real lhs, const Vector2& rhs)
   {
    return Vector2(
     lhs + rhs.x,
     lhs + rhs.y);
   }

   inline friend Vector2 operator - (const Vector2& lhs, const Real rhs)
   {
    return Vector2(
     lhs.x - rhs,
     lhs.y - rhs);
   }

   inline friend Vector2 operator - (const Real lhs, const Vector2& rhs)
   {
    return Vector2(
     lhs - rhs.x,
     lhs - rhs.y);
   }
   // arithmetic updates
   inline Vector2& operator += ( const Vector2& rkVector )
   {
    x += rkVector.x;
    y += rkVector.y;

    return *this;
   }

   inline Vector2& operator += ( const Real fScaler )
   {
    x += fScaler;
    y += fScaler;

    return *this;
   }

   inline Vector2& operator -= ( const Vector2& rkVector )
   {
    x -= rkVector.x;
    y -= rkVector.y;

    return *this;
   }

   inline Vector2& operator -= ( const Real fScaler )
   {
    x -= fScaler;
    y -= fScaler;

    return *this;
   }

   inline Vector2& operator *= ( const Real fScalar )
   {
    x *= fScalar;
    y *= fScalar;

    return *this;
   }

   inline Vector2& operator *= ( const Vector2& rkVector )
   {
    x *= rkVector.x;
    y *= rkVector.y;

    return *this;
   }

   inline Vector2& operator /= ( const Real fScalar )
   {
    assert( fScalar != 0.0 );

    Real fInv = 1.0f / fScalar;

    x *= fInv;
    y *= fInv;

    return *this;
   }

   inline Vector2& operator /= ( const Vector2& rkVector )
   {
    x /= rkVector.x;
    y /= rkVector.y;

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
    return Math::Sqrt( x * x + y * y );
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
    return x * x + y * y;
   }
   /** Returns the distance to another vector.
   @warning
   This operation requires a square root and is expensive in
   terms of CPU operations. If you don't need to know the exact
   distance (e.g. for just comparing distances) use squaredDistance()
   instead.
   */
   inline Real distance(const Vector2& rhs) const
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
   inline Real squaredDistance(const Vector2& rhs) const
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
   inline Real dotProduct(const Vector2& vec) const
   {
    return x * vec.x + y * vec.y;
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
    Real fLength = Math::Sqrt( x * x + y * y);

    // Will also work for zero-sized vectors, but will change nothing
    if ( fLength > 1e-08 )
    {
     Real fInvLength = 1.0f / fLength;
     x *= fInvLength;
     y *= fInvLength;
    }

    return fLength;
   }



   /** Returns a vector at a point half way between this and the passed
   in vector.
   */
   inline Vector2 midPoint( const Vector2& vec ) const
   {
    return Vector2(
     ( x + vec.x ) * 0.5f,
     ( y + vec.y ) * 0.5f );
   }

   /** Returns true if the vector's scalar components are all greater
   that the ones of the vector it is compared against.
   */
   inline bool operator < ( const Vector2& rhs ) const
   {
    if( x < rhs.x && y < rhs.y )
     return true;
    return false;
   }

   /** Returns true if the vector's scalar components are all smaller
   that the ones of the vector it is compared against.
   */
   inline bool operator > ( const Vector2& rhs ) const
   {
    if( x > rhs.x && y > rhs.y )
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
   inline void makeFloor( const Vector2& cmp )
   {
    if( cmp.x < x ) x = cmp.x;
    if( cmp.y < y ) y = cmp.y;
   }

   /** Sets this vector's components to the maximum of its own and the
   ones of the passed in vector.
   @remarks
   'Maximum' in this case means the combination of the highest
   value of x, y and z from both vectors. Highest is taken just
   numerically, not magnitude, so 1 > -3.
   */
   inline void makeCeil( const Vector2& cmp )
   {
    if( cmp.x > x ) x = cmp.x;
    if( cmp.y > y ) y = cmp.y;
   }

   /** Generates a vector perpendicular to this vector (eg an 'up' vector).
   @remarks
   This method will return a vector which is perpendicular to this
   vector. There are an infinite number of possibilities but this
   method will guarantee to generate one of them. If you need more
   control you should use the Quaternion class.
   */
   inline Vector2 perpendicular(void) const
   {
    return Vector2 (-y, x);
   }
   /** Calculates the 2 dimensional cross-product of 2 vectors, which results
   in a single floating point value which is 2 times the area of the triangle.
   */
   inline Real crossProduct( const Vector2& rkVector ) const
   {
    return x * rkVector.y - y * rkVector.x;
   }
   /** Generates a new random vector which deviates from this vector by a
   given angle in a random direction.
   @remarks
   This method assumes that the random number generator has already
   been seeded appropriately.
   @param
   angle The angle at which to deviate in radians
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
   inline Vector2 randomDeviant(
    Real angle) const
   {

    angle *=  Math::UnitRandom() * Math::TWO_PI;
    Real cosa = cos(angle);
    Real sina = sin(angle);
    return  Vector2(cosa * x - sina * y,
     sina * x + cosa * y);
   }

   /** Returns true if this vector is zero length. */
   inline bool isZeroLength(void) const
   {
    Real sqlen = (x * x) + (y * y);
    return (sqlen < (1e-06 * 1e-06));

   }

   /** As normalise, except that this vector is unaffected and the
   normalised vector is returned as a copy. */
   inline Vector2 normalisedCopy(void) const
   {
    Vector2 ret = *this;
    ret.normalise();
    return ret;
   }

   /** Calculates a reflection vector to the plane with the given normal .
   @remarks NB assumes 'this' is pointing AWAY FROM the plane, invert if it is not.
   */
   inline Vector2 reflect(const Vector2& normal) const
   {
    return Vector2( *this - ( 2 * this->dotProduct(normal) * normal ) );
   }
   /// Check whether this vector contains valid values
   inline bool isNaN() const
   {
    return Math::isNaN(x) || Math::isNaN(y);
   }

   // special points
   static const Vector2 ZERO;
   static const Vector2 UNIT_X;
   static const Vector2 UNIT_Y;
   static const Vector2 NEGATIVE_UNIT_X;
   static const Vector2 NEGATIVE_UNIT_Y;
   static const Vector2 UNIT_SCALE;

   /** Function for writing to a stream.
   */
   inline friend std::ostream& operator <<
    ( std::ostream& o, const Vector2& v )
   {
    o << "Vector2(" << v.x << ", " << v.y <<  ")";
    return o;
   }

  }; // class Vector2
  /** @} */
  /** @} */




  // --------------------------------------------------------------------------------


  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Standard 3-dimensional vector.
  @remarks
  A direction in 3D space represented as distances along the 3
  orthogonal axes (x, y, z). Note that positions, directions and
  scaling factors can be represented by a vector, depending on how
  you interpret the values.
  */
  template<typename o_O> class Vector3T
  {
  public:
   Real x, y, z;

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

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
    if( x < rhs.x && y < rhs.y && z < rhs.z )
     return true;
    return false;
   }

   /** Returns true if the vector's scalar components are all smaller
   that the ones of the vector it is compared against.
   */
   inline bool operator > ( const Vector3& rhs ) const
   {
    if( x > rhs.x && y > rhs.y && z > rhs.z )
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
    if( cmp.x < x ) x = cmp.x;
    if( cmp.y < y ) y = cmp.y;
    if( cmp.z < z ) z = cmp.z;
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
    if( cmp.x > x ) x = cmp.x;
    if( cmp.y > y ) y = cmp.y;
    if( cmp.z > z ) z = cmp.z;
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
    if( perp.squaredLength() < fSquareZero )
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
    const Vector3& up = Vector3::ZERO ) const
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

   /** Gets the angle between 2 vectors.
   @remarks
   Vectors do not have to be unit-length but must represent directions.
   */
   inline Radian angleBetween(const Vector3& dest)
   {
    Real lenProduct = length() * dest.length();

    // Divide by zero check
    if(lenProduct < 1e-6f)
     lenProduct = 1e-6f;

    Real f = dotProduct(dest) / lenProduct;

    f = Math::Clamp(f, (Real)-1.0, (Real)1.0);
    return Math::ACos(f);

   }
   /** Gets the shortest arc quaternion to rotate this vector to the destination
   vector.
   @remarks
   If you call this with a dest vector that is close to the inverse
   of this vector, we will rotate 180 degrees around the 'fallbackAxis'
   (if specified, or a generated axis if not) since in this case
   ANY axis of rotation is valid.
   */
   Quaternion getRotationTo(const Vector3& dest,
    const Vector3& fallbackAxis = Vector3::ZERO) const
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
   inline friend std::ostream& operator <<
    ( std::ostream& o, const Vector3& v )
   {
    o << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
    return o;
   }
  };
  /** @} */
  /** @} */






  // -----------------------------------------------------------------------------------



  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** 4-dimensional homogeneous vector.
  */
  template<typename o_O> class Vector4T
  {
  public:
   Real x, y, z, w;

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   inline Vector4()
   {
   }

   inline Vector4( const Real fX, const Real fY, const Real fZ, const Real fW )
    : x( fX ), y( fY ), z( fZ ), w( fW)
   {
   }

   inline explicit Vector4( const Real afCoordinate[4] )
    : x( afCoordinate[0] ),
    y( afCoordinate[1] ),
    z( afCoordinate[2] ),
    w( afCoordinate[3] )
   {
   }

   inline explicit Vector4( const int afCoordinate[4] )
   {
    x = (Real)afCoordinate[0];
    y = (Real)afCoordinate[1];
    z = (Real)afCoordinate[2];
    w = (Real)afCoordinate[3];
   }

   inline explicit Vector4( Real* const r )
    : x( r[0] ), y( r[1] ), z( r[2] ), w( r[3] )
   {
   }

   inline explicit Vector4( const Real scaler )
    : x( scaler )
    , y( scaler )
    , z( scaler )
    , w( scaler )
   {
   }

   inline explicit Vector4(const Vector3& rhs)
    : x(rhs.x), y(rhs.y), z(rhs.z), w(1.0f)
   {
   }

   /** Exchange the contents of this vector with another. 
   */
   inline void swap(Vector4& other)
   {
    std::swap(x, other.x);
    std::swap(y, other.y);
    std::swap(z, other.z);
    std::swap(w, other.w);
   }

   inline Real operator [] ( const size_t i ) const
   {
    assert( i < 4 );

    return *(&x+i);
   }

   inline Real& operator [] ( const size_t i )
   {
    assert( i < 4 );

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
   inline Vector4& operator = ( const Vector4& rkVector )
   {
    x = rkVector.x;
    y = rkVector.y;
    z = rkVector.z;
    w = rkVector.w;

    return *this;
   }

   inline Vector4& operator = ( const Real fScalar)
   {
    x = fScalar;
    y = fScalar;
    z = fScalar;
    w = fScalar;
    return *this;
   }

   inline bool operator == ( const Vector4& rkVector ) const
   {
    return ( x == rkVector.x &&
     y == rkVector.y &&
     z == rkVector.z &&
     w == rkVector.w );
   }

   inline bool operator != ( const Vector4& rkVector ) const
   {
    return ( x != rkVector.x ||
     y != rkVector.y ||
     z != rkVector.z ||
     w != rkVector.w );
   }

   inline Vector4& operator = (const Vector3& rhs)
   {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    w = 1.0f;
    return *this;
   }

   // arithmetic operations
   inline Vector4 operator + ( const Vector4& rkVector ) const
   {
    return Vector4(
     x + rkVector.x,
     y + rkVector.y,
     z + rkVector.z,
     w + rkVector.w);
   }

   inline Vector4 operator - ( const Vector4& rkVector ) const
   {
    return Vector4(
     x - rkVector.x,
     y - rkVector.y,
     z - rkVector.z,
     w - rkVector.w);
   }

   inline Vector4 operator * ( const Real fScalar ) const
   {
    return Vector4(
     x * fScalar,
     y * fScalar,
     z * fScalar,
     w * fScalar);
   }

   inline Vector4 operator * ( const Vector4& rhs) const
   {
    return Vector4(
     rhs.x * x,
     rhs.y * y,
     rhs.z * z,
     rhs.w * w);
   }

   inline Vector4 operator / ( const Real fScalar ) const
   {
    assert( fScalar != 0.0 );

    Real fInv = 1.0f / fScalar;

    return Vector4(
     x * fInv,
     y * fInv,
     z * fInv,
     w * fInv);
   }

   inline Vector4 operator / ( const Vector4& rhs) const
   {
    return Vector4(
     x / rhs.x,
     y / rhs.y,
     z / rhs.z,
     w / rhs.w);
   }

   inline const Vector4& operator + () const
   {
    return *this;
   }

   inline Vector4 operator - () const
   {
    return Vector4(-x, -y, -z, -w);
   }

   inline friend Vector4 operator * ( const Real fScalar, const Vector4& rkVector )
   {
    return Vector4(
     fScalar * rkVector.x,
     fScalar * rkVector.y,
     fScalar * rkVector.z,
     fScalar * rkVector.w);
   }

   inline friend Vector4 operator / ( const Real fScalar, const Vector4& rkVector )
   {
    return Vector4(
     fScalar / rkVector.x,
     fScalar / rkVector.y,
     fScalar / rkVector.z,
     fScalar / rkVector.w);
   }

   inline friend Vector4 operator + (const Vector4& lhs, const Real rhs)
   {
    return Vector4(
     lhs.x + rhs,
     lhs.y + rhs,
     lhs.z + rhs,
     lhs.w + rhs);
   }

   inline friend Vector4 operator + (const Real lhs, const Vector4& rhs)
   {
    return Vector4(
     lhs + rhs.x,
     lhs + rhs.y,
     lhs + rhs.z,
     lhs + rhs.w);
   }

   inline friend Vector4 operator - (const Vector4& lhs, Real rhs)
   {
    return Vector4(
     lhs.x - rhs,
     lhs.y - rhs,
     lhs.z - rhs,
     lhs.w - rhs);
   }

   inline friend Vector4 operator - (const Real lhs, const Vector4& rhs)
   {
    return Vector4(
     lhs - rhs.x,
     lhs - rhs.y,
     lhs - rhs.z,
     lhs - rhs.w);
   }

   // arithmetic updates
   inline Vector4& operator += ( const Vector4& rkVector )
   {
    x += rkVector.x;
    y += rkVector.y;
    z += rkVector.z;
    w += rkVector.w;

    return *this;
   }

   inline Vector4& operator -= ( const Vector4& rkVector )
   {
    x -= rkVector.x;
    y -= rkVector.y;
    z -= rkVector.z;
    w -= rkVector.w;

    return *this;
   }

   inline Vector4& operator *= ( const Real fScalar )
   {
    x *= fScalar;
    y *= fScalar;
    z *= fScalar;
    w *= fScalar;
    return *this;
   }

   inline Vector4& operator += ( const Real fScalar )
   {
    x += fScalar;
    y += fScalar;
    z += fScalar;
    w += fScalar;
    return *this;
   }

   inline Vector4& operator -= ( const Real fScalar )
   {
    x -= fScalar;
    y -= fScalar;
    z -= fScalar;
    w -= fScalar;
    return *this;
   }

   inline Vector4& operator *= ( const Vector4& rkVector )
   {
    x *= rkVector.x;
    y *= rkVector.y;
    z *= rkVector.z;
    w *= rkVector.w;

    return *this;
   }

   inline Vector4& operator /= ( const Real fScalar )
   {
    assert( fScalar != 0.0 );

    Real fInv = 1.0f / fScalar;

    x *= fInv;
    y *= fInv;
    z *= fInv;
    w *= fInv;

    return *this;
   }

   inline Vector4& operator /= ( const Vector4& rkVector )
   {
    x /= rkVector.x;
    y /= rkVector.y;
    z /= rkVector.z;
    w /= rkVector.w;

    return *this;
   }

   /** Calculates the dot (scalar) product of this vector with another.
   @param
   vec Vector with which to calculate the dot product (together
   with this one).
   @returns
   A float representing the dot product value.
   */
   inline Real dotProduct(const Vector4& vec) const
   {
    return x * vec.x + y * vec.y + z * vec.z + w * vec.w;
   }
   /// Check whether this vector contains valid values
   inline bool isNaN() const
   {
    return Math::isNaN(x) || Math::isNaN(y) || Math::isNaN(z) || Math::isNaN(w);
   }
   /** Function for writing to a stream.
   */
   inline friend std::ostream& operator <<
    ( std::ostream& o, const Vector4& v )
   {
    o << "Vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
    return o;
   }
   // special
   static const Vector4 ZERO;
  }; // class Vector4
  /** @} */
  /** @} */








  // -----------------------------------------------------------------------------------



  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Implementation of a Quaternion, i.e. a rotation around an axis.
  */
  template<typename o_O> class QuaternionT
  {
  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   /// Default constructor, initializes to identity rotation (aka 0°)
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
   friend Quaternion operator* (Real fScalar,
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
   inline friend std::ostream& operator <<
    ( std::ostream& o, const Quaternion& q )
   {
    o << "Quaternion(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
    return o;
   }

  };
  /** @} */
  /** @} */







  // ---------------------------------------------------------------------------------------------




  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Representation of a ray in space, i.e. a line with an origin and direction. */
  template<typename o_O> class RayT
  {

  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

  protected:
   Vector3 mOrigin;
   Vector3 mDirection;
  public:

   Ray():mOrigin(Vector3::ZERO), mDirection(Vector3::UNIT_Z) {}
   Ray(const Vector3& origin, const Vector3& direction)
    :mOrigin(origin), mDirection(direction) {}

   /** Sets the origin of the ray. */
   void setOrigin(const Vector3& origin) {mOrigin = origin;} 
   /** Gets the origin of the ray. */
   const Vector3& getOrigin(void) const {return mOrigin;} 

   /** Sets the direction of the ray. */
   void setDirection(const Vector3& dir) {mDirection = dir;} 
   /** Gets the direction of the ray. */
   const Vector3& getDirection(void) const {return mDirection;} 

   /** Gets the position of a point t units along the ray. */
   Vector3 getPoint(Real t) const { 
    return Vector3(mOrigin + (mDirection * t));
   }

   /** Gets the position of a point t units along the ray. */
   Vector3 operator*(Real t) const { 
    return getPoint(t);
   }

   /** Tests whether this ray intersects the given plane. 
   @returns A pair structure where the first element indicates whether
   an intersection occurs, and if true, the second element will
   indicate the distance along the ray at which it intersects. 
   This can be converted to a point in space by calling getPoint().
   */
   std::pair<bool, Real> intersects(const Plane& p) const
   {
    return Math::intersects(*this, p);
   }

#if BETAJAEN_NEEDS_TO_DO
   /** Tests whether this ray intersects the given plane bounded volume. 
   @returns A pair structure where the first element indicates whether
   an intersection occurs, and if true, the second element will
   indicate the distance along the ray at which it intersects. 
   This can be converted to a point in space by calling getPoint().
   */
   std::pair<bool, Real> intersects(const PlaneBoundedVolume& p) const
   {
    return Math::intersects(*this, p.planes, p.outside == Plane::POSITIVE_SIDE);
   }

#endif
   /** Tests whether this ray intersects the given sphere. 
   @returns A pair structure where the first element indicates whether
   an intersection occurs, and if true, the second element will
   indicate the distance along the ray at which it intersects. 
   This can be converted to a point in space by calling getPoint().
   */
   std::pair<bool, Real> intersects(const Sphere& s) const
   {
    return Math::intersects(*this, s);
   }
   /** Tests whether this ray intersects the given box. 
   @returns A pair structure where the first element indicates whether
   an intersection occurs, and if true, the second element will
   indicate the distance along the ray at which it intersects. 
   This can be converted to a point in space by calling getPoint().
   */
   std::pair<bool, Real> intersects(const AxisAlignedBox& box) const
   {
    return Math::intersects(*this, box);
   }

  };
  /** @} */
  /** @} */





  // ----------------------------------------------------------------------------------------





  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** A sphere primitive, mostly used for bounds checking. 
  @remarks
  A sphere in math texts is normally represented by the function
  x^2 + y^2 + z^2 = r^2 (for sphere's centered on the origin). Ogre stores spheres
  simply as a center point and a radius.
  */
  template<typename o_O> class SphereT
  {

  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;


  protected:
   Real mRadius;
   Vector3 mCenter;
  public:

   /** Standard constructor - creates a unit sphere around the origin.*/
   Sphere() : mRadius(1.0), mCenter(Vector3::ZERO) {}
   /** Constructor allowing arbitrary spheres. 
   @param center The center point of the sphere.
   @param radius The radius of the sphere.
   */
   Sphere(const Vector3& center, Real radius)
    : mRadius(radius), mCenter(center) {}

   /** Returns the radius of the sphere. */
   Real getRadius(void) const { return mRadius; }

   /** Sets the radius of the sphere. */
   void setRadius(Real radius) { mRadius = radius; }

   /** Returns the center point of the sphere. */
   const Vector3& getCenter(void) const { return mCenter; }

   /** Sets the center point of the sphere. */
   void setCenter(const Vector3& center) { mCenter = center; }

   /** Returns whether or not this sphere intersects another sphere. */
   bool intersects(const Sphere& s) const
   {
    return (s.mCenter - mCenter).squaredLength() <=
     Math::Sqr(s.mRadius + mRadius);
   }
   /** Returns whether or not this sphere intersects a box. */
   bool intersects(const AxisAlignedBox& box) const
   {
    return Math::intersects(*this, box);
   }
   /** Returns whether or not this sphere intersects a plane. */
   bool intersects(const Plane& plane) const
   {
    return Math::intersects(*this, plane);
   }
   /** Returns whether or not this sphere intersects a point. */
   bool intersects(const Vector3& v) const
   {
    return ((v - mCenter).squaredLength() <= Math::Sqr(mRadius));
   }


  };
  /** @} */
  /** @} */



  // ---------------------------------------------------------------------------------

  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** Defines a plane in 3D space.
  @remarks
  A plane is defined in 3D space by the equation
  Ax + By + Cz + D = 0
  @par
  This equates to a vector (the normal of the plane, whose x, y
  and z components equate to the coefficients A, B and C
  respectively), and a constant (D) which is the distance along
  the normal you have to go to move the plane back to the origin.
  */
  template<typename o_O> class PlaneT
  {
  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   /** The "positive side" of the plane is the half space to which the
   plane normal points. The "negative side" is the other half
   space. The flag "no side" indicates the plane itself.
   */
   enum Side
   {
    NO_SIDE,
    POSITIVE_SIDE,
    NEGATIVE_SIDE,
    BOTH_SIDE
   };

   Plane ()
   {
    normal = Vector3::ZERO;
    d = 0.0;
   }

   Plane (const Plane& rhs)
   {
    normal = rhs.normal;
    d = rhs.d;
   }

   Plane (const Vector3& rkNormal, Real fConstant)
   {
    normal = rkNormal;
    d = -fConstant;
   }

   Plane (Real a, Real b, Real c, Real _d)
    : normal(a, b, c), d(_d)
   {
   }

   Plane (const Vector3& rkNormal, const Vector3& rkPoint)
   {
    redefine(rkNormal, rkPoint);
   }

   Plane (const Vector3& rkPoint0, const Vector3& rkPoint1,
    const Vector3& rkPoint2)
   {
    redefine(rkPoint0, rkPoint1, rkPoint2);
   }

   Real getDistance (const Vector3& rkPoint) const
   {
    return normal.dotProduct(rkPoint) + d;
   }

   Side getSide (const Vector3& rkPoint) const
   {
    Real fDistance = getDistance(rkPoint);

    if ( fDistance < 0.0 )
     return Plane::NEGATIVE_SIDE;

    if ( fDistance > 0.0 )
     return Plane::POSITIVE_SIDE;

    return Plane::NO_SIDE;
   }


   Side getSide (const AxisAlignedBox& box) const
   {
    if (box.isNull()) 
     return NO_SIDE;
    if (box.isInfinite())
     return BOTH_SIDE;

    return getSide(box.getCenter(), box.getHalfSize());
   }

   Side getSide (const Vector3& centre, const Vector3& halfSize) const
   {
    // Calculate the distance between box centre and the plane
    Real dist = getDistance(centre);

    // Calculate the maximise allows absolute distance for
    // the distance between box centre and plane
    Real maxAbsDist = normal.absDotProduct(halfSize);

    if (dist < -maxAbsDist)
     return Plane::NEGATIVE_SIDE;

    if (dist > +maxAbsDist)
     return Plane::POSITIVE_SIDE;

    return Plane::BOTH_SIDE;
   }

   void redefine(const Vector3& rkPoint0, const Vector3& rkPoint1,
    const Vector3& rkPoint2)
   {
    Vector3 kEdge1 = rkPoint1 - rkPoint0;
    Vector3 kEdge2 = rkPoint2 - rkPoint0;
    normal = kEdge1.crossProduct(kEdge2);
    normal.normalise();
    d = -normal.dotProduct(rkPoint0);
   }

   void redefine(const Vector3& rkNormal, const Vector3& rkPoint)
   {
    normal = rkNormal;
    d = -rkNormal.dotProduct(rkPoint);
   }

   Vector3 projectVector(const Vector3& p) const
   {
    // We know plane normal is unit length, so use simple method
    Matrix3 xform;
    xform[0][0] = 1.0f - normal.x * normal.x;
    xform[0][1] = -normal.x * normal.y;
    xform[0][2] = -normal.x * normal.z;
    xform[1][0] = -normal.y * normal.x;
    xform[1][1] = 1.0f - normal.y * normal.y;
    xform[1][2] = -normal.y * normal.z;
    xform[2][0] = -normal.z * normal.x;
    xform[2][1] = -normal.z * normal.y;
    xform[2][2] = 1.0f - normal.z * normal.z;
    return xform * p;

   }

   Real normalise(void)
   {
    Real fLength = normal.length();

    // Will also work for zero-sized vectors, but will change nothing
    if (fLength > 1e-08f)
    {
     Real fInvLength = 1.0f / fLength;
     normal *= fInvLength;
     d *= fInvLength;
    }

    return fLength;
   }

   std::ostream& operator<< (std::ostream& o)
   {
    o << "Plane(normal=" << p.normal << ", d=" << p.d << ")";
    return o;
   }


   Vector3 normal;
   Real d;

   /// Comparison operator
   bool operator==(const Plane& rhs) const
   {
    return (rhs.d == d && rhs.normal == normal);
   }
   bool operator!=(const Plane& rhs) const
   {
    return (rhs.d != d && rhs.normal != normal);
   }

  };

  typedef vector< PlaneT<> >::type PlaneList;
  /** @} */
  /** @} */






  // ------------------------------------------------------------------------------------



  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */

  /** A 3D box aligned with the x/y/z axes.
  @remarks
  This class represents a simple box which is aligned with the
  axes. Internally it only stores 2 points as the extremeties of
  the box, one which is the minima of all 3 axes, and the other
  which is the maxima of all 3 axes. This class is typically used
  for an axis-aligned bounding box (AABB) for collision and
  visibility determination.
  */
  template<typename o_O> class AxisAlignedBoxT
  {
  public:


   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

   enum Extent
   {
    EXTENT_NULL,
    EXTENT_FINITE,
    EXTENT_INFINITE
   };
  protected:

   Vector3 mMinimum;
   Vector3 mMaximum;
   Extent mExtent;
   mutable Vector3* mpCorners;

  public:
   /*
   1-----2
   /|    /|
   / |   / |
   5-----4  |
   |  0--|--3
   | /   | /
   |/    |/
   6-----7
   */
   typedef enum {
    FAR_LEFT_BOTTOM = 0,
    FAR_LEFT_TOP = 1,
    FAR_RIGHT_TOP = 2,
    FAR_RIGHT_BOTTOM = 3,
    NEAR_RIGHT_BOTTOM = 7,
    NEAR_LEFT_BOTTOM = 6,
    NEAR_LEFT_TOP = 5,
    NEAR_RIGHT_TOP = 4
   } CornerEnum;
   inline AxisAlignedBox() : mMinimum(Vector3::ZERO), mMaximum(Vector3::UNIT_SCALE), mpCorners(0)
   {
    // Default to a null box 
    setMinimum( -0.5, -0.5, -0.5 );
    setMaximum( 0.5, 0.5, 0.5 );
    mExtent = EXTENT_NULL;
   }
   inline AxisAlignedBox(Extent e) : mMinimum(Vector3::ZERO), mMaximum(Vector3::UNIT_SCALE), mpCorners(0)
   {
    setMinimum( -0.5, -0.5, -0.5 );
    setMaximum( 0.5, 0.5, 0.5 );
    mExtent = e;
   }

   inline AxisAlignedBox(const AxisAlignedBox & rkBox) : mMinimum(Vector3::ZERO), mMaximum(Vector3::UNIT_SCALE), mpCorners(0)

   {
    if (rkBox.isNull())
     setNull();
    else if (rkBox.isInfinite())
     setInfinite();
    else
     setExtents( rkBox.mMinimum, rkBox.mMaximum );
   }

   inline AxisAlignedBox( const Vector3& min, const Vector3& max ) : mMinimum(Vector3::ZERO), mMaximum(Vector3::UNIT_SCALE), mpCorners(0)
   {
    setExtents( min, max );
   }

   inline AxisAlignedBox(
    Real mx, Real my, Real mz,
    Real Mx, Real My, Real Mz ) : mMinimum(Vector3::ZERO), mMaximum(Vector3::UNIT_SCALE), mpCorners(0)
   {
    setExtents( mx, my, mz, Mx, My, Mz );
   }

   AxisAlignedBox& operator=(const AxisAlignedBox& rhs)
   {
    // Specifically override to avoid copying mpCorners
    if (rhs.isNull())
     setNull();
    else if (rhs.isInfinite())
     setInfinite();
    else
     setExtents(rhs.mMinimum, rhs.mMaximum);

    return *this;
   }

   ~AxisAlignedBoxT<o_O>()
   {
    if (mpCorners)
     OGREMATH_FREE(mpCorners, MEMCATEGORY_SCENE_CONTROL);
   }


   /** Gets the minimum corner of the box.
   */
   inline const Vector3& getMinimum(void) const
   { 
    return mMinimum; 
   }

   /** Gets a modifiable version of the minimum
   corner of the box.
   */
   inline Vector3& getMinimum(void)
   { 
    return mMinimum; 
   }

   /** Gets the maximum corner of the box.
   */
   inline const Vector3& getMaximum(void) const
   { 
    return mMaximum;
   }

   /** Gets a modifiable version of the maximum
   corner of the box.
   */
   inline Vector3& getMaximum(void)
   { 
    return mMaximum;
   }


   /** Sets the minimum corner of the box.
   */
   inline void setMinimum( const Vector3& vec )
   {
    mExtent = EXTENT_FINITE;
    mMinimum = vec;
   }

   inline void setMinimum( Real x, Real y, Real z )
   {
    mExtent = EXTENT_FINITE;
    mMinimum.x = x;
    mMinimum.y = y;
    mMinimum.z = z;
   }

   /** Changes one of the components of the minimum corner of the box
   used to resize only one dimension of the box
   */
   inline void setMinimumX(Real x)
   {
    mMinimum.x = x;
   }

   inline void setMinimumY(Real y)
   {
    mMinimum.y = y;
   }

   inline void setMinimumZ(Real z)
   {
    mMinimum.z = z;
   }

   /** Sets the maximum corner of the box.
   */
   inline void setMaximum( const Vector3& vec )
   {
    mExtent = EXTENT_FINITE;
    mMaximum = vec;
   }

   inline void setMaximum( Real x, Real y, Real z )
   {
    mExtent = EXTENT_FINITE;
    mMaximum.x = x;
    mMaximum.y = y;
    mMaximum.z = z;
   }

   /** Changes one of the components of the maximum corner of the box
   used to resize only one dimension of the box
   */
   inline void setMaximumX( Real x )
   {
    mMaximum.x = x;
   }

   inline void setMaximumY( Real y )
   {
    mMaximum.y = y;
   }

   inline void setMaximumZ( Real z )
   {
    mMaximum.z = z;
   }

   /** Sets both minimum and maximum extents at once.
   */
   inline void setExtents( const Vector3& min, const Vector3& max )
   {
    assert( (min.x <= max.x && min.y <= max.y && min.z <= max.z) &&
     "The minimum corner of the box must be less than or equal to maximum corner" );

    mExtent = EXTENT_FINITE;
    mMinimum = min;
    mMaximum = max;
   }

   inline void setExtents(
    Real mx, Real my, Real mz,
    Real Mx, Real My, Real Mz )
   {
    assert( (mx <= Mx && my <= My && mz <= Mz) &&
     "The minimum corner of the box must be less than or equal to maximum corner" );

    mExtent = EXTENT_FINITE;

    mMinimum.x = mx;
    mMinimum.y = my;
    mMinimum.z = mz;

    mMaximum.x = Mx;
    mMaximum.y = My;
    mMaximum.z = Mz;

   }

   /** Returns a pointer to an array of 8 corner points, useful for
   collision vs. non-aligned objects.
   @remarks
   If the order of these corners is important, they are as
   follows: The 4 points of the minimum Z face (note that
   because Ogre uses right-handed coordinates, the minimum Z is
   at the 'back' of the box) starting with the minimum point of
   all, then anticlockwise around this face (if you are looking
   onto the face from outside the box). Then the 4 points of the
   maximum Z face, starting with maximum point of all, then
   anticlockwise around this face (looking onto the face from
   outside the box). Like this:
   <pre>
   1-----2
   /|    /|
   / |   / |
   5-----4  |
   |  0--|--3
   | /   | /
   |/    |/
   6-----7
   </pre>
   @remarks as this implementation uses a static member, make sure to use your own copy !
   */
   inline const Vector3* getAllCorners(void) const
   {
    assert( (mExtent == EXTENT_FINITE) && "Can't get corners of a null or infinite AAB" );

    // The order of these items is, using right-handed co-ordinates:
    // Minimum Z face, starting with Min(all), then anticlockwise
    //   around face (looking onto the face)
    // Maximum Z face, starting with Max(all), then anticlockwise
    //   around face (looking onto the face)
    // Only for optimization/compatibility.
    if (!mpCorners)
     mpCorners = OGRE_ALLOC_T(Vector3, 8, MEMCATEGORY_SCENE_CONTROL);

    mpCorners[0] = mMinimum;
    mpCorners[1].x = mMinimum.x; mpCorners[1].y = mMaximum.y; mpCorners[1].z = mMinimum.z;
    mpCorners[2].x = mMaximum.x; mpCorners[2].y = mMaximum.y; mpCorners[2].z = mMinimum.z;
    mpCorners[3].x = mMaximum.x; mpCorners[3].y = mMinimum.y; mpCorners[3].z = mMinimum.z;            

    mpCorners[4] = mMaximum;
    mpCorners[5].x = mMinimum.x; mpCorners[5].y = mMaximum.y; mpCorners[5].z = mMaximum.z;
    mpCorners[6].x = mMinimum.x; mpCorners[6].y = mMinimum.y; mpCorners[6].z = mMaximum.z;
    mpCorners[7].x = mMaximum.x; mpCorners[7].y = mMinimum.y; mpCorners[7].z = mMaximum.z;

    return mpCorners;
   }

   /** gets the position of one of the corners
   */
   Vector3 getCorner(CornerEnum cornerToGet) const
   {
    switch(cornerToGet)
    {
    case FAR_LEFT_BOTTOM:
     return mMinimum;
    case FAR_LEFT_TOP:
     return Vector3(mMinimum.x, mMaximum.y, mMinimum.z);
    case FAR_RIGHT_TOP:
     return Vector3(mMaximum.x, mMaximum.y, mMinimum.z);
    case FAR_RIGHT_BOTTOM:
     return Vector3(mMaximum.x, mMinimum.y, mMinimum.z);
    case NEAR_RIGHT_BOTTOM:
     return Vector3(mMaximum.x, mMinimum.y, mMaximum.z);
    case NEAR_LEFT_BOTTOM:
     return Vector3(mMinimum.x, mMinimum.y, mMaximum.z);
    case NEAR_LEFT_TOP:
     return Vector3(mMinimum.x, mMaximum.y, mMaximum.z);
    case NEAR_RIGHT_TOP:
     return mMaximum;
    default:
     return Vector3();
    }
   }

   friend std::ostream& operator<<( std::ostream& o, const AxisAlignedBox aab )
   {
    switch (aab.mExtent)
    {
    case EXTENT_NULL:
     o << "AxisAlignedBox(null)";
     return o;

    case EXTENT_FINITE:
     o << "AxisAlignedBox(min=" << aab.mMinimum << ", max=" << aab.mMaximum << ")";
     return o;

    case EXTENT_INFINITE:
     o << "AxisAlignedBox(infinite)";
     return o;

    default: // shut up compiler
     assert( false && "Never reached" );
     return o;
    }
   }

   /** Merges the passed in box into the current box. The result is the
   box which encompasses both.
   */
   void merge( const AxisAlignedBox& rhs )
   {
    // Do nothing if rhs null, or this is infinite
    if ((rhs.mExtent == EXTENT_NULL) || (mExtent == EXTENT_INFINITE))
    {
     return;
    }
    // Otherwise if rhs is infinite, make this infinite, too
    else if (rhs.mExtent == EXTENT_INFINITE)
    {
     mExtent = EXTENT_INFINITE;
    }
    // Otherwise if current null, just take rhs
    else if (mExtent == EXTENT_NULL)
    {
     setExtents(rhs.mMinimum, rhs.mMaximum);
    }
    // Otherwise merge
    else
    {
     Vector3 min = mMinimum;
     Vector3 max = mMaximum;
     max.makeCeil(rhs.mMaximum);
     min.makeFloor(rhs.mMinimum);

     setExtents(min, max);
    }

   }

   /** Extends the box to encompass the specified point (if needed).
   */
   inline void merge( const Vector3& point )
   {
    switch (mExtent)
    {
    case EXTENT_NULL: // if null, use this point
     setExtents(point, point);
     return;

    case EXTENT_FINITE:
     mMaximum.makeCeil(point);
     mMinimum.makeFloor(point);
     return;

    case EXTENT_INFINITE: // if infinite, makes no difference
     return;
    }

    assert( false && "Never reached" );
   }

   /** Transforms the box according to the matrix supplied.
   @remarks
   By calling this method you get the axis-aligned box which
   surrounds the transformed version of this box. Therefore each
   corner of the box is transformed by the matrix, then the
   extents are mapped back onto the axes to produce another
   AABB. Useful when you have a local AABB for an object which
   is then transformed.
   */
   inline void transform( const Matrix4& matrix )
   {
    // Do nothing if current null or infinite
    if( mExtent != EXTENT_FINITE )
     return;

    Vector3 oldMin, oldMax, currentCorner;

    // Getting the old values so that we can use the existing merge method.
    oldMin = mMinimum;
    oldMax = mMaximum;

    // reset
    setNull();

    // We sequentially compute the corners in the following order :
    // 0, 6, 5, 1, 2, 4 ,7 , 3
    // This sequence allows us to only change one member at a time to get at all corners.

    // For each one, we transform it using the matrix
    // Which gives the resulting point and merge the resulting point.

    // First corner 
    // min min min
    currentCorner = oldMin;
    merge( matrix * currentCorner );

    // min,min,max
    currentCorner.z = oldMax.z;
    merge( matrix * currentCorner );

    // min max max
    currentCorner.y = oldMax.y;
    merge( matrix * currentCorner );

    // min max min
    currentCorner.z = oldMin.z;
    merge( matrix * currentCorner );

    // max max min
    currentCorner.x = oldMax.x;
    merge( matrix * currentCorner );

    // max max max
    currentCorner.z = oldMax.z;
    merge( matrix * currentCorner );

    // max min max
    currentCorner.y = oldMin.y;
    merge( matrix * currentCorner );

    // max min min
    currentCorner.z = oldMin.z;
    merge( matrix * currentCorner ); 
   }

   /** Transforms the box according to the affine matrix supplied.
   @remarks
   By calling this method you get the axis-aligned box which
   surrounds the transformed version of this box. Therefore each
   corner of the box is transformed by the matrix, then the
   extents are mapped back onto the axes to produce another
   AABB. Useful when you have a local AABB for an object which
   is then transformed.
   @note
   The matrix must be an affine matrix. @see Matrix4::isAffine.
   */
   void transformAffine(const Matrix4& m)
   {
    assert(m.isAffine());

    // Do nothing if current null or infinite
    if ( mExtent != EXTENT_FINITE )
     return;

    Vector3 centre = getCenter();
    Vector3 halfSize = getHalfSize();

    Vector3 newCentre = m.transformAffine(centre);
    Vector3 newHalfSize(
     Math::Abs(m[0][0]) * halfSize.x + Math::Abs(m[0][1]) * halfSize.y + Math::Abs(m[0][2]) * halfSize.z, 
     Math::Abs(m[1][0]) * halfSize.x + Math::Abs(m[1][1]) * halfSize.y + Math::Abs(m[1][2]) * halfSize.z,
     Math::Abs(m[2][0]) * halfSize.x + Math::Abs(m[2][1]) * halfSize.y + Math::Abs(m[2][2]) * halfSize.z);

    setExtents(newCentre - newHalfSize, newCentre + newHalfSize);
   }

   /** Sets the box to a 'null' value i.e. not a box.
   */
   inline void setNull()
   {
    mExtent = EXTENT_NULL;
   }

   /** Returns true if the box is null i.e. empty.
   */
   inline bool isNull(void) const
   {
    return (mExtent == EXTENT_NULL);
   }

   /** Returns true if the box is finite.
   */
   bool isFinite(void) const
   {
    return (mExtent == EXTENT_FINITE);
   }

   /** Sets the box to 'infinite'
   */
   inline void setInfinite()
   {
    mExtent = EXTENT_INFINITE;
   }

   /** Returns true if the box is infinite.
   */
   bool isInfinite(void) const
   {
    return (mExtent == EXTENT_INFINITE);
   }

   /** Returns whether or not this box intersects another. */
   inline bool intersects(const AxisAlignedBox& b2) const
   {
    // Early-fail for nulls
    if (this->isNull() || b2.isNull())
     return false;

    // Early-success for infinites
    if (this->isInfinite() || b2.isInfinite())
     return true;

    // Use up to 6 separating planes
    if (mMaximum.x < b2.mMinimum.x)
     return false;
    if (mMaximum.y < b2.mMinimum.y)
     return false;
    if (mMaximum.z < b2.mMinimum.z)
     return false;

    if (mMinimum.x > b2.mMaximum.x)
     return false;
    if (mMinimum.y > b2.mMaximum.y)
     return false;
    if (mMinimum.z > b2.mMaximum.z)
     return false;

    // otherwise, must be intersecting
    return true;

   }

   /// Calculate the area of intersection of this box and another
   inline AxisAlignedBox intersection(const AxisAlignedBox& b2) const
   {
    if (this->isNull() || b2.isNull())
    {
     return AxisAlignedBox();
    }
    else if (this->isInfinite())
    {
     return b2;
    }
    else if (b2.isInfinite())
    {
     return *this;
    }

    Vector3 intMin = mMinimum;
    Vector3 intMax = mMaximum;

    intMin.makeCeil(b2.getMinimum());
    intMax.makeFloor(b2.getMaximum());

    // Check intersection isn't null
    if (intMin.x < intMax.x &&
     intMin.y < intMax.y &&
     intMin.z < intMax.z)
    {
     return AxisAlignedBox(intMin, intMax);
    }

    return AxisAlignedBox();
   }

   /// Calculate the volume of this box
   Real volume(void) const
   {
    switch (mExtent)
    {
    case EXTENT_NULL:
     return 0.0f;

    case EXTENT_FINITE:
     {
      Vector3 diff = mMaximum - mMinimum;
      return diff.x * diff.y * diff.z;
     }

    case EXTENT_INFINITE:
     return Math::POS_INFINITY;

    default: // shut up compiler
     assert( false && "Never reached" );
     return 0.0f;
    }
   }

   /** Scales the AABB by the vector given. */
   inline void scale(const Vector3& s)
   {
    // Do nothing if current null or infinite
    if (mExtent != EXTENT_FINITE)
     return;

    // NB assumes centered on origin
    Vector3 min = mMinimum * s;
    Vector3 max = mMaximum * s;
    setExtents(min, max);
   }

   /** Tests whether this box intersects a sphere. */
   bool intersects(const Sphere& s) const
   {
    return Math::intersects(s, *this); 
   }
   /** Tests whether this box intersects a plane. */
   bool intersects(const Plane& p) const
   {
    return Math::intersects(p, *this);
   }
   /** Tests whether the vector point is within this box. */
   bool intersects(const Vector3& v) const
   {
    switch (mExtent)
    {
    case EXTENT_NULL:
     return false;

    case EXTENT_FINITE:
     return(v.x >= mMinimum.x  &&  v.x <= mMaximum.x  && 
      v.y >= mMinimum.y  &&  v.y <= mMaximum.y  && 
      v.z >= mMinimum.z  &&  v.z <= mMaximum.z);

    case EXTENT_INFINITE:
     return true;

    default: // shut up compiler
     assert( false && "Never reached" );
     return false;
    }
   }
   /// Gets the centre of the box
   Vector3 getCenter(void) const
   {
    assert( (mExtent == EXTENT_FINITE) && "Can't get center of a null or infinite AAB" );

    return Vector3(
     (mMaximum.x + mMinimum.x) * 0.5f,
     (mMaximum.y + mMinimum.y) * 0.5f,
     (mMaximum.z + mMinimum.z) * 0.5f);
   }
   /// Gets the size of the box
   Vector3 getSize(void) const
   {
    switch (mExtent)
    {
    case EXTENT_NULL:
     return Vector3::ZERO;

    case EXTENT_FINITE:
     return mMaximum - mMinimum;

    case EXTENT_INFINITE:
     return Vector3(
      Math::POS_INFINITY,
      Math::POS_INFINITY,
      Math::POS_INFINITY);

    default: // shut up compiler
     assert( false && "Never reached" );
     return Vector3::ZERO;
    }
   }
   /// Gets the half-size of the box
   Vector3 getHalfSize(void) const
   {
    switch (mExtent)
    {
    case EXTENT_NULL:
     return Vector3::ZERO;

    case EXTENT_FINITE:
     return (mMaximum - mMinimum) * 0.5;

    case EXTENT_INFINITE:
     return Vector3(
      Math::POS_INFINITY,
      Math::POS_INFINITY,
      Math::POS_INFINITY);

    default: // shut up compiler
     assert( false && "Never reached" );
     return Vector3::ZERO;
    }
   }

   /** Tests whether the given point contained by this box.
   */
   bool contains(const Vector3& v) const
   {
    if (isNull())
     return false;
    if (isInfinite())
     return true;

    return mMinimum.x <= v.x && v.x <= mMaximum.x &&
     mMinimum.y <= v.y && v.y <= mMaximum.y &&
     mMinimum.z <= v.z && v.z <= mMaximum.z;
   }

   /** Tests whether another box contained by this box.
   */
   bool contains(const AxisAlignedBox& other) const
   {
    if (other.isNull() || this->isInfinite())
     return true;

    if (this->isNull() || other.isInfinite())
     return false;

    return this->mMinimum.x <= other.mMinimum.x &&
     this->mMinimum.y <= other.mMinimum.y &&
     this->mMinimum.z <= other.mMinimum.z &&
     other.mMaximum.x <= this->mMaximum.x &&
     other.mMaximum.y <= this->mMaximum.y &&
     other.mMaximum.z <= this->mMaximum.z;
   }

   /** Tests 2 boxes for equality.
   */
   bool operator== (const AxisAlignedBox& rhs) const
   {
    if (this->mExtent != rhs.mExtent)
     return false;

    if (!this->isFinite())
     return true;

    return this->mMinimum == rhs.mMinimum &&
     this->mMaximum == rhs.mMaximum;
   }

   /** Tests 2 boxes for inequality.
   */
   bool operator!= (const AxisAlignedBox& rhs) const
   {
    return !(*this == rhs);
   }

   // special values
   static const AxisAlignedBox BOX_NULL;
   static const AxisAlignedBox BOX_INFINITE;


  };

  /** @} */
  /** @} */






  // ---------------------------------------------------------------------------------


  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
  /** A 3x3 matrix which can represent rotations around axes.
  @note
  <b>All the code is adapted from the Wild Magic 0.2 Matrix
  library (http://www.geometrictools.com/).</b>
  @par
  The coordinate system is assumed to be <b>right-handed</b>.


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

  */
  template<typename o_O> class Matrix3T
  {



  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;

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
   Vector3 GetColumn (size_t iCol) const
   {
    assert( 0 <= iCol && iCol < 3 );
    return Vector3(m[0][iCol],m[1][iCol],
     m[2][iCol]);
   }

   void SetColumn(size_t iCol, const Vector3& vec)
   {
    assert( 0 <= iCol && iCol < 3 );
    m[0][iCol] = vec.x;
    m[1][iCol] = vec.y;
    m[2][iCol] = vec.z;
   }

   void FromAxes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis)
   {
    SetColumn(0,xAxis);
    SetColumn(1,yAxis);
    SetColumn(2,zAxis);
   }

   // assignment and comparison
   inline Matrix3& operator= (const Matrix3& rkMatrix)
   {
    memcpy(m,rkMatrix.m,9*sizeof(Real));
    return *this;
   }

   bool operator== (const Matrix3& rkMatrix) const
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

   inline bool operator!= (const Matrix3& rkMatrix) const
   {
    return !operator==(rkMatrix);
   }

   // arithmetic operations
   Matrix3 operator+ (const Matrix3& rkMatrix) const
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

   Matrix3 operator- (const Matrix3& rkMatrix) const
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

   Matrix3 operator* (const Matrix3& rkMatrix) const
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

   Matrix3 operator- () const
   {
    Matrix3 kNeg;
    for (size_t iRow = 0; iRow < 3; iRow++)
    {
     for (size_t iCol = 0; iCol < 3; iCol++)
      kNeg[iRow][iCol] = -m[iRow][iCol];
    }
    return kNeg;
   }


   // matrix * vector [3x3 * 3x1 = 3x1]
   Vector3 operator* (const Vector3& rkVector) const
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


   // vector * matrix [1x3 * 3x3 = 1x3]
   friend Vector3 operator* (const Vector3& rkVector,
    const Matrix3& rkMatrix)
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


   // matrix * scalar
   Matrix3 operator* (Real fScalar) const
   {
    Matrix3 kProd;
    for (size_t iRow = 0; iRow < 3; iRow++)
    {
     for (size_t iCol = 0; iCol < 3; iCol++)
      kProd[iRow][iCol] = fScalar*m[iRow][iCol];
    }
    return kProd;
   }


   // scalar * matrix
   friend Matrix3 operator* (Real fScalar, const Matrix3& rkMatrix)
   {
    Matrix3 kProd;
    for (size_t iRow = 0; iRow < 3; iRow++)
    {
     for (size_t iCol = 0; iCol < 3; iCol++)
      kProd[iRow][iCol] = fScalar*rkMatrix.m[iRow][iCol];
    }
    return kProd;
   }


   // utilities
   Matrix3 Transpose () const
   {
    Matrix3 kTranspose;
    for (size_t iRow = 0; iRow < 3; iRow++)
    {
     for (size_t iCol = 0; iCol < 3; iCol++)
      kTranspose[iRow][iCol] = m[iCol][iRow];
    }
    return kTranspose;
   }

   bool Inverse (Matrix3& rkInverse, Real fTolerance = 1e-06) const
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

   Matrix3 Inverse (Real fTolerance = 1e-06) const
   {
    Matrix3 kInverse = Matrix3::ZERO;
    Inverse(kInverse,fTolerance);
    return kInverse;
   }

   Real Determinant () const
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


   // singular value decomposition
   void SingularValueDecomposition (Matrix3& rkL, Vector3& rkS,
    Matrix3& rkR) const
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

   void SingularValueComposition (const Matrix3& rkL,
    const Vector3& rkS, const Matrix3& rkR)
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


   // Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
   void Orthonormalize ()
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


   // orthogonal Q, diagonal D, upper triangular U stored as (u01,u02,u12)
   void QDUDecomposition (Matrix3& rkQ, Vector3& rkD,
    Vector3& rkU) const
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


   Real SpectralNorm () const
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


   // matrix must be orthonormal
   void ToAxisAngle (Vector3& rkAxis, Radian& rfAngle) const
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


   inline void ToAxisAngle (Vector3& rkAxis, Degree& rfAngle) const {
    Radian r;
    ToAxisAngle ( rkAxis, r );
    rfAngle = r;
   }

   void FromAxisAngle (const Vector3& rkAxis, const Radian& fRadians)
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


   // The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
   // where yaw is rotation about the Up vector, pitch is rotation about the
   // Right axis, and roll is rotation about the Direction axis.
   bool ToEulerAnglesXYZ (Radian& rfYAngle, Radian& rfPAngle,
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



   bool ToEulerAnglesXZY (Radian& rfYAngle, Radian& rfPAngle,
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


   bool ToEulerAnglesYXZ (Radian& rfYAngle, Radian& rfPAngle,
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


   bool ToEulerAnglesYZX (Radian& rfYAngle, Radian& rfPAngle,
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


   bool ToEulerAnglesZXY (Radian& rfYAngle, Radian& rfPAngle,
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


   bool ToEulerAnglesZYX (Radian& rfYAngle, Radian& rfPAngle,
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



   void FromEulerAnglesXYZ (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   void FromEulerAnglesXZY (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   void FromEulerAnglesYXZ (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   void FromEulerAnglesYZX (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   void FromEulerAnglesZXY (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   void FromEulerAnglesZYX (const Radian& fYAngle, const Radian& fPAngle, const Radian& fRAngle)
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


   // eigensolver, matrix must be symmetric
   void EigenSolveSymmetric (Real afEigenvalue[3],
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


   static void TensorProduct (const Vector3& rkU, const Vector3& rkV,
    Matrix3& rkProduct)
   {
    for (size_t iRow = 0; iRow < 3; iRow++)
    {
     for (size_t iCol = 0; iCol < 3; iCol++)
      rkProduct[iRow][iCol] = rkU[iRow]*rkV[iCol];
    }
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

   /** Function for writing to a stream.
   */
   inline friend std::ostream& operator <<
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
   void Tridiagonal (Real afDiag[3], Real afSubDiag[3])
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


   bool QLAlgorithm (Real afDiag[3], Real afSubDiag[3])
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


   // support for singular value decomposition
   static const Real ms_fSvdEpsilon;
   static const unsigned int ms_iSvdMaxIterations;
   static void Bidiagonalize (Matrix3& kA, Matrix3& kL,
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


   static void GolubKahanStep (Matrix3& kA, Matrix3& kL,
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


   // support for spectral norm
   static Real MaxCubicRoot (Real afCoeff[3])
   {
    // Spectral norm is for A^T*A, so characteristic polynomial
    // P(x) = c[0]+c[1]*x+c[2]*x^2+x^3 has three positive real roots.
    // This yields the assertions c[0] < 0 and c[2]*c[2] >= 3*c[1].

    // quick out for uniform scale (triple root)
    const Real fOneThird = 1.0/3.0;
    const Real fEpsilon = 1e-06;
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


   Real m[3][3];

   // for faster access
   friend class Matrix4;
  };
  /** @} */
  /** @} */




  // ---------------------------------------------------------------------------------



  /** \addtogroup Core
  *  @{
  */
  /** \addtogroup Math
  *  @{
  */
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
  template<typename o_O> class Matrix4T
  {
  protected:
   /// The matrix entries, indexed by [row][col].
   union {
    Real m[4][4];
    Real _m[16];
   };
  public:


   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector4T<o_O>         Vector4;
   typedef Private::QuaternionT<o_O>      Quaternion;
   typedef Private::RayT<o_O>             Ray;
   typedef Private::SphereT<o_O>          Sphere;
   typedef Private::PlaneT<o_O>           Plane;
   typedef Private::AxisAlignedBoxT<o_O>  AxisAlignedBox;
   typedef Private::Matrix3T<o_O>         Matrix3;
   typedef Private::Matrix4T<o_O>         Matrix4;


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
   inline Vector4 operator * (const Vector4& v) const
   {
    return Vector4(
     m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w, 
     m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w,
     m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w,
     m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w
     );
   }
   inline Plane operator * (const Plane& p) const
   {
    Plane ret;
    Matrix4 invTrans = inverse().transpose();
    Vector4 v4( p.normal.x, p.normal.y, p.normal.z, p.d );
    v4 = invTrans * v4;
    ret.normal.x = v4.x; 
    ret.normal.y = v4.y; 
    ret.normal.z = v4.z;
    ret.d = v4.w / ret.normal.normalise();

    return ret;
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
    if( 
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
    if( 
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
    m[0][0] = mat3.m[0][0]; m[0][1] = mat3.m[0][1]; m[0][2] = mat3.m[0][2];
    m[1][0] = mat3.m[1][0]; m[1][1] = mat3.m[1][1]; m[1][2] = mat3.m[1][2];
    m[2][0] = mat3.m[2][0]; m[2][1] = mat3.m[2][1]; m[2][2] = mat3.m[2][2];
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
    m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = v.x;
    m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = v.y;
    m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = v.z;
    m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0;
   }

   inline void makeTrans( Real tx, Real ty, Real tz )
   {
    m[0][0] = 1.0; m[0][1] = 0.0; m[0][2] = 0.0; m[0][3] = tx;
    m[1][0] = 0.0; m[1][1] = 1.0; m[1][2] = 0.0; m[1][3] = ty;
    m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 1.0; m[2][3] = tz;
    m[3][0] = 0.0; m[3][1] = 0.0; m[3][2] = 0.0; m[3][3] = 1.0;
   }

   /** Gets a translation matrix.
   */
   inline static Matrix4 getTrans( const Vector3& v )
   {
    Matrix4 r;

    r.m[0][0] = 1.0; r.m[0][1] = 0.0; r.m[0][2] = 0.0; r.m[0][3] = v.x;
    r.m[1][0] = 0.0; r.m[1][1] = 1.0; r.m[1][2] = 0.0; r.m[1][3] = v.y;
    r.m[2][0] = 0.0; r.m[2][1] = 0.0; r.m[2][2] = 1.0; r.m[2][3] = v.z;
    r.m[3][0] = 0.0; r.m[3][1] = 0.0; r.m[3][2] = 0.0; r.m[3][3] = 1.0;

    return r;
   }

   /** Gets a translation matrix - variation for not using a vector.
   */
   inline static Matrix4 getTrans( Real t_x, Real t_y, Real t_z )
   {
    Matrix4 r;

    r.m[0][0] = 1.0; r.m[0][1] = 0.0; r.m[0][2] = 0.0; r.m[0][3] = t_x;
    r.m[1][0] = 0.0; r.m[1][1] = 1.0; r.m[1][2] = 0.0; r.m[1][3] = t_y;
    r.m[2][0] = 0.0; r.m[2][1] = 0.0; r.m[2][2] = 1.0; r.m[2][3] = t_z;
    r.m[3][0] = 0.0; r.m[3][1] = 0.0; r.m[3][2] = 0.0; r.m[3][3] = 1.0;

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
    r.m[0][0] = v.x; r.m[0][1] = 0.0; r.m[0][2] = 0.0; r.m[0][3] = 0.0;
    r.m[1][0] = 0.0; r.m[1][1] = v.y; r.m[1][2] = 0.0; r.m[1][3] = 0.0;
    r.m[2][0] = 0.0; r.m[2][1] = 0.0; r.m[2][2] = v.z; r.m[2][3] = 0.0;
    r.m[3][0] = 0.0; r.m[3][1] = 0.0; r.m[3][2] = 0.0; r.m[3][3] = 1.0;

    return r;
   }

   /** Gets a scale matrix - variation for not using a vector.
   */
   inline static Matrix4 getScale( Real s_x, Real s_y, Real s_z )
   {
    Matrix4 r;
    r.m[0][0] = s_x; r.m[0][1] = 0.0; r.m[0][2] = 0.0; r.m[0][3] = 0.0;
    r.m[1][0] = 0.0; r.m[1][1] = s_y; r.m[1][2] = 0.0; r.m[1][3] = 0.0;
    r.m[2][0] = 0.0; r.m[2][1] = 0.0; r.m[2][2] = s_z; r.m[2][3] = 0.0;
    r.m[3][0] = 0.0; r.m[3][1] = 0.0; r.m[3][2] = 0.0; r.m[3][3] = 1.0;

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
   inline friend std::ostream& operator <<
    ( std::ostream& o, const Matrix4& mat )
   {
    o << "Matrix4(";
    for (size_t i = 0; i < 4; ++i)
    {
     o << " row" << (unsigned)i << "{";
     for(size_t j = 0; j < 4; ++j)
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

   /** 4-D Vector transformation specially for an affine matrix.
   @note
   The matrix must be an affine matrix. @see Matrix4::isAffine.
   */
   inline Vector4 transformAffine(const Vector4& v) const
   {
    assert(isAffine());

    return Vector4(
     m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w, 
     m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w,
     m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w,
     v.w);
   }
  };

  /* Removed from Vector4 and made a non-member here because otherwise
  OgreMatrix4.h and OgreVector4.h have to try to include and inline each 
  other, which frankly doesn't work ;)
  */
  inline Vector4T<> operator * (const Vector4T<>& v, const Matrix4T<>& mat)
  {
   return Vector4T<>(
    v.x*mat[0][0] + v.y*mat[1][0] + v.z*mat[2][0] + v.w*mat[3][0],
    v.x*mat[0][1] + v.y*mat[1][1] + v.z*mat[2][1] + v.w*mat[3][1],
    v.x*mat[0][2] + v.y*mat[1][2] + v.z*mat[2][2] + v.w*mat[3][2],
    v.x*mat[0][3] + v.y*mat[1][3] + v.z*mat[2][3] + v.w*mat[3][3]
   );
  }
  /** @} */
  /** @} */



 } // namespace OgreMath::Private


 typedef Private::MathT<>            Math;
 typedef Private::DegreeT<>          Degree;
 typedef Private::RadianT<>          Radian;
 typedef Private::Vector2T<>         Vector2;
 typedef Private::Vector3T<>         Vector3;
 typedef Private::Vector4T<>         Vector4;
 typedef Private::QuaternionT<>      Quaternion;
 typedef Private::RayT<>             Ray;
 typedef Private::SphereT<>          Sphere;
 typedef Private::PlaneT<>           Plane;
 typedef Private::AxisAlignedBoxT<>  AxisAlignedBox;
 typedef Private::Matrix3T<>         Matrix3;
 typedef Private::Matrix4T<>         Matrix4;


} // namespace OgreMath

#endif
