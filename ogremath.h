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

  typedef vector< PlaneT<void> >::type PlaneList;
  /** @} */
  /** @} */


 } // namespace Private


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
