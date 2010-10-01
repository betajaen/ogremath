/*
-----------------------------------------------------------------------------

OgreMath.h

This source file was part of OGRE
(Object-oriented Graphics Rendering Engine)
For the latest info, see http://www.ogre3d.org/

Copyright (c) 2000-2009 Torus Knot Software Ltd

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
  template<int=0> class MathT;
  template<int=0> class DegreeT;
  template<int=0> class RadianT;
  template<int=0> class Vector2T;
  template<int=0> class Vector3T;
  template<int=0> class Vector4T;
  template<int=0> class QuaternionT;
  template<int=0> class RayT;
  template<int=0> class SphereT;
  template<int=0> class PlaneT;
  template<int=0> class AxisAlignedBoxT;
  template<int=0> class Matrix3T;
  template<int=0> class Matrix4T;
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
  template<int o_O> class MathT
  {
  public:

   typedef Private::MathT<o_O>            Math;
   typedef Private::DegreeT<o_O>          Degree;
   typedef Private::RadianT<o_O>          Radian;
   typedef Private::Vector2T<o_O>         Vector2;
   typedef Private::Vector3T<o_O>         Vector3;
   typedef Private::Vector3T<o_O>         Vector4;
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

  };


 } // namespace Private


 typedef Private::MathT<> Math;

} // namespace OgreMath

#endif
