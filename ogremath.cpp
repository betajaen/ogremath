/*
-----------------------------------------------------------------------------

OgreMath.cpp

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
-----------------------------------------------------------------------------
*/


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

  const Vector2 Vector2::ZERO( 0, 0);

  const Vector2 Vector2::UNIT_X( 1, 0);
  const Vector2 Vector2::UNIT_Y( 0, 1);
  const Vector2 Vector2::NEGATIVE_UNIT_X( -1,  0);
  const Vector2 Vector2::NEGATIVE_UNIT_Y(  0, -1);
  const Vector2 Vector2::UNIT_SCALE(1, 1);

  const Vector3 Vector3::ZERO( 0, 0, 0 );

  const Vector3 Vector3::UNIT_X( 1, 0, 0 );
  const Vector3 Vector3::UNIT_Y( 0, 1, 0 );
  const Vector3 Vector3::UNIT_Z( 0, 0, 1 );
  const Vector3 Vector3::NEGATIVE_UNIT_X( -1,  0,  0 );
  const Vector3 Vector3::NEGATIVE_UNIT_Y(  0, -1,  0 );
  const Vector3 Vector3::NEGATIVE_UNIT_Z(  0,  0, -1 );
  const Vector3 Vector3::UNIT_SCALE(1, 1, 1);

  const Vector4 Vector4::ZERO( 0, 0, 0, 0 );

 }
}