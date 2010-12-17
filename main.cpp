#include "OgreMathLib.h"

void main()
{
 OgreMathLib::Radian r;
 OgreMathLib::Matrix4 m;
 m.setTrans(OgreMathLib::Quaternion(1.0,0.0,0.0,0.0) * OgreMathLib::Vector3(1,2,3));
}
