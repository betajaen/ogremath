# ogremath
ogremath是三维几何计算库，ogremath中的算法提取自OGRE。
## ogremath特性
- Radian  存储弧度值。
- Degree 存储角度值。
- Vector3 三维向量，表示坐标或方向。提供了求向量长度、两坐标距离、点乘、叉乘等向量运算。
- Quaternion 四元数。提供了点乘、求逆、旋转Vector3、等四元素运算。
- Matrix3 三维矩阵，3x3。
- Matrix4 四维矩阵，4x4。
- AxisAlignedBox 对齐到坐标轴的包围盒。提供了和Plane、Sphere、Ray的碰撞检测算法。
- Plane 平面。提供了映射点、点到平面距离、点位置判断等算法。
- Sphere 球。提供了和AxisAlignedBox、Plane、Vector3的碰撞检测算法。
- Ray 射线。提供了和AxisAlignedBox、Plane、Sphere的碰撞检测算法。
## 使用方法
直接把OgreMathLib.h和OgreMathLib.cpp添加到工程中，在适当的地方包含OgreMathLib.h即可。


------------


# ogremath
ogremath is a 3D geometric computing library, and the algorithms in ogremath are extracted from OGRE.
## ogremath特性
- **Radian**  Store radian values.
- **Degree** Store the Angle value.
- **Vector3** A three-dimensional vector representing coordinates or directions. Vector operations such as length of vector, distance between two coordinates, dot product and cross product are provided.
- **Quaternion** Provide point multiplication, inversion, rotation Vector3, and other four element operations.
- **Matrix3** Three-dimensional matrix, 3x3.
- **Matrix4** Four-dimensional matrix, 4x4.
- **AxisAlignedBox**  It provides collision detection algorithms with Plane, Sphere and Ray.
- **Plane** Algorithms such as mapping point, distance from point to plane and point position judgment are provided.
- **Sphere** Provides collision detection algorithms with AxisAlignedBox, Plane, and Vector3.
- **Ray** Provides collision detection algorithms with AxisAlignedBox, Plane, Sphere.
## Usage
Add OgreMathLib.h and OgreMathLib.cpp directly to the project, including OgreMathLib.h where appropriate.
