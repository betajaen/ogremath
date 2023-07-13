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
