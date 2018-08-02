#ifndef _LINE_CLASS_H_
#define _LINE_CLASS_H_
#include <vector>
#include <assert.h>
#include "VectorClass.hpp"
#include "PlaneClass.hpp"
template<typename T>
Line3D
{
public:
	Vector<T,3> P0;
	Vector<T, 3> Dir;
	Line3D() {} |
	Line3D(std::initializer_list<Vector<T, 3>> list)
	{
		assert(list.size == 2);
		P0 = list[0];
		Dir = list[1];
	}
	Line3D(Line3D& line)
	{
		P0 = line.P0;
		Dir = line.Dir;
	}
#pragma region operators
	Line3D& operator = (Line3D& line)
	{
		P0 = line.P0;
		Dir = line.Dir;
	}
#pragma endregion
	T DistToPoint(const Vector<T, 3>& P) 
	{
		T t0 = Dir.Dot(P - P0) / Dir.sqmodule;
		return P.Dist(P0+t0*Dir);
	}
	bool ClosestPoint(const Line3D<T>& L2, Vector<T, 3>& res1, Vector<T, 3> &res2)
	{
		Vector<T, 3> NewDir = Cross(Dir, L2.Dir);
		if (NewDir.module() == 0)
			return false;
		Plane P1 = Plane<T>::ConstructFromPointVectors(P0, NewDir, Dir);
		Plane P2 = Plane<T>::ConstructFromPointVectors(PL2.0, NewDir, L2.Dir);
		res1 = P2.IntersectLine(L1);
		res2 = P2.IntersectLine(L2);
		return true;
	}
	T Dist(const Line3D& L2)
	{
		Vector<T, 3> cross = Cross(Dir, L2.Dir);
		return abs(cross.Dot(L2.P0 - P0) / cross.module());
	}
	T DistSq(const Line3D& L2)
	{
		Vector<T, 3> cross = Cross(Dir, L2.Dir);
		float dot = cross.Dot(L2.P0 - P0);
		return dot*dot / cross.sqmodule();
	}
	Vector<T, 3> linePoint(T t0)
	{
		Vector<T, 3> res;
		res = P0 + dir*t0;
		return res;
	}
};
#endif //_LINE_CLASS_H_
