#ifndef _PLANE_CLASS_H_
#define _PLANE_CLASS_H_
#include <vector>
#include <assert.h>
#include "VectorClass.hpp"

template <typename T>
class Plane
{
public:
	Plane<T>() {}
	Plane<T>(const T* v) 
	{
		a = v[0];
		b = v[1];
		c = v[2];
		d = v[3];
	}
	Plane<T>(Vector<T, 4>& v) 
	{
		a = v[0];
		b = v[1];
		c = v[2];
		d = v[3];
	}
	Plane<T>(T _a, T _b, T _c, T _d) 
	{
		a = _a;
		b = _b;
		c = _c;
		d = _d;
	}
	Plane<T>(Vector<T, 3>& v, T _d) 
	{
		a = v[0];
		b = v[1];
		c = v[2];
		d = _d;
	}
	Line3D(std::initializer_list<Vector<T, 3>> list)
	{
		assert(list.size == 4);
		a = list[0];
		b = list[1];
		c = list[2];
		d = list[3];
	}

	//loads the plane from a point on the surface and a normal vector
	static Plane ConstructFromPointNormal(Vector<T, 3>&, Vector<T, 3>&)
	{}

	//loads the plane from a point on the surface and two vectors in the plane
	static Plane ConstructFromPointVectors(Vector<T, 3>&, Vector<T, 3>&, Vector<T, 3>&)
	{}

	//loads the plane from 3 points on the surface
	static Plane ConstructFromPoints(Vector<T, 3>&, Vector<T, 3>&, Vector<T, 3>&)
	{}

	// Normalization
	Plane<T> Normalize()
	{}

	T UnsignedDistance(Vector<T, 3>&)
	{}
	T SignedDistance(Vector<T, 3>&)
	{}
	Vector<T, 3> ClosestPoint(Vector<T, 3>&)
	{}

	//determines the intersect of the line defined by the points V1 and V2 with the plane.
	//If there is no intersection, Hit will be false.
	Vector<T, 3> IntersectLine(Vector<T, 3>&, Vector<T, 3>&, bool&)
	{}

	//determines the intersect of the line defined by the points V1 and V2 with the plane.
	//Returns the point of intersection.  Origin is returned if no intersection exists.
	Vector<T, 3> IntersectLine(Vector<T, 3>&, Vector<T, 3>&)
	{}

	//Paramaterize the line with the variable t such that t = 0 is V1 and t = 1 is V2.
	//returns the t for this line that lies on this plane.
	T IntersectLineRatio(Vector<T, 3>&, Vector<T, 3>&)
	{}

	Vector<T, 3> IntersectLine(Line3D<T>&)
	{}



	//dot product of a plane and a 4D vector
	static T Dot(Plane&, Vector<T, 4>&)
	{}
	//dot product of a plane and a 3D coordinate
	static T DotCoord(Plane&, Vector<T, 3>&)
	{}
	//dot product of a plane and a 3D normal
	static T DotNormal(Plane&, Vector<T, 3>&)
	{}

	static bool PlanePlaneIntersection(Plane<T>&, Plane<T>&, Line<T>&)
	{}


	__forceinline Vector<T, 3> Plane::Normal() const
	{
		return Vector<T, 3>(a, b, c);
	}
	__forceinline Plane Flip()
	{
		Plane<T> result;
		result.a = -a;
		result.b = -b;
		result.c = -c;
		result.d = -d;
		return result;
	}

#pragma region operators
	bool operator == (const Plane<T>& p) const
	{}

	bool operator != (const Plane<T>& p) const
	{}


	friend Plane operator *(T k, const Plane<T>& p)
	{}


	Plane operator * (T) const
	{}
	Plane operator / (T) const
	{}

	Plane& operator *= (T)
	{}
	Plane& operator /= (T)
	{}

	operator const T*() const
	{}
	operator T*()
	{}
#pragma endregion
	//the (a, b, c, d) in a*x + b*y + c*z + d = 0.
	T a, b, c, d;
};
#endif // !_PLANE_CLASS_H_