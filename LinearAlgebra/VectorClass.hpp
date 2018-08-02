#ifndef VECTOR_CLASS_H
#define VECTOR_CLASS_H
#include <math.h>
#include <assert.h>
#include <iostream>
template<typename T, size_t Size>
class Vector
{
public:
	T* v;

	Vector()
	{
		v = new T[Size];
		for (size_t i = 0; i < Size; ++i)
			v[i] = (T)0;
	}

	Vector(const Vector& vec)
	{
		v = new T[Size];
		for (size_t i = 0; i < Size; ++i)
			v[i] = vec.v[i];
	}
	Vector(std::initializer_list<T> list)
	{
		assert(list.size() <= Size);
		v = new T[Size];
		for(size_t i = 0;i<Size;++i)
			v[i] = (T)0;
		size_t i = 0;
		for (const auto& l : list)
		{
			v[i] = l;
			i++;
		}
	}
	Vector(const Vector<T, Size - 1>& vec, T val)
	{
		assert(Size - 1 > 0);
		v = new T[Size];
		for (size_t i = 0; i < Size-1; ++i)
			v[i] = vec.v[i];
		v[Size - 1] = val;
	}
	Vector(const Vector<T, Size - 2>& vec, T val1, T val2)
	{
		assert(Size - 2 > 0);
		v = new T[Size];
		for (size_t i = 0; i < Size - 2; ++i)
			v[i] = vec.v[i];
		v[Size - 2] = val1;
		v[Size - 1] = val2;
	}
	Vector(const T* values)
	{
		v = new T[Size];
		for (size_t i = 0; i < Size; ++i)
			v[i] = values[i];
	}
	~Vector()
	{
		if (v != nullptr)
			delete[] v;
	}
#pragma region operators
	T& operator[] (size_t i)
	{
		return v[i];
	}

	Vector<T, Size>& operator =(const Vector<T, Size>& vec)
	{
		for (size_t i = 0; i < Size; ++i)
			v[i] = vec.v[i];
		return *this;
	}

	Vector<T, Size>& operator =(const T* vec)
	{
		for (size_t i = 0; i < Size; ++i)
			v[i] = vec[i];
		return *this;
	}

	Vector<T, Size> operator +(const Vector<T, Size>& v2)
	{
		Vector<T, Size> res;
		for (size_t i = 0; i < Size; ++i)
			res[i] = v[i] + v2.v[i];
		return res;
	}
	Vector<T, Size> operator -(const Vector<T, Size>& v2)
	{
		Vector<T, Size> res;
		for (size_t i = 0; i < Size; ++i)
			res[i] = v[i] - v2.v[i];
		return res;
	}
	Vector<T, Size>& operator +=(const Vector<T, Size>& v2)
	{
		for (size_t i = 0; i < Size; ++i)
			v[i] += v2.v[i];
		return *this;
	}
	Vector<T, Size>& operator -=(const Vector<T, Size>& v2)
	{
		for (size_t i = 0; i < Size; ++i)
			v[i] -= v2.v[i];
		return *this;
	}
	bool operator ==(const Vector<T, Size>& vec)
	{
		for (size_t i = 0; i < Size; ++i)
			if (v[i] != vec.v[i])
				return false;
		return true;
	}
	bool operator !=(const Vector<T, Size>& vec)
	{
		for (size_t i = 0; i < Size; ++i)
			if (v[i] != vec.v[i])
				return true;
		return false;
	}
	Vector<T, Size> operator-()
	{
		Vector<T, Size> res;
		for (size_t i = 0; i < Size; ++i)
			res[i] = -v[i];
		return res;
	}
	Vector<T, Size> operator *(T k)
	{
		Vector<T, Size> res;
		for (size_t i = 0; i < Size; ++i)
			res[i] = v[i] * k;
		return res;
	}
	Vector<T, Size>& operator *=(T k)
	{
		for (size_t i = 0; i < Size; ++i)
			v[i] *= k;
		return *this;
	}
	friend Vector<T, Size> operator*(T k, const Vector<T, Size>& vec)
	{
		Vector<T, Size> res;
		for (size_t i = 0; i < Size; ++i)
			res[i] = k * vec.v[i];
		return res;
	}
	friend std::ostream& operator<<(std::ostream& os, Vector<T, Size> vec)
	{
		os << "(";
		for (size_t i = 0; i < Size; ++i)
		{
			os << vec.v[i];
			if (i < Size - 1)
				os << ", ";
		}
		os << ")" << std::endl;
		return os;
	}
#pragma endregion
	//sqrt(x1*x1+x2*x2+...+xn*xn)
	T module()
	{
		T res = 0;
		for (size_t i = 0; i < Size; ++i)
			res += v[i] * v[i];
		res = sqrt(res);
		return res;
	}
	// x1*x1+x2*x2+...+xn*xn
	T sqmodule()
	{
		T res = 0;
		for (size_t i = 0; i < Size; ++i)
			res += v[i] * v[i];
		return res;
	}
	T Dist(const Vector<T, Size>& vec)
	{
		T res;
		for (size_t i; i < Size; ++i)
		{
			float aux = vec.v[i] - v[i];
			res += aux*aux;
		}
		res = sqrt(res);
		return res;
	}
	Vector<T, Size> normalize()
	{
		Vector<T, Size> res;
		T fmod = module();
		if (fmod == (T)0)
			return res;
		for (size_t i = 0; i < Size; ++i)
			res[i] /= fmod;
		return res;
	}
#pragma region DirectionCosines
	T cosa()
	{
		return v[0] / module();
	}
	T cosb()
	{
		if (Size > 1)
			return v[1] / module();
		return 0;
	}
	T cosg()
	{
		if (Size > 2)
			return v[2] / module();
		return 0;
	}
#pragma endregion
	T Dot(const Vector<T,Size>& vec)
	{
		T res = 0;
		for (size_t i = 0; i < Size; ++i)
			res += v[i] * vec.v[i];
		return res;
	}
	//decompose vector v in vectors p & q, which p is parallel to the object vector and q perpendicular
	void Decompose(Vector<T, Size>& V, Vector<T, Size>& p, Vector<T, Size>& q)
	{
		T uv = Dot(v);
		T uu = sqmodule();
		p = (uv / uu)*(*this);
		q = v - p;
	}
#pragma region Subvectors
	inline T x()
	{
		return v[0];
	}
	inline T y()
	{
		assert(Size >= 2);
		return v[1];
	}
	inline T z()
	{
		assert(Size >= 3)
		return v[2];
	}
	inline T w()
	{
		assert(Size >= 4);
		return v[3];
	}
	Vector<T, 2> xy()
	{
		assert(Size >= 2);
		Vector<T, 2> res;
		res[0] = v[0];
		res[1] = v[1];
		return res;
	}
	Vector<T, 2> xz()
	{
		assert(Size >= 3);
		Vector<T, 2> res;
		res[0] = v[0];
		res[1] = v[2];
		return res;
	}
	Vector<T, 2> yz()
	{
		assert(Size >= 3);
		Vector<T, 2> res;
		res[0] = v[1];
		res[1] = v[2];
		return res;
	}
	Vector<T, 2> xw()
	{
		assert(Size >= 4);
		Vector<T, 2> res;
		res[0] = v[2];
		res[1] = v[3];
		return res;
	}
	Vector<T, 3> xyz()
	{
		assert(Size >= 3);
		Vector<T, 3> res;
		res[0] = v[0];
		res[1] = v[1];
		res[2] = v[2];
		return res;
	}
#pragma endregion

};

typedef Vector<float, 2> vec2f;
typedef Vector<float, 3> vec3f;
typedef Vector<float, 4> vec4f;

typedef Vector<double, 2> vec2d;
typedef Vector<double, 3> vec3d;
typedef Vector<double, 4> vec4d;

template<typename T>
Vector<T,3> Cross(Vector<T, 3>& u, Vector<T, 3>& v)
{
	vec3f res;
	float t1, t2, t3, t4;
	t1 = u[0] - u[1];	t2 = v[1] + v[2];	t3 = u[0] * v[2];	t4 = t1*t2 - t3;
	//res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	//res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	//res[2] = v1[0] * v2[1] - v1[1] * v2[0];
	res[0] = v[1] * (t1 - u[2]) - t4;
	res[1] = u[2] * v[0] - t3;
	res[2] = t4 - u[1] * (v[0] - t2);
	return res;
}


template<typename T>
T SignedVolume(Vector<T, 3>& u, Vector<T, 3>& v, Vector<T, 3>& w)
{
	return u.Dot(Cross(v, w));
}


//compute barycentric coordinates u,v,w for
//P with respect to triangle(A,B,C);
template<typename T>
void Barycentric(Vector<T, 2>& A, Vector<T, 2>& B, Vector<T, 2>& C, Vector<T, 2>& P, T& u, T& v, T& w)
{
	vec2f v0 = B - A, v1 = C - A, v2 = P - A;
	T d00 = v0.sqmodule(), d01 = v0.Dot(v1), d11 = v1.sqmodule();
	T d20 = v2.Dot(v0), d21 = v2.Dot(v1);
	T det = d00*d11 - d01*d01;
	v = (d11*d20 - d01*d21) / det;
	w = (d00*d21 - d01*d20) / det;
	u = T(1.0) - v - w;
}

//compute barycentric coordinates u,v,w for all indicated
//Points with respect to triangle(A,B,C);
template<typename T>
void Barycentric(Vector<T, 2>& A, Vector<T, 2>& B, Vector<T, 2>& C, std::vector<Vector<T, 2>>& Points,
	std::vector<T>& u, std::vector<T>& v, std::vector<T>& w)
{
	vec2f v0 = B - A, v1 = C - A;
	T d00 = v0.sqmodule(), d01 = v0.Dot(v1), d11 = v1.sqmodule();
	T det = d00*d11 - d01*d01;
	for (uint32_t i = 0; i < Points.size(); ++i)
	{
		vec2f v2 = Points[i] - A;
		T d20 = v2.Dot(v0), d21 = v2.Dot(v1);
		v[i] = (d11*d20 - d01*d21) / det;
		w[i] = (d00*d21 - d01*d20) / det;
		u[i] = 1.0f - v[i] - w[i];
	}
}

template<typename T>
inline T TriArea2D(Vector<T, 2> v1, Vector<T, 2> v2, Vector<T, 2> v3)
{
	return (v1.x() - v2.x())*(v2.y() - v3.y()) - (v2.x() - v3.x())*(v1.y()*v2.y());
}
//compute barycentric coordinates (u,v,w) for point p
//with respect to triangle (a,b,c)
template<typename T>
void Barycentric(Vector<T, 3>& a, Vector<T, 3>& b, Vector<T, 3>& c, Vector<T, 3>& p,
	T& u, T& v, T& w)
{
	//Unnormalized triangle normal
	vec3f m = Cross(b - a, c - a);
	// Nominators and one-over-denominator for u and v ratios
	T nu, nv, ood;
	//absolute component for determining projection plane
	T x = abs(m.x()), y = abs(m.y()), z = abs(m.z());
	
	//compute areas in plane of largest projection
	if (x >= y && x >= z)
	{
		//x is the largest, project to the yz plane
		nu = TriArea2D(p.yz(), b.yz(), c.yz());//Area of PBC in yz plane
		nv = TriArea2D(p.yz(), c.yz(), a.yz());//Area of PCA in yz plane
		ood = 1.0f / m.x();						//1/(2*area of ABC in yz plane)
	}
	else if(y >=x&&y>=z)
	{
		//y is the largest, project to the xz plane
		nu = TriArea2D(p.xz(), b.xz(), c.xz());
		nv = TriArea2D(p.xz(), c.xz(), a.xz());
		ood = 1.0f / -m.y;
	}
	else 
	{
		//z is the largest, project to xy plane
		nu = TriArea2D(p.xy(), b.xy(), c.xy());
		nv = TriArea2D(p.xy(), c.xy(), a.xy());
		ood = 1.0f / m.z;
	}
	u = nu*ood;
	v = nv*ood;
	w = T(1.0) - u - v;
}
//Test if point p is contained in triangle (a,b,c)
template<typename T>
bool TestPointTriangle(Vector<T, 3>& p, Vector<T, 3>& a, Vector<T, 3>& b, Vector<T, 3>& c)
{
	float u, v, w;
	Barycentric(a, b, c, p, u, v, w);
	return v >= 0.0f && w >= 0.0f && (v + w) <= 1.0f;
}
#endif // VECTOR_CLASS_H

