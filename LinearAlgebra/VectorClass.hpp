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
	T x()
	{
		return v[0];
	}
	T y()
	{
		assert(Size >= 2);
		return v[1];
	}
	T z()
	{
		assert(Size >= 3)
		return v[2];
	}
	T w()
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


vec3f Cross(vec3f& u, vec3f& v)
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

vec3d Cross(vec3d& v1, vec3d& v2)
{
	vec3d res;
	double t1, t2, t3, t4;
	t1 = u[0] - u[1];	t2 = v[1] + v[2];	t3 = u[0] * v[2];	t4 = t1*t2 - t3;
	//res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	//res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	//res[2] = v1[0] * v2[1] - v1[1] * v2[0];
	res[0] = v[1] * (t1 - u[2]) - t4;
	res[1] = u[2] * v[0] - t3;
	res[2] = t4 - u[1] * (v[0] - t2);
	return res;
}

float SignedVolume(vec3f& u, vec3f& v, vec3f& w)
{
	return u.Dot(Cross(v, w));
}


//compute barycentric coordinates u,v,w for
//P with respect to triangle(A,B,C);
void Barycentric(vec2f& A, vec2f& B, vec2f& C, vec2f& P, float& u, float& v, float& w)
{
	vec2f v0 = B - A, v1 = C - A, v2 = P - A;
	float d00 = v0.sqmodule(), d01 = v0.Dot(v1), d11 = v1.sqmodule();
	float d20 = v2.Dot(v0), d21 = v2.Dot(v1);
	float det = d00*d11 - d01*d01;
	v = (d11*d20 - d01*d21) / det;
	w = (d00*d21 - d01*d20) / det;
	u = 1.0f - v - w;
}
//compute barycentric coordinates u,v,w for all indicated
//Points with respect to triangle(A,B,C);
void Barycentric(vec2f& A, vec2f& B, vec2f& C, std::vector<vec2f>& Points, std::vector<float>& u, std::vector<float>& v, std::vector<float>& w)
{
	vec2f v0 = B - A, v1 = C - A;
	float d00 = v0.sqmodule(), d01 = v0.Dot(v1), d11 = v1.sqmodule();
	float det = d00*d11 - d01*d01;
	for (uint32_t i = 0; i < Points.size(); ++i)
	{
		vec2f v2 = Points[i] - A;
		float d20 = v2.Dot(v0), d21 = v2.Dot(v1);
		v[i] = (d11*d20 - d01*d21) / det;
		w[i] = (d00*d21 - d01*d20) / det;
		u[i] = 1.0f - v[i] - w[i];
	}
}
#endif // VECTOR_CLASS_H

