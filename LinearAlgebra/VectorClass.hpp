#ifndef VECTOR_CLASS_H
#define VECTOR_CLASS_H
#include <math.h>
#include <assert.h>
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
#pragma endregion
	T module()
	{
		T res = 0;
		for (size_t i = 0; i < Size; ++i)
			res += v[i] * v[i];
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
	T Dot(const Vector& vec)
	{
		T res = 0;
		for (size_t i = 0; i < Size; ++i)
			res += v[i] * vec.v[i];
		return res;
	}
	
};

typedef Vector<float, 2> vec2f;
typedef Vector<float, 3> vec3f;
typedef Vector<float, 4> vec4f;

typedef Vector<double, 2> vec2d;
typedef Vector<double, 3> vec3d;
typedef Vector<double, 4> vec4d;


vec3f Cross(vec3f& v1, vec3f& v2)
{
	vec3f res;
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return res;
}

vec3d Cross(vec3d& v1, vec3d& v2)
{
	vec3d res;
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return res;
}
#endif // VECTOR_CLASS_H

