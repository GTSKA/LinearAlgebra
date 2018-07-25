#ifndef _MATRIX_CLASS_H
#define _MATRIX_CLASS_H
#include <assert.h>
#include "VectorClass.hpp"
template<typename T, size_t m, size_t n>
class Matrix
{
public:
	T** mat;
	Matrix()
	{
		mat = new T*[m];
		for (size_t i = 0; i < m; ++i)
		{
			mat[i] = new T[n];
			for (size_t j = 0; j < n; ++j)
				mat[i][j] = 0;
		}
	}
	Matrix(const Matrix<T, m, n>& mat2)
	{
		mat = new T*[m];
		for (size_t i = 0; i < m; ++i)
		{
			mat[i] = new T[n];
			for (size_t j = 0; j < n; ++j)
				mat[i][j] = mat2.mat[i][j];
		}
	}
	~Matrix()
	{
		if (mat != nullptr)
		{
			for (size_t i = 0; i < m; ++i)
				delete[] mat[i];
			delete[] mat;
		}
	}
	size_t rowSize()
	{
		return m;
	}
	size_t colSize()
	{
		return n;
	}
#pragma region operators
	T* operator [](size_t i)
	{
		return mat[i];
	}
	Matrix<T, m, n>& operator =(const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				mat[i][j] = mat2.mat[i][j];
		}
		return *this;
	}
	Matrix<T, m, n>& operator +=(const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				mat[i][j] += mat2[i][j];
		}
		return *this;
	}
	Matrix<T, m, n> operator +(const Matrix<T, m, n>& mat2)
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j] + mat2.mat[i][j];
		return res;
	}
	Matrix<T, m, n>& operator -=(const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
		{			
			for (size_t j = 0; j < n; ++j)
				mat[i][j] -= mat2[i][j];
		}
		return *this;
	}
	Matrix<T, m, n> operator -(const Matrix<T, m, n>& mat2)
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j] - mat2[i][j];
		return res;
	}
	Matrix<T, m, n>& operator *=(const T k)
	{
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				mat[i][j] *= k;
		}
		return *this;
	}
	Matrix<T, m, n> operator * (const T& k)
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j] * k;
		}
		return res;
	}
	friend Matrix<T, m, n> operator *(const T& k,const Matrix<T, m, n>& mat2)
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				res[i][j] = k*mat2[i][j];
		}
		return res;
	}
	Vector<T, m> operator * (const Vector<T, m>& v)
	{
		assert(m == n);
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += mat[i][j] * v[j];
		return res;
	}
	Vector<T, m> operator * (Vector<T, m>& v)
	{
		assert(m == n);
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += mat[i][j] * v[j];
		return res;
	}
	friend Vector<T, m> operator *(const Vector<T, m>& v, const Matrix<T, m, n>& mat2)
	{
		assert(m == n);
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += v[j] * mat[j][i];
		return res;
	}
	friend Vector<T, m> operator *(Vector<T, m>& v, Matrix<T, m, n>& mat2)
	{
		assert(m == n);
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += v[j] * mat[j][i];
		return res;
	}
	template<size_t n2>
	Matrix<T, m, n2> operator *(const Matrix<T, n, n2>& mat2)
	{
		Matrix<T, m, n2> res;
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n2; ++j)
			{
				for (size_t k = 0; k < n; ++k)
				{
					res[i][j] += mat[i][k] * mat2.mat[k][j];
				}
			}
		}
		return res;
	}
	bool operator ==(const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				if (mat[i][j] != mat[i][j])
					return false;
		return true;
	}
	bool operator !=(const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				if (mat[i][j] != mat[i][j])
					return true;
		return false;
	}

#pragma endregion
	void SetIdentity()
	{
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				if (i == j)
					mat[i][j] = (T)1;
				else
					mat[i][j] = (T)0;
	}
	Matrix<T, n, m> Transpose()
	{
		Matrix<T, n, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[j][i] = mat[i][j];
		return res;
	}
	Matrix<T, m, n> SetRotationX(float angle)
	{
		assert(n >= 3);
		assert(n == m);
		Matrix<T, m, n> res;
		res.SetIdentity();
		float cosa = cosf(angle);
		float sina = sinf(angle);
		res[1][1] = cosa;
		res[1][2] = -sina;
		res[2][1] = sina;
		res[2][2] = cosa;
		return res;
	}
	Matrix<T, m, n> SetRotationY(float angle)
	{
		assert(n >= 3);
		assert(n == m);
		Matrix<T, m, n> res;
		res.SetIdentity();
		float cosa = cosf(angle);
		float sina = sinf(angle);
		res[0][0] = cosa;
		res[0][2] = sina;
		res[2][0] = -sina;
		res[2][2] = cosa;
		return res;
	}
	Matrix<T, m, n> SetRotationZ(float angle)//can rotate 2D coordinates;
	{
		assert(n >= 2);
		assert(n == m);
		Matrix<T, n, m> res;
		res.SetIdentity();
		float cosa = cosf(angle);
		float sina = sinf(angle);
		res[0][0] = cosa;
		res[0][1] = -sina;
		res[1][0] = sina;
		res[1][1] = cosa;
		return res;
	}
	Matrix<T, m, n> SetRotationAxis(float angle, Vector<T, 3>& v)
	{
		Vector<T, 3> u = v.normalize();
		assert(u.module());
		assert(n >= 3);
		assert(n == m);
		Matrix<T, m, n> res;
		float cosa = cosf(angle);
		float sina = sinf(angle);
		res.SetIdentity();
		res[0][0] = cosa + u[0]*u[0]*(1-cosa);				res[0][1] = u[0] * u[1] * (1 - cosa) - u[2]*sina;	res[0][2] = u[0] * u[2] * (1 - cosa) + u[1] * sina;
		res[1][0] = u[1] * u[0] * (1 - cosa) + u[2] * sina; res[1][1] = cosa + u[1] * u[1] * (1 - cosa);		res[1][2] = u[1] * u[2] * (1 - cosa) - u[0] * sina;
		res[2][0] = u[2] * u[0] * (1 - cosa) - u[1] * sina;	res[2][1] = u[2] * u[1] * (1 - cosa) + u[0] * sina;	res[2][2] = cosa + u[2] * u[2] * (1 - cosa);
		return res;
	}

	Matrix<T, m, n> SetScale(Vector<T, m>& v)
	{
		assert(m == n);
		Matrix<T, m, n> res;
		res.SetIdentity();
		for (size_t i = 0; i < m; ++i)
			res[i][i] *= v[i];
		return res;
	}
	Matrix<T, m, n> SetTranslation(Vector<T, m-1> p)
	{
		assert(m == n);
		Matrix<T, m, n> res;
		res.SetIdentity();
		for (size_t i = 0; i < m - 1; ++i)
			res[i][m - 1] = p[i];
		return res;
	}

};
typedef Matrix<float, 4, 4> Mat44f;
typedef Matrix<float, 3, 3> Mat33f;
Mat44f lookAtRH(vec3f& pos, vec3f& target, vec3f& up)
{
	Mat44f lookAt;
	vec3f zaxis, xaxis, yaxis;
	zaxis = (pos - target).normalize();
	xaxis = Cross(up, zaxis).normalize();
	yaxis = Cross(zaxis, xaxis);
	lookAt.SetIdentity();
	lookAt[0][0] = xaxis[0];			lookAt[0][1] = yaxis[0];			lookAt[0][2] = zaxis[0];
	lookAt[1][0] = xaxis[1];			lookAt[1][1] = yaxis[1];			lookAt[1][2] = zaxis[1];
	lookAt[2][0] = xaxis[2];			lookAt[2][1] = yaxis[2];			lookAt[2][2] = zaxis[2];
	lookAt[3][0] = -(xaxis.Dot(pos));	lookAt[3][1] = -(yaxis.Dot(pos));	lookAt[3][2] = -(zaxis.Dot(pos));
	return lookAt;
}
Mat44f lookAtLH(vec3f& pos, vec3f& target, vec3f& up)
{
	Mat44f lookAt;
	vec3f zaxis, xaxis, yaxis;
	zaxis = (target-pos).normalize();
	xaxis = Cross(up, zaxis).normalize();
	yaxis = Cross(zaxis, xaxis);
	lookAt.SetIdentity();
	lookAt[0][0] = xaxis[0];			lookAt[0][1] = yaxis[0];			lookAt[0][2] = zaxis[0];
	lookAt[1][0] = xaxis[1];			lookAt[1][1] = yaxis[1];			lookAt[1][2] = zaxis[1];
	lookAt[2][0] = xaxis[2];			lookAt[2][1] = yaxis[2];			lookAt[2][2] = zaxis[2];
	lookAt[3][0] = -(xaxis.Dot(pos));	lookAt[3][1] = -(yaxis.Dot(pos));	lookAt[3][2] = -(zaxis.Dot(pos));
	return lookAt;
}
Mat44f perspectiveRH(float w, float h, float near, float far)
{
	assert(w && h);
	assert(far != near);
	Mat44f res;
	res[0][0] = 2 * near / w;
	res[1][1] = 2 * near / h;
	res[2][2] = far / (near - far); res[2][3] = -1;
	res[3][2] = near*far / (near - far);
	return res;
}
Mat44f perspectiveLH(float w, float h, float near, float far)
{
	assert(w && h);
	assert(far != near);
	Mat44f res;
	res[0][0] = 2 * near / w;
	res[1][1] = 2 * near / h;
	res[2][2] = far / (far - near); res[2][3] = 1;
	res[3][2] = near*far / (near - far);
	return res;
}
Mat44f perspectiveFovRH(float fov, float aspectRatio, float near, float far)
{
	assert(fov);
	assert(far != near);
	float h = 1.0f / tan(fov * 0.5f);
	float w = h * aspectRatio;
	Mat44f res;
	res[0][0] = w;
	res[1][1] = h;
	res[2][2] = far / (near-far); res[2][3] = -1;
	res[3][2] = near*far / (near - far);
	return res;
}
Mat44f perspectiveFovLH(float fov, float aspectRatio, float near, float far)
{
	assert(fov);
	assert(far != near);
	float h = 1.0f / tan(fov * 0.5f);
	float w = h * aspectRatio;
	Mat44f res;
	res[0][0] = w;
	res[1][1] = h;
	res[2][2] = far / (far-near); res[2][3] = 1;
	res[3][2] = near*far / (near - far);
	return res;
}
Mat44f perspectiveOffCenterRH(float left, float right, float bottom, float top, float near, float far)
{
	assert(right != left);
	assert(top != bottom);
	assert(far != near);
	Mat44f res;
	res[0][0] = 2 * near / (right - left);
	res[1][1] = 2 * near / (top - bottom);
	res[2][0] = (left + right) / (right - left); res[2][1] = (top + bottom) / (top - bottom); res[2][2] = far / (near - far); res[2][3] = -1;
	res[3][2] = near*far / (near - far);
	return res;
}
Mat44f perspectiveOffCenterLH(float left, float right, float bottom, float top, float near, float far)
{
	assert(right != left);
	assert(top != bottom);
	assert(far != near);
	Mat44f res;
	res[0][0] = 2 * near / (right - left);
	res[1][1] = 2 * near / (top - bottom);
	res[2][0] = (left + right) / (right - left); res[2][1] = (top + bottom) / (top - bottom); res[2][2] = far / (far - near); res[2][3] = 1;
	res[3][2] = near*far / (near - far);
	return res;
}
typedef Matrix<double, 4, 4> Mat44d;
typedef Matrix<double, 3, 3> Mat33d;
#endif // !_MATRIX_CLASS_H

