#ifndef _MATRIX_CLASS_H
#define _MATRIX_CLASS_H
#include <assert.h>
#include "VectorClass.hpp"
#include <iostream>
#include <vector>
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
	Matrix(std::initializer_list<T> list)
	{
		assert(list.size() == m*n);
		mat = new T*[m];
		for (size_t i = 0; i < m; ++i)
		{
			mat[i] = new T[n];
		}
		size_t i = 0, j = 0;
		for (const auto& l : list)
		{
			mat[i][j] = l;
			j++;
			if (j == n)
			{
				j = 0;
				i++;
			}
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
				mat[i][j] -= mat2.mat[i][j];
		}
		return *this;
	}
	Matrix<T, m, n> operator -(const Matrix<T, m, n>& mat2)
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j] - mat2.mat[i][j];
		return res;
	}
	Matrix<T, m, n> operator-()
	{
		Matrix<T, m, n> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i][j] = -mat[i][j];
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
	Vector<T, m> operator * (Vector<T, n>& v)
	{
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i] += mat[i][j] * v[j];
		return res;
	}
	friend Vector<T, n> operator *(const Vector<T, m>& v, const Matrix<T, m, n>& mat2)
	{
		Vector<T, n> res;
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += v[j] * mat2[j][i];
		return res;
	}
	friend Vector<T, n> operator *(Vector<T, m>& v, Matrix<T, m, n>& mat2)
	{
		Vector<T, n> res;
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i] += v[j] * mat2[j][i];
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

	friend std::ostream& operator<<(std::ostream& os, const Matrix<T, m, n>& mat2)
	{
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				os << mat2.mat[i][j] << "\t";
			os << std::endl;
		}
		return os;
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
	void Transpose()
	{
		Matrix<T, n, m> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[j][i] = mat[i][j];
		return res;
	}
	void SetRotationX(float angle)
	{
		assert(n >= 3);
		assert(n == m);
		float cosa = cosf(angle);
		float sina = sinf(angle);
		SetIdentity();
		mat[1][1] = cosa;
		mat[1][2] = -sina;
		mat[2][1] = sina;
		mat[2][2] = cosa;
	}
	void SetRotationY(float angle)
	{
		assert(n >= 3);
		assert(n == m);
		SetIdentity();
		float cosa = cosf(angle);
		float sina = sinf(angle);
		mat[0][0] = cosa;
		mat[0][2] = sina;
		mat[2][0] = -sina;
		mat[2][2] = cosa;
	}
	void SetRotationZ(float angle)//can rotate 2D coordinates;
	{
		assert(n >= 2);
		assert(n == m);
		SetIdentity();
		float cosa = cosf(angle);
		float sina = sinf(angle);
		mat[0][0] = cosa;
		mat[0][1] = -sina;
		mat[1][0] = sina;
		mat[1][1] = cosa;
	}
	void SetRotationAxis(float angle, const Vector<T, 3>& v)
	{
		Vector<T, 3> u = v.normalize();
		assert(u.module());
		assert(n >= 3);
		assert(n == m);
		float cosa = cosf(angle);
		float sina = sinf(angle);
		SetIdentity();
		mat[0][0] = cosa + u[0]*u[0]*(1-cosa);				mat[0][1] = u[0] * u[1] * (1 - cosa) - u[2]*sina;	mat[0][2] = u[0] * u[2] * (1 - cosa) + u[1] * sina;
		mat[1][0] = u[1] * u[0] * (1 - cosa) + u[2] * sina; mat[1][1] = cosa + u[1] * u[1] * (1 - cosa);		mat[1][2] = u[1] * u[2] * (1 - cosa) - u[0] * sina;
		mat[2][0] = u[2] * u[0] * (1 - cosa) - u[1] * sina;	mat[2][1] = u[2] * u[1] * (1 - cosa) + u[0] * sina;	mat[2][2] = cosa + u[2] * u[2] * (1 - cosa);
	}

	void SetScale(const Vector<T, m>& v)
	{
		assert(m == n);
		SetIdentity();
		for (size_t i = 0; i < m; ++i)
			mat[i][i] *= v[i];
	}
	void SetTranslation(const Vector<T, m-1>& p)
	{
		assert(m == n);
		res.SetIdentity();
		for (size_t i = 0; i < m - 1; ++i)
			mat[i][m - 1] = p[i];
	}
	Matrix<T, m, n + 1> AdjuntCol(Vector<T, m>& v)
	{
		Matrix<T, m, n + 1> res;
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j];
			res[i][n] = v[i];
		}
		return res;
	}
	Matrix<T, m + 1, n> AdjuntRow(Vector<T, n>& v)
	{
		Matrix<T, m + 1, n> res;
		for (size_t i = 0; i < m; ++i)
			for (size_t j = 0; j < n; ++j)
				res[i][j] = mat[i][j];
		for (size_t j = 0; j < n; ++j)
			res[m][j] = v[j];
		return res;
	}
	Vector<T, n> getRow(size_t i)
	{
		assert(i >= 0 && i < m);
		Vector<T, n> res;
		for (size_t j = 0; j < n; ++j)
			res[j] = mat[i][j];
		return res;
	}
	Vector<T, m> getCol(size_t j)
	{
		assert(i >= 0 && i < n);
		Vector<T, m> res;
		for (size_t i = 0; i < m; ++i)
			res[i] = mat[i][j];
		return res;
	}
	T determinant()
	{
		assert(m == n);
		if (m == 1)
			return mat[0][0];
		T det;
		std::vector<bool> col;
		col.resize(m);
		for (size_t i = 0; i < m; ++i)
			col[i] = true;
		det = detAux(col);
		return det;
	}
	private:
	T detAux(std::vector<bool>& col)
	{
		size_t cont, cont2;
		T aux, det = T(0);
		cont = 0;
		for (size_t i = 0; i < col.size(); ++i)
			if (col[i])
				++cont;
		if (cont > 2)
		{
			cont2 = 0;
			for (size_t i = 0; i < col.size(); ++i)
			{
				if (col[i])
				{
					if (mat[col.size() - cont][i] != 0)
					{
						std::vector<bool> cols;
						cols.resize(col.size());
						for (size_t j = 0; j < col.size(); ++j)
							cols[j] = col[j];
						cols[i] = false;
						T powMinusOne = (cont2 % 2) ? (T)(-1.0) : (T)1.0;
						aux = mat[col.size() - cont][i] * detAux(cols) * powMinusOne;
						det += aux;
					}
					++cont2;
				}
			}
			return det;
		}
		else
		{
			T mataux[2][2];
			cont2 = 0;
			for(size_t i =0;i<col.size();++i)
				if (col[i])
				{
					mataux[0][cont2] = mat[col.size() - 2][i];
					mataux[1][cont2] = mat[col.size() - 1][i];
					++cont2;
				}
			return mataux[0][0] * mataux[1][1] - mataux[1][0] * mataux[0][1];
		}
	}
	public:
	Matrix<T, m, n> Inverse()
	{
		assert(m == n);
		T det = determinant();
		if (det == 0)
			return *this;
		if (m == 1)
		{
			Matrix<T, m, n> res;
			res[0][0] = T(1) / mat[0][0];
			return res;
		}
		if (m == 2)
		{
			Matrix<T, m, n> res;
			
			res[0][0] = mat[1][1]/det;	res[0][1] = -mat[0][1]/det;
			res[1][0] = -mat[1][0]/det;	res[1][1] = mat[0][0]/det;
			return res;
		}
		if (m == 3)
		{
			T a, b, c, d, e, f, g, h, i;
			a = mat[0][0]; b = mat[0][1]; c = mat[0][2]; 
			d = mat[1][0]; e = mat[1][1]; f = mat[1][2];
			g = mat[2][0]; h = mat[2][1]; i = mat[2][2];
			Matrix<T, m, n> res;
			res[0][0] = e*i - f*h;		res[0][1] = -(b*i - c*h);	res[0][2] = b*f - c*e;
			res[1][0] = -(d*i - f*g);	res[1][1] = a*i - c*g;		res[1][2] = -(a*f - c*d);
			res[2][0] = d*h - e*g;		res[2][1] = -(a*h - b*g);	res[2][2] = a*c - b*d;
			res *= T(1) / det;
			return res;
		}
		if (m == 4)
		{
			Matrix<T, 2, 2> A, B, C, D, A2,B2,C2,D2,DI,A_BDIC_1;//(A-BD.Inverse()*C).Inverse()
			
			for (size_t i = 0; i < 2; ++i)
			{
				for (size_t j = 0; j < 2; ++j)
				{
					A[i][j] = mat[i][j];
					B[i][j] = mat[i][j + 2];
					C[i][j] = mat[i + 2][j];
					D[i][j] = mat[i + 2][j + 2];
				}
			}
			A_BDIC_1 = (A - B*(D.Inverse())*C).Inverse();
			DI = D.Inverse();
			A2 = A_BDIC_1;
			B2 = -A_BDIC_1*B*DI;
			C2 = -DI*C*A_BDIC_1;
			D2 = DI + DI*C*A_BDIC_1*B*DI;
			Matrix<T, m, n> res;
			for (size_t i = 0; i < 2; ++i)
			{
				for (size_t j = 0; j < 2; ++j)
				{
					res[i][j] = A2[i][j];
					res[i][j + 2] = B2[i][j];
					res[i + 2][j] = C2[i][j];
					res[i + 2][j + 2] = D2[i][j];
				}
			}
			return res;
		}
		if (m > 4)
		{
			//TO DO, still this class is recomendated for low size matrix only
			return *this;
		}
	}
};
typedef Matrix<float, 4, 4> Mat44f;
typedef Matrix<float, 3, 3> Mat33f;
typedef Matrix<float, 2, 2> Mat22f;
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
//true:oriented counterclockwise(Right hand rule) 
//false:oriented clockwise(Left handed)
enum Orientation
{
	CW,//clockwise
	CCW,//counterclockwise
	COLLINEAR,//Colineal
	COPLANAR,//coplanar
	ABOVE,//above the plane
	BELOW, //below the plane
	INSIDE,
	OUTSIDE,
	COCIRCULAR,
	COSPHERICAL,
	ERROR
};
Orientation Orient2D(vec2f& a, vec2f& b, vec2f& c)
{
	Mat33f m{ a[0], a[1], 1.0f,
			  b[0], b[1], 1.0f,
			  c[0], c[1], 1.0f };
	float det = m.determinant();
	if (det > 0) return Orientation::CCW;
	if (det < 0) return Orientation::CW;
	return Orientation::COLLINEAR;
}
float SignedArea(vec2f& a, vec2f& b, vec2f& c)
{
	Mat33f m{ a[0], a[1], 1.0f,
		b[0], b[1], 1.0f,
		c[0], c[1], 1.0f };
	return m.determinant();
}
Orientation Orient3D(vec3f& a, vec3f& b, vec3f& c, vec3f& d)
{
	Mat44f m{ a[0],a[1],a[2],1.0f,
			  b[0],b[1],b[2],1.0f,
			  c[0],c[1],c[2],1.0f,
			  d[0],d[1],d[2],1.0f };
	float det = m.determinant();
	if (det < 0)return Orientation::ABOVE;
	if (det > 0)return Orientation::BELOW;
	return Orientation::COPLANAR;
}
Orientation incircle2D(vec2f&a, vec2f& b, vec2f& c, vec2f& d)
{
	Orientation p3o;//first three point orientation
	p3o = Orient2D(a, b, c);
	if (p3o == Orientation::COLLINEAR)
		return Orientation::ERROR;
	Mat44f m{ a[0],a[1],a.sqmodule(),1.0f,
				b[0],b[1],b.sqmodule(),1.0f,
				c[0],c[1],c.sqmodule(),1.0f,
				d[0],d[1],d.sqmodule(),1.0f };
	float det = m.determinant();
	if(p3o>0)
	{
		if (det < 0) return Orientation::OUTSIDE;
		if (det > 0) return Orientation::INSIDE;	
	}
	else
	{
		if (det < 0) return Orientation::INSIDE;
		if (det > 0) return Orientation::OUTSIDE;
	}
	return Orientation::COCIRCULAR;
}
Orientation insphere(vec3f& a, vec3f& b, vec3f& c, vec3f& d, vec3f& e)
{
	Orientation p4o;//first four points orientation
	p4o = Orient3D(a, b, c, d);
	if (p4o == Orientation::COPLANAR)
		return Orientation::ERROR;
	Matrix<float, 5, 5> m{ a[0],a[1],a[2],a.sqmodule(),1.0f,
							b[0],b[1],b[2],b.sqmodule(),1.0f,
							c[0],c[1],c[2],c.sqmodule(),1.0f,
							d[0],d[1],d[2],d.sqmodule(),1.0f,
							e[0],e[1],e[2],e.sqmodule(),1.0f, };
	float det = m.determinant();
	if (p4o > 0)
	{
		if (det > 0) return Orientation::INSIDE;
		if (det < 0) return Orientation::OUTSIDE;
	}
	else 
	{
		if (det > 0) return Orientation::OUTSIDE;
		if (det < 0) return Orientation::INSIDE;
	}
	return Orientation::COSPHERICAL;
}
typedef Matrix<double, 4, 4> Mat44d;
typedef Matrix<double, 3, 3> Mat33d;
typedef Matrix<double, 2, 2> Mat22d;

#endif // !_MATRIX_CLASS_H

