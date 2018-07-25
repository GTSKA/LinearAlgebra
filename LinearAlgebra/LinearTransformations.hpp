#ifndef _LINEAR_TRANSFORMATIONS_H_
#define _LINEAR_TRANSFORMATIONS_H_
#include "VectorClass.hpp"
#include "MatrixClass.hpp"
vec3f TransformPoint(Mat44f& m,const vec3f& p)
{
	vec4f p = m*vec4f(p,1.0f);
}

#endif // !_LINEAR_TRANSFORMATIONS_H_

