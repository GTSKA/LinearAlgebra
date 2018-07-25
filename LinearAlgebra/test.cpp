#include <stdio.h>
#include "LinearAlgebra.hpp"
int main(void)
{
	Matrix<float, 4, 3> a;
	Matrix<float, 3, 4> b;
	Matrix<float, 4, 4> c;
	a.SetIdentity();
	b.SetIdentity();
	c = a*b;
	for (size_t i = 0; i < c.rowSize(); ++i)
	{
		for (size_t j = 0; j < c.colSize(); ++j)
			printf("%f\t", c[i][j]);
		printf("\n");
	}
	getchar();


}