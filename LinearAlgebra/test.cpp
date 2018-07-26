#include <stdio.h>
#include "LinearAlgebra.hpp"
#include <iostream>
#include <vector>
int main(void)
{

	Mat44f a{1,2,3,4,
			0,2,3,4,
			0,0,3,4,
			0,0,0,4};
	Mat44f b = a.Inverse();
	std::cout << b << std::endl;
	std::cout << b*a << std::endl;
	getchar();


}