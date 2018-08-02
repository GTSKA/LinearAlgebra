#include <stdio.h>
#include "LinearAlgebra.hpp"
#include <iostream>
#include <vector>
int main(void)
{

	vec3f u{ 1,4,-2 };
	vec3f v{ 2,-3,-1 };
	std::cout << Cross(u, v) << std::endl;
	getchar();


}