#pragma once

//#include "Vector3.h"

class Ray
{
public:
	Vector3 origin;
	Vector3 direction;

	Ray() 
	{
		origin = { 0.0f, 0.0f, 0.0f };
		direction = { 0.0f, 0.0f, 0.0f };
	}

	// Without both parameters, origin defaults to 0,0,0
	Ray(Vector3 nDirection)
	{
		origin = { 0.0f, 0.0f, 0.0f };
		direction = direction;
	}

	Ray(Vector3 nOrigin, Vector3 nDirection)
	{
		origin = origin;
		direction = direction;
	}
};