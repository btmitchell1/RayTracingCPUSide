#pragma once

//#include "Vector3.h"

class Ray
{
public:
	Vector3 mOrigin;
	Vector3 mDirection;

	Ray() {}

	// Without both parameters, origin defaults to 0,0,0
	Ray(Vector3 direction)
	{
		mOrigin = { 0.0f, 0.0f, 0.0f };
		mDirection = direction;
	}

	Ray(Vector3 origin, Vector3 direction)
	{
		mOrigin = origin;
		mDirection = direction;
	}
};