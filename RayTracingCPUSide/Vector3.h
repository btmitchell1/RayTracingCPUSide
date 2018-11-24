#pragma once
#include <math.h>
#include <cmath>

const float kEpsilon = 0.5e-6f;

class Vector3
{
public:
	float x;
	float y;
	float z;

	// Default constructor
	Vector3() {}

	Vector3(float iX, float iY, float iZ)
	{
		x = iX;
		y = iY;
		z = iZ;
	}

	void Set(float iX, float iY, float iZ)
	{
		x = iX;
		y = iY;
		z = iZ;
	}

	void Normalise();

	Vector3 Cross(const Vector3&);
	float Dot(const Vector3& v);
};

//------------OPERATOR OVERLOADS------------//
inline Vector3 operator+
(
	const Vector3& v1,
	const Vector3& v2
	)
{
	return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3 operator-
(
	const Vector3& v1,
	const Vector3& v2
	)
{
	return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector3 operator*
(
	const float  s,
	const Vector3& v
	)
{
	return Vector3(v.x*s, v.y*s, v.z*s);
}
////////////////////////////////////////////////
