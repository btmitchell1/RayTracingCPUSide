#include "Vector3.h"

void Vector3::Normalise()
{
	float lengthSq = x * x + y * y + z * z;

	if (abs(lengthSq) < kEpsilon)
	{
		x = y = z = 0.0f;
	}
	else
	{
		float invLength = (float) 1 / sqrt(lengthSq);
		x *= invLength;
		y *= invLength;
		z *= invLength;
	}
}

Vector3 Vector3::Cross(const Vector3& v)
{
	return Vector3(y*v.z - z * v.y, z*v.x - x * v.z, x*v.y - y * v.x);
}
float Vector3::Dot(const Vector3& v)
{
	return x * v.x + y * v.y + z * v.z;
}