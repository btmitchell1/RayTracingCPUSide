#include "Vector3.h"
#include "Ray.h"
#include <fstream>
#include <algorithm>

using namespace std;

float clamp(const float &lo, const float &hi, const float &v)
{
	return max(lo, std::min(hi, v));
}

bool RayTriangleIntersection
(
	const Vector3 &origin, const Vector3 &direction,
	const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
	float &t, float &u, float &v
)
{
	/*
		To start, we must first find P

		P = O + tR
		where
		P = the pont where the ray intersects the plane of the triangle
		O = the origin of the ray (Vecotr3 origin)
		t = the distance from the ray's origin to the point of intersection
		R = direction of the ray (Vector3 direction)
	*/
	// Calculate plane normal
	Vector3 v0v1 = v1 - v0;
	Vector3 v0v2 = v2 - v0;
	Vector3 N = v0v1.Cross(v0v2);
	//N.Normalise();
	float denom = N.Dot(N);

	// If N dot R = 0 then ray is parallel to the plane and don't intersect
	// float therefore can't equal 0, use upper bound of floating point
	float NDotR = N.Dot(direction); // This is used later on if ray and plane are parallel
	if (fabs(NDotR) < kEpsilon)
	{
		return false;
	}

	/*
		Ax + By + Cz + D = 0
		where
		A, B, C = coordinates of normal to the plane
		x, y, z = coordinates of any point on the plane (v0 is a vertex of the triangle, therefore is on the plane)
		D       = distance from the origin to the plane

		The equation can be rearranged to
		D = N(A, B, C) dot v0
	*/
	float D = N.Dot(v0);

	/*
		Now that we have D, we can now calculate t using the equation by substituting
		P = O + tR
		   into
		Ax + By + Cz + D = 0

		This rearranges to
		t = - (N(A, B, C) dot O) + D  /  N(A, B, C) dot R
	*/
	t = (N.Dot(origin)) + D / NDotR;

	// Check if the triangle is behind the ray
	// If so, it shouldn't be visible
	if (t < 0.0f)
	{
		return false;
	}

	// We can now complete the original equation
	Vector3 P = origin + (t * direction);

	/*
		We must now test whether or not P is outside of the triangle. If it is not, then the ray
		does not intersect the triangle.

		Calculate one edge of the triangle using 2 vertices
		Calculate the vector between one vertex and the point P
		Calculate the cross product between both of these, the dot product that with the normal
		If the point P is inside the triangle, the dot product will be positive, (to the left of the edge)
		If the point P is outside the triangle, it will be negative, and to the right of the edge

		Do this for each edge on the triangle
	*/
	// Edge 0
	Vector3 edge0 = v1 - v0;
	Vector3 vp0 = P - v0;
	Vector3 C = edge0.Cross(vp0);
	if (N.Dot(C) < 0)
	{
		return false;
	}

	// Edge 1
	Vector3 edge1 = v2 - v1;
	Vector3 vp1 = P - v1;
	C = edge1.Cross(vp1);
	if (N.Dot(C) < 0)
	{
		return false;
	}
	u = N.Dot(C);

	// Edge 2
	Vector3 edge2 = v0 - v2;
	Vector3 vp2 = P - v2;
	C = edge2.Cross(vp2);
	if (N.Dot(C) < 0)
	{
		return false;
	}
	v = N.Dot(C);

	u /= denom;
	v /= denom;

	return true;
}

int main()
{
	// Ray origin starts at the camera origin
	Vector3 origin = { 0.0f, 0.0f, 0.0f };

	Vector3 v0(-1, -1, -3); // Bottom left
	Vector3 v1(1, -1, -3);  // Bottom right
	Vector3 v2(0, 1, -3);   // Middle

	// Determines image size
	const int screenWidth = 640;
	const int screenHeight = 480;

	float scale = 1.0f;
	float aspectRatio = (float)screenWidth / (float)screenHeight;

	// Colours (R, G, B)
	Vector3 triColours[3] = { {1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f}, {0.1f, 0.3f, 0.7f} };
	Vector3 backgroundColour = { 0.0f, 0.0f, 0.0f };

	Vector3 *framebuffer = new Vector3[screenWidth * screenHeight];
	Vector3 *pixels = framebuffer;

	// For each pixel, decide what colour it should be
	for (int i = 0; i < screenHeight; ++i)
	{
		for (int j = 0; j < screenWidth; ++j)
		{
			float x = (2 * (j + 0.5f) / (float)screenWidth - 1) * aspectRatio * scale;
			float y = (1 - 2 * (i + 0.5f) / (float)screenHeight) * scale;

			Vector3 direction(x, y, -1);
			direction.Normalise();

			float t;
			float u, v;
			if (RayTriangleIntersection(origin, direction, v0, v1, v2, t, u, v))
			{
				/*
					Using Barycentric coordinates, we can calculate a blend of the 3 vertices' colours depending
					on the position of the point

					P = uA + vB +wC
					where
					P       = pixel colour
					A, B, C = vertices colours

					u + v + w = 1

					w can be easily calculated so isn't stored as a variable
				*/
				*pixels = u * triColours[0] + v * triColours[1] + (1 - u - v) * triColours[2];
			}
			else
			{
				*pixels = backgroundColour;
			}

			++pixels;
		}
	}

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./output.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << screenWidth << " " << screenHeight << "\n255\n";
	for (uint32_t i = 0; i < screenHeight * screenWidth; ++i) {
		char r = (char)(255 * (0, 1, framebuffer[i].x));
		char g = (char)(255 * (0, 1, framebuffer[i].y));
		char b = (char)(255 * (0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}

	ofs.close();

	delete[] framebuffer;

	return 0;
}