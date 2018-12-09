#include "Vector3.h"
#include "Vector2.h"
#include "Ray.h"
#include <fstream>
#include <algorithm>
#include <time.h>

using namespace std;

// Determines image size
const int kScreenWidth = 1920;
const int kScreenHeight = 1080;

const float aspectRatio = (float)kScreenWidth / (float)kScreenHeight;

const Vector3 backgroundColour(0.0f, 0.0f, 0.0f);

// Forward Declaration
class Triangle;
bool RayTriangleIntersectionMT
(
	Triangle* triangle, const Ray* ray,
	float &t, float &u, float &v,
	bool culling
);

class Triangle
{
public:
	Vector3 v0;
	Vector3 v1;
	Vector3 v2;

	Triangle(const Vector3 iV0, const Vector3 iV1, const Vector3 iV2)
	{
		v0 = iV0;
		v1 = iV1;
		v2 = iV2;
	}
};

class Mesh
{
private:
	int mNumTriangles;
	Vector3* mVertices;
	int* mIndexBuffer;
	Vector3* mVertexNormals;

public:
	Mesh(const int numFaces, const int* faceIndex, const int* verticesIndex,
		const Vector3* vertices, Vector2* st)
	{
		int k = 0;
		int maxVertexIndex = 0;

		// find out how many triangles we need to create for this mesh
		for (int i = 0; i < numFaces; ++i) {
			mNumTriangles += faceIndex[i] - 2;
			for (int j = 0; j < faceIndex[i]; ++j)
				if (verticesIndex[k + j] > maxVertexIndex)
					maxVertexIndex = verticesIndex[k + j];
			k += faceIndex[i];
		}
		maxVertexIndex += 1;

		// allocate memory to store the position of the mesh vertices
		mVertices = new Vector3[maxVertexIndex];
		for (int i = 0; i < maxVertexIndex; ++i) {
			mVertices[i] = vertices[i];
		}

		// allocate memory to store triangle indices
		mIndexBuffer = new int[mNumTriangles * 3];
		int l = 0;
		// generate the triangle index array
		for (int i = 0, k = 0; i < numFaces; ++i) { // for each  face 
			for (int j = 0; j < faceIndex[i] - 2; ++j) { // for each triangle in the face 
				mIndexBuffer[l] = verticesIndex[k];
				mIndexBuffer[l + 1] = verticesIndex[k + j + 1];
				mIndexBuffer[l + 2] = verticesIndex[k + j + 2];
				l += 3;
			}
			k += faceIndex[i];
		}
	}

	Vector3 GetNormal(int index)
	{
		return mVertexNormals[index];
	}

	bool RayMeshIntersection(const Ray* ray, float &tNear, int &triangleIndex, float u, float v)
	{
		int j = 0;

		// For each triangle in the mesh, check if the ray intersects
		for (int i = 0; i < mNumTriangles; ++i)
		{
			Vector3 v0 = mVertices[mIndexBuffer[j]];
			Vector3 v1 = mVertices[mIndexBuffer[j + 1]];
			Vector3 v2 = mVertices[mIndexBuffer[j + 2]];

			Triangle* triangle = new Triangle(v0, v1, v2);

			float t;
			if (RayTriangleIntersectionMT(triangle, ray, t, u, v, false) && t < tNear)
			{
				tNear = t;

				triangleIndex = i;

				return true;
			}
			j += 3;
		}

		return false;
	}
};

//bool RayTriangleIntersection
//(
//	Triangle* triangle, const Vector3 &origin, const Vector3 &direction,
//	float &t, float &u, float &v
//)
//{
//	/*
//		To start, we must first find P
//
//		P = O + tR
//		where
//		P = the pont where the ray intersects the plane of the triangle
//		O = the origin of the ray (Vecotr3 origin)
//		t = the distance from the ray's origin to the point of intersection
//		R = direction of the ray (Vector3 direction)
//	*/
//	// Calculate plane normal
//	Vector3 v0v1 = triangle->v1 - triangle->v0;
//	Vector3 v0v2 = triangle->v2 - triangle->v0;
//	Vector3 N = v0v1.Cross(v0v2);
//	//N.Normalise();
//	float denom = N.Dot(N);
//
//	// If N dot R = 0 then ray is parallel to the plane and don't intersect
//	// float therefore can't equal 0, use upper bound of floating point
//	float NDotR = N.Dot(direction); // This is used later on if ray and plane are parallel
//	if (fabs(NDotR) < kEpsilon && fabs(NDotR) > -kEpsilon)
//	{
//		return false;
//	}
//
//	/*
//		Ax + By + Cz + D = 0
//		where
//		A, B, C = coordinates of normal to the plane
//		x, y, z = coordinates of any point on the plane (v0 is a vertex of the triangle, therefore is on the plane)
//		D       = distance from the origin to the plane
//
//		The equation can be rearranged to
//		D = N(A, B, C) dot v0
//	*/
//	float D = N.Dot(triangle->v0);
//
//	/*
//		Now that we have D, we can now calculate t using the equation by substituting
//		P = O + tR
//		   into
//		Ax + By + Cz + D = 0
//
//		This rearranges to
//		t = - (N(A, B, C) dot O) + D  /  N(A, B, C) dot R
//	*/
//	t = (N.Dot(origin)) + D / NDotR;
//
//	// Check if the triangle is behind the ray
//	// If so, it shouldn't be visible
//	if (t < 0.0f)
//	{
//		return false;
//	}
//
//	// We can now complete the original equation
//	Vector3 P = origin + (t * direction);
//
//	/*
//		We must now test whether or not P is outside of the triangle. If it is not, then the ray
//		does not intersect the triangle.
//
//		Calculate one edge of the triangle using 2 vertices
//		Calculate the vector between one vertex and the point P
//		Calculate the cross product between both of these, the dot product that with the normal
//		If the point P is inside the triangle, the dot product will be positive, (to the left of the edge)
//		If the point P is outside the triangle, it will be negative, and to the right of the edge
//
//		Do this for each edge on the triangle
//	*/
//	// Edge 0
//	Vector3 edge0 = triangle->v1 - triangle->v0;
//	Vector3 vp0 = P - triangle->v0;
//	Vector3 C = edge0.Cross(vp0);
//	if (N.Dot(C) < 0)
//	{
//		return false;
//	}
//
//	// Edge 1
//	Vector3 edge1 = triangle->v2 - triangle->v1;
//	Vector3 vp1 = P - triangle->v1;
//	C = edge1.Cross(vp1);
//	if (N.Dot(C) < 0)
//	{
//		return false;
//	}
//	u = N.Dot(C);
//
//	// Edge 2
//	Vector3 edge2 = triangle->v0 - triangle->v2;
//	Vector3 vp2 = P - triangle->v2;
//	C = edge2.Cross(vp2);
//	if (N.Dot(C) < 0)
//	{
//		return false;
//	}
//	v = N.Dot(C);
//
//	u /= denom;
//	v /= denom;
//
//	// If P is within triangle, then ray and triangle intersect
//	return true;
//}

bool RayTriangleIntersectionMT
(
	Triangle* triangle, const Ray* ray,
	float &t, float &u, float &v,
	bool culling
)
{
	Vector3 edge1, edge2, tvec, pvec, qvec;
	float det, detInverse;

	// Find edges that share v0
	edge1 = triangle->v1 - triangle->v0;
	edge2 = triangle->v2 - triangle->v0;

	pvec = ray->direction.Cross(edge2);

	det = edge1.Dot(pvec);

	if (culling)
	{
		if (det < kEpsilon)
		{
			return false;
		}

		// tvec = distance from v0 to origin
		tvec = ray->origin - triangle->v0;
		u = tvec.Dot(pvec);
		if (u < 0.0f || u > det)
		{
			return false;
		}

		qvec = tvec.Cross(edge1);
		v = ray->direction.Dot(qvec);
		if (v < 0.0f || (u + v) > det)
		{
			return false;
		}

		t = edge2.Dot(qvec);
		detInverse = 1.0f / det;
		t *= detInverse;
		u *= detInverse;
		v *= detInverse;
	}
	else
	{
		if (det > -kEpsilon && det < kEpsilon)
		{
			return false;
		}
		detInverse = 1.0f / det;

		tvec = ray->origin - triangle->v0;
		
		u = tvec.Dot(pvec) * detInverse;
		if (u < 0.0f || u > 1.0f)
		{
			return false;
		}

		qvec = tvec.Cross(edge1);
		v = ray->direction.Dot(qvec) * detInverse;
		if (v < 0.0f || (u + v) > 1.0f)
		{
			return false;
		}

		t = edge2.Dot(qvec) * detInverse;
	}

	return true;
}

bool traceRay
(
	Ray* ray, Mesh* mesh,
	float &tNear, float &u, float &v,
	int& index
)
{
	if (mesh->RayMeshIntersection(ray, tNear, index, u, v))
	{
		return true;
	}

	return false;
}

Vector3 castRay
(
	Ray* ray, Mesh* mesh
)
{
	Vector3 hitColour = backgroundColour;
	float tNear, u, v;
	int index;

	if (traceRay(ray, mesh, tNear, u, v, index))
	{
		Vector3 hitPoint = tNear * ray->direction + ray->origin;
		Vector3 hitNormal = mesh->GetNormal(index);
		hitColour = Vector3(u, v, 0.0f);
	}

	return hitColour;
}

int main()
{
	// Colours (R, G, B)
	Vector3 triColours[3] = { {1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f}, {0.1f, 0.3f, 0.7f} };

	Vector3 *frameBuffer = new Vector3[kScreenWidth * kScreenHeight];
	Vector3 *pixels = frameBuffer;


	// Ray origin starts at the camera origin
	Vector3 origin(0.0f, 0.0f, 0.0f);

	Ray* newRay = new Ray();

	Vector3 v0(-1.0f, -1.0f, -3.0f); // Bottom left
	Vector3 v1(1.0f, -1.0f, -3.0f);  // Bottom right
	Vector3 v2(0.0f, 1.0f, -3.0f);   // Middle

	Triangle* triangle1 = new Triangle(v0, v1, v2);

	// For each pixel, decide what colour it should be
	for (int i = 0; i < kScreenHeight; ++i)
	{
		for (int j = 0; j < kScreenWidth; ++j)
		{
			float x = (2 * (j + 0.5f) / (float)kScreenWidth - 1) * aspectRatio;
			float y = (1 - 2 * (i + 0.5f) / (float)kScreenHeight);

			newRay->direction = Vector3(x, y, -1);
			newRay->direction.Normalise();

			float t;
			float u, v;
			if (RayTriangleIntersectionMT(triangle1, newRay, t, u, v, false))
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

	// Output to PPM file
	ofstream ofs("./output.ppm", ios::binary);
	ofs << "P6" << kScreenWidth << " " << kScreenHeight << "\n255\n";
	for (int i = 0; i < kScreenHeight * kScreenWidth; ++i) {
		char r = (char)(255 * (0, 1, frameBuffer[i].x));
		char g = (char)(255 * (0, 1, frameBuffer[i].y));
		char b = (char)(255 * (0, 1, frameBuffer[i].z));
		ofs << r << g << b;
	}
	ofs.close();

	// MEMORY DEALLOCATION //
	delete[] frameBuffer;

	return 0;
}