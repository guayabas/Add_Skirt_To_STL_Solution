#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <GLFW/glfw3.h>

std::ostream & operator<<(std::ostream & os, glm::vec3 & v) 
{
	os << v.x << " " << v.y << " " << v.z;
	return os;
}

struct Vertex
{
	glm::vec3 position;
	glm::vec3 normal;
};

class Mesh
{
public:
	Mesh(const std::vector<Vertex> &, const std::vector<unsigned int> &);
	const std::vector<Vertex> & getVertices() const;
	const std::vector<unsigned int> & getIndices() const;

private:
	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
};

const std::vector<Vertex> & Mesh::getVertices() const
{
	return vertices;
}

const std::vector<unsigned int> & Mesh::getIndices() const
{
	return indices;
}

Mesh::Mesh(const std::vector<Vertex> & vertices, const std::vector<unsigned int> & indices)
	: vertices(vertices)
	, indices(indices)
{

}

struct Edge
{
	Edge(glm::vec3 v0, glm::vec3 v1, glm::vec3 n0, glm::vec3 n1)
		: v0(v0)
		, v1(v1)
		, n0(n0)
		, n1(n1)
	{

	}
	//bool boundary = true;
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 n0;
	glm::vec3 n1;
};

struct HashForPoint
{
	// https://en.cppreference.com/w/cpp/utility/hash
    std::size_t operator()(glm::vec3 const & s) const noexcept
    {
		std::string stringContent = std::to_string(s.x) + "_" + std::to_string(s.y) + "_" + std::to_string(s.z); 
		return std::hash<std::string>{}(stringContent);
    }
};

void printEdgeInformation(const Edge & edge)
{
	printf("\t(%3.5f, %3.5f, %3.5f) - (%3.5f, %3.5f, %3.5f)\n", 
		edge.v0.x, edge.v0.y, edge.v0.z,
		edge.v1.x, edge.v1.y, edge.v1.z
	);
}

struct ComparatorForPoint
{
	bool operator()(const glm::vec3 & lhs, const glm::vec3 & rhs) const
	{
		return ((lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z));
	}
};

inline std::string getStringFromPoint(const glm::vec3 & v)
{
	return 
		std::to_string(v.x) + "_" + 
		std::to_string(v.y) + "_" + 
		std::to_string(v.z); 
}

struct HashForEdge
{
	std::size_t operator()(Edge const & edge) const
	{
		std::string s0 = getStringFromPoint(edge.v0);
		std::string s1 = getStringFromPoint(edge.v1);
		return std::hash<std::string>{}(s0 + "_" + s1);
	}
};

struct ComparatorForEdge
{
	bool operator()(const Edge & lhs, const Edge & rhs) const
	{
		auto pointComparator = ComparatorForPoint();
		auto order1 = 
			pointComparator.operator()(lhs.v0, rhs.v0) && 
			pointComparator.operator()(lhs.v1, rhs.v1);
		auto order2 = 
			pointComparator.operator()(lhs.v0, rhs.v1) && 
			pointComparator.operator()(lhs.v1, rhs.v0);
		return (order1 || order2);
	}
};

struct BoundaryEdge
{
	BoundaryEdge(glm::vec3 v0, glm::vec3 v1, glm::vec3 n0, glm::vec3 n1)
		: v0(v0)
		, v1(v1)
		, n0(n0)
		, n1(n1)
	{

	}

	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 n0;
	glm::vec3 n1;
};

struct Skirt
{
	std::vector<std::vector<glm::vec3>> vertices;
	unsigned int numberOfLevels = 0;
	float lengthScale = 0.0f;
	float fallingAngle = 0.0f;

	void setNumberOfLevels(unsigned int levels)
	{
		numberOfLevels = levels;
	}

	void setFallingAngle(float angle)
	{
		fallingAngle = angle;
	}

	void setLengthScale(float length)
	{
		lengthScale = length;
	}
};

struct BoundaySegment
{
	std::unordered_set<glm::vec3, HashForPoint, ComparatorForPoint> vertices;
	std::vector<Edge> edges;
	glm::vec3 normal = glm::vec3(0.0f);
	glm::vec3 direction = glm::vec3(0.0f);
	glm::vec3 min = glm::vec3(+std::numeric_limits<float>::infinity());
	glm::vec3 max = glm::vec3(-std::numeric_limits<float>::infinity());
};

class PerformanceClock
{
public:
	inline constexpr long long getTimeInMillisecond() const
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	}
	inline void initCount()
	{
		t1 = std::chrono::high_resolution_clock::now();
	}
	inline void stopCount()
	{
		t2 = std::chrono::high_resolution_clock::now();
	}

private:
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
};

typedef std::unordered_map<Edge, bool, HashForEdge, ComparatorForEdge> EdgeHash;

class CADModel
{
public:
	void load(const std::string &);
	const Mesh & getMesh(unsigned int) const;
	unsigned int getNumberOfMeshes() const;
	const std::vector<Edge> & getBoundary() const;
	const Skirt & getSkirt(unsigned int) const;
	glm::vec3 getCenter();
	float getScaleFactor();
	std::pair<glm::vec3, glm::vec3> getMinMax() const;
	glm::vec3 getAverageNormal();

	const std::unordered_map<std::string, float> & getMetrics() const;

	void generateSkirt(unsigned int, float, float);

	std::vector<BoundaySegment> segments;

private:
	Mesh processsMesh(aiMesh *);
	void processNode(aiNode *, const aiScene *);
	std::vector<Mesh> meshes;
	
	//std::vector<Edge> edgesContainer;
	EdgeHash edgesHash;
	std::vector<Edge> boundary;

	glm::vec3 center;
	glm::vec3 min;
	glm::vec3 max;
	glm::vec3 averageNormal = glm::vec3(0.0f);
	float scaleFactor;
	Skirt skirt;
	bool computeBoundary = true;
	bool computeSegments = true;
	bool computeSkirt = true;
	std::unordered_map<std::string, float> metrics;
	PerformanceClock clock;
};

const std::unordered_map<std::string, float> & CADModel::getMetrics() const
{
	return metrics;
}

const Skirt & CADModel::getSkirt(unsigned int index) const
{
	return skirt;
}

glm::vec3 CADModel::getAverageNormal()
{
	return averageNormal;
}

std::pair<glm::vec3, glm::vec3> CADModel::getMinMax() const
{
	return std::make_pair(min, max);
}

const std::vector<Edge> & CADModel::getBoundary() const
{
	return boundary;
}

float CADModel::getScaleFactor()
{
	return scaleFactor;
}

glm::vec3 CADModel::getCenter()
{
	return center;
}

unsigned int CADModel::getNumberOfMeshes() const
{
	return (unsigned int)meshes.size();
}

const Mesh & CADModel::getMesh(unsigned int i) const
{
	return meshes[i];
}

inline aiVector3D glmVec3ToAssimpVec3(const glm::vec3 & v)
{
	return aiVector3D(v.x, v.y, v.z);
}

inline glm::vec3 AssimpVec3ToglmVec3(const aiVector3D & v)
{
	return glm::vec3(v.x, v.y, v.z);
}

Mesh CADModel::processsMesh(aiMesh * mesh)
{
	std::vector<Vertex> vertices(mesh->mNumVertices);
	for (unsigned int i = 0; i < mesh->mNumVertices; i++)
	{
		Vertex vertex;
		vertex.position = AssimpVec3ToglmVec3(mesh->mVertices[i]);
		vertex.normal = AssimpVec3ToglmVec3(mesh->mNormals[i]);
		vertex.normal = glm::normalize(vertex.normal);
		vertices[i] = vertex;
	}
	size_t numberOfIndices = ((size_t)mesh->mNumFaces * 3);
	std::vector<unsigned int> indices(numberOfIndices);
	for (size_t i = 0; i < mesh->mNumFaces; i++)
	{
		auto face = mesh->mFaces[i];
		auto faceIndex = (3 * i);
		indices[faceIndex + 0] = face.mIndices[0];
		indices[faceIndex + 1] = face.mIndices[1];
		indices[faceIndex + 2] = face.mIndices[2];
	}
	return Mesh(vertices, indices);
}

void CADModel::processNode(aiNode * node, const aiScene * scene)
{
	for (unsigned int i = 0; i < node->mNumMeshes; i++)
	{
		auto mesh = scene->mMeshes[node->mMeshes[i]];
		meshes.push_back(processsMesh(mesh));
	}
	for (unsigned int i = 0; i < node->mNumChildren; i++)
	{
		processNode(node->mChildren[i], scene);
	}
}

bool areVerticesEqual(const std::pair<glm::vec3, glm::vec3> & edge1, const std::pair<glm::vec3, glm::vec3> & edge2)
{
	auto c1 = (edge1.first.x == edge2.first.x && edge1.first.y == edge2.first.y && edge1.first.z == edge2.first.z);
	auto c2 = (edge1.second.x == edge2.second.x && edge1.second.y == edge2.second.y && edge1.second.z == edge2.second.z);
	auto c3 = (edge1.first.x == edge2.second.x && edge1.first.y == edge2.second.y && edge1.first.z == edge2.second.z);
	auto c4 = (edge1.second.x == edge2.first.x && edge1.second.y == edge2.first.y && edge1.second.z == edge2.first.z);
	auto equal = (c1 && c2) || (c3 && c4);
	return equal;
}

void addEdgeToEdgesHash(EdgeHash & edges, Edge & edge)
{
	auto dualEdge = Edge(edge.v1, edge.v0, edge.n1, edge.n0);
	if (edges.find(edge) != edges.end() || edges.find(dualEdge) != edges.end())
	{
		edges[edge] = false;
		edges[dualEdge] = false;
	}
	else
	{
		edges[edge] = true;
	}
}

void CADModel::generateSkirt(unsigned int levels, float length, float angle)
{
	for (auto & skirtSegment : skirt.vertices)
	{
		skirtSegment.clear();
	}
	skirt.vertices.clear();
	skirt.setNumberOfLevels(levels);
	skirt.setLengthScale(length);
	skirt.setFallingAngle(angle);

	// Skirt
	clock.initCount();
	float fallingDirectionLength = skirt.lengthScale * std::tanf(glm::radians(skirt.fallingAngle));
	if (computeSkirt && !segments.empty())
	{
		for (auto & segment : segments)
		{
			for (auto & boundaryedge : segment.edges)
			{
				skirt.vertices.push_back(std::vector<glm::vec3>());
				skirt.vertices.back().push_back(boundaryedge.v0);
				skirt.vertices.back().push_back(boundaryedge.v1);
				for (int n = 1; n <= skirt.numberOfLevels; n++)
				{
					float normalLength = ((n / float(skirt.numberOfLevels)) * skirt.lengthScale);
					float slantValue = fallingDirectionLength * (n / float(skirt.numberOfLevels));
					skirt.vertices.back().push_back(boundaryedge.v0 - normalLength * segment.normal + slantValue * segment.direction);
					skirt.vertices.back().push_back(boundaryedge.v1 - normalLength * segment.normal + slantValue * segment.direction);
				}
			}
		}

		// Fix ordering of max/min
		skirt.vertices.push_back(std::vector<glm::vec3>());
		for (int n = 0; n <= skirt.numberOfLevels; n++)
		{
			float normalLength = ((n / float(skirt.numberOfLevels)) * skirt.lengthScale);
			float slantValue = fallingDirectionLength * (n / float(skirt.numberOfLevels));
			skirt.vertices.back().push_back(segments[4].max - normalLength * segments[4].normal + slantValue * segments[4].direction);
			skirt.vertices.back().push_back(segments[0].min - normalLength * segments[0].normal + slantValue * segments[0].direction);
			//skirt.vertices.back().push_back(segments[0].min - normalLength * segments[0].normal + slantValue * segments[0].direction);
			//skirt.vertices.back().push_back(segments[3].max - normalLength * segments[3].normal + slantValue * segments[3].direction);
		}
		skirt.vertices.push_back(std::vector<glm::vec3>());
		for (int n = 0; n <= skirt.numberOfLevels; n++)
		{
			float normalLength = ((n / float(skirt.numberOfLevels)) * skirt.lengthScale);
			float slantValue = fallingDirectionLength * (n / float(skirt.numberOfLevels));
			skirt.vertices.back().push_back(segments[0].max - normalLength * segments[0].normal + slantValue * segments[0].direction);
			skirt.vertices.back().push_back(segments[5].max - normalLength * segments[5].normal + slantValue * segments[5].direction);
			//skirt.vertices.back().push_back(segments[0].max - normalLength * segments[0].normal + slantValue * segments[0].direction);
			//skirt.vertices.back().push_back(segments[2].min - normalLength * segments[2].normal + slantValue * segments[2].direction);
		}
		skirt.vertices.push_back(std::vector<glm::vec3>());
		for (int n = 0; n <= skirt.numberOfLevels; n++)
		{
			float normalLength = ((n / float(skirt.numberOfLevels)) * skirt.lengthScale);
			float slantValue = fallingDirectionLength * (n / float(skirt.numberOfLevels));
			skirt.vertices.back().push_back(segments[5].min - normalLength * segments[5].normal + slantValue * segments[5].direction);
			skirt.vertices.back().push_back(segments[3].max - normalLength * segments[3].normal + slantValue * segments[3].direction);
			//skirt.vertices.back().push_back(segments[2].max - normalLength * segments[2].normal + slantValue * segments[2].direction);
			//skirt.vertices.back().push_back(segments[1].min - normalLength * segments[1].normal + slantValue * segments[1].direction);
		}
		skirt.vertices.push_back(std::vector<glm::vec3>());
		for (int n = 0; n <= skirt.numberOfLevels; n++)
		{
			float normalLength = ((n / float(skirt.numberOfLevels)) * skirt.lengthScale);
			float slantValue = fallingDirectionLength * (n / float(skirt.numberOfLevels));
			skirt.vertices.back().push_back(segments[3].min - normalLength * segments[3].normal + slantValue * segments[3].direction);
			skirt.vertices.back().push_back(segments[4].min - normalLength * segments[4].normal + slantValue * segments[4].direction);
			//skirt.vertices.back().push_back(segments[1].max - normalLength * segments[1].normal + slantValue * segments[1].direction);
			//skirt.vertices.back().push_back(segments[3].min - normalLength * segments[3].normal + slantValue * segments[3].direction);
		}
	}
	clock.stopCount();
	metrics["skirt"] = clock.getTimeInMillisecond();
}

void checkMinMaxVertexPosition(const Vertex & vertex, glm::vec3 & min, glm::vec3 & max)
{
	if (vertex.position.x < min.x) min.x = vertex.position.x;
	if (vertex.position.y < min.y) min.y = vertex.position.y;
	if (vertex.position.z < min.z) min.z = vertex.position.z;
	if (vertex.position.x > max.x) max.x = vertex.position.x;
	if (vertex.position.y > max.y) max.y = vertex.position.y;
	if (vertex.position.z > max.z) max.z = vertex.position.z;
}

void CADModel::load(const std::string & filename)
{
	Assimp::Importer importer;

	auto processingFlags = (aiProcess_JoinIdenticalVertices | aiProcess_Triangulate);
	clock.initCount();
	auto scene = importer.ReadFile(filename.c_str(), processingFlags);
	clock.stopCount();
	if (!scene)
	{
		printf("Error loading file %s\n", filename.c_str());
		return;
	}
	metrics["scene"] = clock.getTimeInMillisecond();
	meshes.clear();
	clock.initCount();
	processNode(scene->mRootNode, scene);
	clock.stopCount();

	// Find axis-aligned bounding box
	min = glm::vec3(+std::numeric_limits<float>::infinity());
	max = glm::vec3(-std::numeric_limits<float>::infinity());
	for (auto & mesh : meshes)
	{
		auto & vertices = mesh.getVertices();
		for (auto & vertex : vertices)
		{
			checkMinMaxVertexPosition(vertex, min, max);
		}
	}
	glm::vec3 dimensions = glm::vec3((max.x - min.x), (max.y - min.y), (max.z - min.z));
	scaleFactor = std::max(dimensions.x, std::max(dimensions.y, dimensions.z));
	center = (max + min) * 0.5f;
	printf("Scene %s imported in %llu ms\n", filename.c_str(), clock.getTimeInMillisecond());
	printf("\tMeshes: %d\n", scene->mNumMeshes);
	printf("\tMeshes processed in %llu ms\n", clock.getTimeInMillisecond());
	for (unsigned int i = 0; i < scene->mNumMeshes; i++)
	{
		printf("\t\t Mesh %u has %lu vertices and %lu faces\n", i, scene->mMeshes[i]->mNumVertices, scene->mMeshes[i]->mNumFaces);
	}
	printf("\tMin: (%f, %f, %f)\n", min.x, min.y, min.z);
	printf("\tMax: (%f, %f, %f)\n", max.x, max.y, max.z);
	printf("\tMax Dimension: %f\n", scaleFactor);
	printf("\tCenter: (%f, %f, %f)\n", center.x, center.y, center.z);

	// Boundary
	clock.initCount();
	if (computeBoundary)
	{
		auto indices = meshes[0].getIndices();
		auto vertices = meshes[0].getVertices();
		for (size_t i = 0; i < indices.size(); i += 3)
		{
			size_t i0 = (i + 0);
			size_t i1 = (i + 1);
			size_t i2 = (i + 2);
			glm::vec3 v0 = vertices[indices[i0]].position;
			glm::vec3 n0 = vertices[indices[i0]].normal;
			glm::vec3 v1 = vertices[indices[i1]].position;
			glm::vec3 n1 = vertices[indices[i1]].normal;
			glm::vec3 v2 = vertices[indices[i2]].position;
			glm::vec3 n2 = vertices[indices[i2]].normal;
			addEdgeToEdgesHash(edgesHash, Edge(v0, v1, n0, n1));
			addEdgeToEdgesHash(edgesHash, Edge(v1, v2, n1, n2));
			addEdgeToEdgesHash(edgesHash, Edge(v2, v0, n2, n0));
		}
		for (auto & edge : edgesHash)
		{
			if (edge.second)
			{
				boundary.push_back(edge.first);
			}
		}
		printf("Boundary has %zu edges\n", boundary.size());
	}
	clock.stopCount();
	metrics["boundary"] = clock.getTimeInMillisecond();

	// Segments
	clock.initCount();
	if (computeSegments)
	{
		segments.resize(6);
		for (size_t n = 0; n < boundary.size(); n++)
		{
			auto & edge = boundary[n];
			auto & normal = glm::normalize(edge.n0 + edge.n1);
			glm::vec3 direction = glm::cross(glm::normalize(edge.v1 - edge.v0), normal);
			int segmentIndex = -1;
			auto maxDirectionValue = std::max(std::abs(direction.x), std::max(std::abs(direction.y), std::abs(direction.z)));
			if (maxDirectionValue == direction.x)
			{
				segmentIndex = 0;
			}
			else if (maxDirectionValue == -direction.x)
			{
				segmentIndex = 1;
			}
			else if (maxDirectionValue == direction.y)
			{
				segmentIndex = 2;
			}
			else if (maxDirectionValue == -direction.y)
			{
				segmentIndex = 3;
			}
			else if (maxDirectionValue == direction.z)
			{
				segmentIndex = 4;
			}
			else
			{
				segmentIndex = 5;
			}
			segments[segmentIndex].edges.push_back(edge);
			segments[segmentIndex].vertices.insert(edge.v0);
			segments[segmentIndex].vertices.insert(edge.v1);
		}
		//printf("Number of edges in +x segment %zu\n", segments[0].edges.size());
		//printf("Number of edges in -x segment %zu\n", segments[1].edges.size());
		//printf("Number of edges in +y segment %zu\n", segments[2].edges.size());
		//printf("Number of edges in -y segment %zu\n", segments[3].edges.size());
		//printf("Number of edges in +z segment %zu\n", segments[4].edges.size());
		//printf("Number of edges in -z segment %zu\n", segments[5].edges.size());
		for (auto & segment : segments)
		{
			float minDistance = std::numeric_limits<float>::max();
			float maxDistance = std::numeric_limits<float>::min();
			for (auto & vertex : segment.vertices)
			{
				float distanceToOrigin = glm::dot(vertex, vertex);
				if (distanceToOrigin < minDistance)
				{
					minDistance = distanceToOrigin;
					segment.min = vertex;
				}
				if (distanceToOrigin > maxDistance)
				{
					maxDistance = distanceToOrigin;
					segment.max = vertex;
				}
				// small hack ... 
				if (distanceToOrigin == maxDistance && distanceToOrigin == minDistance)
				{
					ComparatorForPoint pointComparator;
					if (!pointComparator.operator()(segment.max, vertex))
					{
						segment.max = vertex;
					}
				}
			}
		}

		//for (int n = 0; n < 6; n++)
		//{
		//	printf("%d : (%3.5f, %3.5f, %3.5f) - (%3.5f, %3.5f, %3.5f)\n", n,
		//		segments[n].min.x, segments[n].min.y, segments[n].min.z,
		//		segments[n].max.x, segments[n].max.y, segments[n].max.z
		//	);
		//}

		segments[0].direction = glm::vec3(+1.0f, 0.0f, 0.0f);
		segments[1].direction = glm::vec3(-1.0f, 0.0f, 0.0f);
		segments[2].direction = glm::vec3(0.0f, +1.0f, 0.0f);
		segments[3].direction = glm::vec3(0.0f, -1.0f, 0.0f);
		segments[4].direction = glm::vec3(0.0f, 0.0f, +1.0f);
		segments[5].direction = glm::vec3(0.0f, 0.0f, -1.0f);
		for (auto & segment : segments)
		{
			glm::vec3 averageNormal = glm::vec3(0.0f);
			for (auto & boundaryedge : segment.edges)
			{
				averageNormal += glm::normalize(boundaryedge.n0 + boundaryedge.n1);
			}
			segment.normal = glm::normalize(averageNormal);
		}
	}
	clock.stopCount();
	metrics["segments"] = clock.getTimeInMillisecond();

	//Assimp::Exporter exporter;
	//aiScene outputScene;
	//outputScene.mRootNode = new aiNode();
	//outputScene.mMaterials = new aiMaterial*[1];
	//outputScene.mNumMaterials = 1;
	//outputScene.mMaterials[0] = new aiMaterial();
	//auto numMeshesInSkirt = skirt.vertices.size();
	//outputScene.mMeshes = new aiMesh*[numMeshesInSkirt];
	//outputScene.mNumMeshes = numMeshesInSkirt;
	//outputScene.mRootNode->mMeshes = new unsigned int[numMeshesInSkirt];
	//outputScene.mRootNode->mNumMeshes = numMeshesInSkirt;
	//for (size_t m = 0; m < numMeshesInSkirt; m++)
	//{
	//	outputScene.mMeshes[m] = new aiMesh();
	//	outputScene.mMeshes[m]->mMaterialIndex = 0;
	//	outputScene.mRootNode->mMeshes[m] = m;
	//	auto pMesh = outputScene.mMeshes[m];
	//	auto numVerticesInSkirt = skirt.vertices[m].size();
	//	pMesh->mVertices = new aiVector3D[numVerticesInSkirt];
	//	pMesh->mNumVertices = numVerticesInSkirt;
	//	for (size_t n = 0; n < numVerticesInSkirt; n++)
	//	{
	//		pMesh->mVertices[n] = glmVec3ToAssimpVec3(skirt.vertices[m][n]);
	//	}
	//	auto numFacesInSkirt = skirt.vertices[m].size() - 2;
	//	pMesh->mFaces = new aiFace[numFacesInSkirt];
	//	pMesh->mNumFaces = numFacesInSkirt;
	//	for (size_t n = 0; n < numFacesInSkirt; n++)
	//	{
	//		aiFace & face = pMesh->mFaces[n];
	//		face.mIndices = new unsigned int[3];
	//		face.mNumIndices = 3;
	//		face.mIndices[0] = (n + 0);
	//		if ((n % 2) == 0)
	//		{
	//			face.mIndices[1] = (n + 1);
	//			face.mIndices[2] = (n + 2);
	//		}
	//		else
	//		{
	//			face.mIndices[1] = (n + 2);
	//			face.mIndices[2] = (n + 1);
	//		}
	//	}
	//}
	//printf("Saving output.stl\n");
	//auto exportResult = exporter.Export(&outputScene, "stl", "./output.stl");
	//if (exportResult != aiReturn_SUCCESS)
	//{
	//	printf("Error exporting output.stl\n");
	//}
}

struct ApplicationParams
{
	void resetValues();
	float cameraTranslateX = 0.0f;
	float cameraTranslateY = 0.0f;
	float cameraTranslateZ = 0.0f;
	float cameraAngleAroundX = 0.0f;
	float cameraAngleAroundY = 0.0f;
	float cursorWindowPositionX = 0.0f;
	float cursorWindowPositionY = 0.0f;
	bool leftMouseButtonPressed = false;
	bool rightMouseButtonPressed = false;
	int windowFullScreenWidth = 0;
	int windowFullScreenHeight = 0;
	int windowScreenWidth = 0;
	int windowScreenHeight = 0;
	int windowScreenPositionX = 0;
	int windowScreenPositionY = 0;
	bool windowFullScreen = false;
	bool openglWireframeModeEnabled = false;
	bool openglMultisampleEnabled = false;
	bool drawWorldAxes = false;
	bool drawingShowMesh = true;
	bool drawingShowAxisAlignedBoundingBox = false;
	bool drawingShowBoundary = false;
	bool drawingShowSegments = false;
	bool drawingShowSkirt = true;
	bool showGUIBottomMenu = true;
	int skirtLevels = 1;
	float skirtLength = 1.0f;
	float skirtFallingAngle = 45.0f;
};

void ApplicationParams::resetValues()
{
	cameraTranslateX = +0.0f;
	cameraTranslateY = +0.0f;
	cameraTranslateZ = -1.0f;
	cameraAngleAroundX = 0.0f;
	cameraAngleAroundY = 0.0f;
}

ApplicationParams params;

static void keyCallback(GLFWwindow * window, int key, int, int action, int)
{
	if (action == GLFW_PRESS)
	{
		if (key == GLFW_KEY_Q)
		{
			glfwSetWindowShouldClose(window, true);
		}
		if (key == GLFW_KEY_R)
		{
			params.resetValues();
		}
		if (key == GLFW_KEY_F)
		{
			params.windowFullScreen ^= 1;
			if (params.windowFullScreen)
			{
				glfwSetWindowPos(window, 0, 0);
				glfwSetWindowAttrib(window, GLFW_DECORATED, GLFW_FALSE);
				glfwSetWindowSize(window, params.windowFullScreenWidth, params.windowFullScreenHeight);
			}
			else
			{
				glfwSetWindowPos(window, params.windowScreenPositionX, params.windowScreenPositionY);
				glfwSetWindowAttrib(window, GLFW_DECORATED, GLFW_TRUE);
				glfwSetWindowSize(window, params.windowScreenWidth, params.windowScreenHeight);
			}
		}
	}
}

static void scrollCallback(GLFWwindow * window, double xoffset, double yoffset)
{
	params.cameraTranslateZ += (0.1f * (float)yoffset);
	if (params.cameraTranslateZ > -0.1f)
	{
		params.cameraTranslateZ -= (0.1f * (float)yoffset);
	}
	//printf("%f\n", params.cameraTranslateZ);
}

static void resizeCallback(GLFWwindow * window, int width, int height)
{
    glViewport(0, 0, width, height);
}

static void cursorPositionCallback(GLFWwindow * window, double xpos, double ypos)
{
	if (ImGui::IsAnyItemActive()) return;
	float dx = (params.cursorWindowPositionX - (float)xpos);
	float dy = (params.cursorWindowPositionY - (float)ypos);
	if (params.leftMouseButtonPressed)
	{
		params.cameraAngleAroundY += dx;
		if (params.cameraAngleAroundY > +360.0f) params.cameraAngleAroundY = 0.0f;
		if (params.cameraAngleAroundY < -360.0f) params.cameraAngleAroundY = 0.0f;
		params.cameraAngleAroundX -= dy;
		if (params.cameraAngleAroundX > +360.0f) params.cameraAngleAroundX = 0.0f;
		if (params.cameraAngleAroundX < -360.0f) params.cameraAngleAroundX = 0.0f;
	}
	if (params.rightMouseButtonPressed)
	{
		params.cameraTranslateX -= 0.01f * dx;
		params.cameraTranslateY += 0.01f * dy;
		//printf("(%f, %f)\n", params.cameraTranslateX, params.cameraTranslateY);
	}
	params.cursorWindowPositionX = (float)xpos;
	params.cursorWindowPositionY = (float)ypos;
}

static void mouseButtonCallback(GLFWwindow * window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		if (!ImGui::IsAnyItemActive())
		{
			double cursorPositionX, cursorPositionY;
			glfwGetCursorPos(window, &cursorPositionX, &cursorPositionY);
			params.cursorWindowPositionX = (float)cursorPositionX;
			params.cursorWindowPositionY = (float)cursorPositionY;
			if (button == GLFW_MOUSE_BUTTON_LEFT)
			{
				params.leftMouseButtonPressed = true;
			}
			if (button == GLFW_MOUSE_BUTTON_RIGHT)
			{
				params.rightMouseButtonPressed = true;
			}
		}
	}
	if (action == GLFW_RELEASE)
	{
		params.leftMouseButtonPressed = false;
		params.rightMouseButtonPressed = false;
	}    
}

glm::vec3 colors[] = 
{
	glm::vec3(1.0f, 0.0f, 0.0f), // red
	glm::vec3(0.0f, 1.0f, 0.0f), // green
	glm::vec3(0.0f, 0.0f, 1.0f), // blue
	glm::vec3(1.0f, 1.0f, 0.0f), // yellow 
	glm::vec3(1.0f, 0.0f, 1.0f), // magenta
	glm::vec3(0.0f, 1.0f, 1.0f), // cyan
};

namespace LegacyOpenGL
{
	void renderWorldAxes()
	{
		glLineWidth(3.0f);
		glBegin(GL_LINES);
		glColor3fv(glm::value_ptr(glm::vec3(1.0f, 0.0f, 0.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(0.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(1.0f, 0.0f, 0.0f)));
		glColor3fv(glm::value_ptr(glm::vec3(0.0f, 1.0f, 0.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(0.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(0.0f, 1.0f, 0.0f)));
		glColor3fv(glm::value_ptr(glm::vec3(0.0f, 0.0f, 1.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(0.0f)));
		glVertex3fv(glm::value_ptr(glm::vec3(0.0f, 0.0f, 1.0f)));
		glEnd();
		glLineWidth(1.0f);
	}

	void renderSegments(const CADModel & model)
	{
		glLineWidth(5.0f);
		glBegin(GL_LINES);
		for (int m = 0; m < model.segments.size(); m++)
		{
			glColor3fv(glm::value_ptr(colors[m]));
			for (int n = 0; n < model.segments[m].edges.size(); n++)
			{
				glVertex3fv(glm::value_ptr(model.segments[m].edges[n].v0));
				glVertex3fv(glm::value_ptr(model.segments[m].edges[n].v1));
			}
		}
		glEnd();
		glLineWidth(1.0f);
	}
	
	void renderSkirt(const CADModel & model)
	{
		// Skirt (has to be per mesh later)
		auto & skirt = model.getSkirt(0);
		glBegin(GL_TRIANGLES);
		for (size_t ribbonIndex = 0; ribbonIndex < skirt.vertices.size(); ribbonIndex++)
		{
			for (size_t i = 0; i < skirt.vertices[ribbonIndex].size() - 2; i += 2)
			{
				auto levelIndex = (((i / 2) + 1) / (float)skirt.numberOfLevels);
				//printf("%f\n", levelIndex);
				//glColor3f(levelIndex, 0.32f, 0.28f);
				float color = 0.9f * (1.0f - levelIndex);
				glColor3f(color, color, color);
				//glColor3f(0.0f, 0.0f, 0.0f);
				//glColor3f(levelIndex, levelIndex, levelIndex);
				glm::vec3 v0 = skirt.vertices[ribbonIndex][i + 0];
				glm::vec3 v1 = skirt.vertices[ribbonIndex][i + 1];
				glm::vec3 v2 = skirt.vertices[ribbonIndex][i + 2];
				glm::vec3 v3 = skirt.vertices[ribbonIndex][i + 3];
				//std::cout << v0;
				//std::cout << v1;
				//std::cout << v2;
				//std::cout << v3;
				glVertex3fv(glm::value_ptr(v0));
				glVertex3fv(glm::value_ptr(v1));
				glVertex3fv(glm::value_ptr(v2));
				glVertex3fv(glm::value_ptr(v1));
				glVertex3fv(glm::value_ptr(v3));
				glVertex3fv(glm::value_ptr(v2));
			}
		}
		glEnd();
		//std::cin.get();
	}

	void renderMesh(const CADModel & model)
	{
		// Mesh - triangles
		auto numberOfMeshes = model.getNumberOfMeshes();
		glColor3f(0.5f, 0.5f, 0.5f);
		glBegin(GL_TRIANGLES);
		for (unsigned int n = 0; n < numberOfMeshes; ++n)
		{
			auto & mesh = model.getMesh(n);
			auto & indices = mesh.getIndices();
			auto & vertices = mesh.getVertices();
			for (size_t i = 0; i < indices.size(); i += 3)
			{
				glm::vec3 v0 = vertices[indices[i + 0]].position;
				glm::vec3 v1 = vertices[indices[i + 1]].position;
				glm::vec3 v2 = vertices[indices[i + 2]].position;
				glVertex3fv(glm::value_ptr(v0));
				glVertex3fv(glm::value_ptr(v1));
				glVertex3fv(glm::value_ptr(v2));
			}
		}
		glEnd();

		// Mesh - lines
		glLineWidth(1.0f);
		//glPolygonOffset(1.0f, 1.0f);
		glColor3f(0.0f, 0.0f, 0.0f);
		glBegin(GL_LINES);
		for (unsigned int n = 0; n < numberOfMeshes; ++n)
		{
			auto & mesh = model.getMesh(n);
			auto & indices = mesh.getIndices();
			auto & vertices = mesh.getVertices();
			for (size_t i = 0; i < indices.size(); i += 3)
			{
				glm::vec3 v0 = vertices[indices[i + 0]].position;
				glm::vec3 v1 = vertices[indices[i + 1]].position;
				glm::vec3 v2 = vertices[indices[i + 2]].position;
				glVertex3fv(glm::value_ptr(v0));
				glVertex3fv(glm::value_ptr(v1));
				glVertex3fv(glm::value_ptr(v1));
				glVertex3fv(glm::value_ptr(v2));
				glVertex3fv(glm::value_ptr(v2));
				glVertex3fv(glm::value_ptr(v0));
			}
		}
		glEnd();
		glLineWidth(1.0f);
	}

	void renderBoundary(const CADModel & model)
	{
		glLineWidth(5.0f);
		glColor3f(0.63f, 0.12f, 0.94f);
		glBegin(GL_LINES);
		auto & boundary = model.getBoundary();
		for (auto & edge : boundary)
		{
			glVertex3fv(glm::value_ptr(edge.v0));
			glVertex3fv(glm::value_ptr(edge.v1));
		}
		glEnd();
		glLineWidth(1.0f);
	}

	void renderAxisAlignedBoundingBox(const CADModel & model)
	{
		// Axis aligned bounding box
		glLineWidth(5.0f);
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		auto minmax = model.getMinMax();
		glm::vec3 v0 = minmax.first;
		glm::vec3 v1 = glm::vec3(minmax.second.x, minmax.first.y, minmax.first.z);
		glm::vec3 v2 = glm::vec3(minmax.second.x, minmax.second.y, minmax.first.z);
		glm::vec3 v3 = glm::vec3(minmax.first.x, minmax.second.y, minmax.first.z);
		glm::vec3 v4 = glm::vec3(minmax.first.x, minmax.first.y, minmax.second.z);
		glm::vec3 v5 = glm::vec3(minmax.second.x, minmax.first.y, minmax.second.z);
		glm::vec3 v6 = minmax.second;
		glm::vec3 v7 = glm::vec3(minmax.first.x, minmax.second.y, minmax.second.z);
		glVertex3fv(glm::value_ptr(v0));
		glVertex3fv(glm::value_ptr(v1));
		glVertex3fv(glm::value_ptr(v1));
		glVertex3fv(glm::value_ptr(v2));
		glVertex3fv(glm::value_ptr(v2));
		glVertex3fv(glm::value_ptr(v3));
		glVertex3fv(glm::value_ptr(v3));
		glVertex3fv(glm::value_ptr(v0));
		glVertex3fv(glm::value_ptr(v4));
		glVertex3fv(glm::value_ptr(v5));
		glVertex3fv(glm::value_ptr(v5));
		glVertex3fv(glm::value_ptr(v6));
		glVertex3fv(glm::value_ptr(v6));
		glVertex3fv(glm::value_ptr(v7));
		glVertex3fv(glm::value_ptr(v7));
		glVertex3fv(glm::value_ptr(v4));
		glVertex3fv(glm::value_ptr(v0));
		glVertex3fv(glm::value_ptr(v4));
		glVertex3fv(glm::value_ptr(v3));
		glVertex3fv(glm::value_ptr(v7));
		glVertex3fv(glm::value_ptr(v1));
		glVertex3fv(glm::value_ptr(v5));
		glVertex3fv(glm::value_ptr(v2));
		glVertex3fv(glm::value_ptr(v6));
		glEnd();
		glLineWidth(1.0f);
	}
}

static void createWindowIcon(std::unique_ptr<GLFWimage> & icon, unsigned int width, unsigned int height)
{
	icon->width = width;
	icon->height = height;
	auto size = (width * height * 4);
	icon->pixels = new unsigned char[size];
	for (unsigned int j = 0; j < height; ++j)
	{
		for (unsigned int i = 0; i < width; ++i)
		{
			unsigned int n = ((i + j * width) * 4);
			icon->pixels[n + 0] = 0;
			icon->pixels[n + 1] = 0;
			icon->pixels[n + 2] = 0;
			icon->pixels[n + 3] = 255;
		}
	}
}

static void destroyWindowIcon(std::unique_ptr<GLFWimage> & icon)
{
	delete[] icon->pixels;
}

int main(int argc, char ** argv)
{
	CADModel CADmodel;
	//CADmodel.load("../Data/cube.stl");
	//CADmodel.load("../Data/bunny.stl");
	CADmodel.load("../Data/part.stl");
	//CADmodel.load("../Data/plane.stl");
	//CADmodel.load("../Data/disk.stl");
	CADmodel.generateSkirt(params.skirtLevels, params.skirtLength, params.skirtFallingAngle);

	glfwInit();
	auto mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	params.windowFullScreenWidth = mode->width;
	params.windowFullScreenHeight = mode->height;
	params.windowScreenWidth = 1300;
	params.windowScreenHeight = 800;
	params.windowScreenPositionX = 50;
	params.windowScreenPositionY = 50;
	auto window = glfwCreateWindow(params.windowScreenWidth, params.windowScreenHeight, "Alejandro Guayaquil - 2023", nullptr, nullptr);
	glfwSetWindowPos(window, 50, 50);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	constexpr int windowIconWidth = 3;
	constexpr int windowIconHeight = 3;
	std::unique_ptr<GLFWimage> windowIcon = std::make_unique<GLFWimage>();
	createWindowIcon(windowIcon, windowIconWidth, windowIconHeight);
	glfwSetWindowIcon(window, 1, windowIcon.get());
	glfwSetKeyCallback(window, keyCallback);
	glfwSetScrollCallback(window, scrollCallback);
	glfwSetFramebufferSizeCallback(window, resizeCallback);
	glfwSetCursorPosCallback(window, cursorPositionCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO & io = ImGui::GetIO();
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL2_Init();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	params.resetValues();
	while (!glfwWindowShouldClose(window))
	{
		if (!params.windowFullScreen) glfwGetWindowPos(window, &params.windowScreenPositionX, &params.windowScreenPositionY);
		int windowBufferWidth = -1;
		int windowBufferHeight = -1;
		glfwGetFramebufferSize(window, &windowBufferWidth, &windowBufferHeight);
		(params.openglWireframeModeEnabled) ? glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) : glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if(params.openglMultisampleEnabled)
		{
			glEnable(GL_POLYGON_SMOOTH);
			glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_POINT_SMOOTH);
		}
		else
		{
			glDisable(GL_POLYGON_SMOOTH);
			glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_POINT_SMOOTH);
		}
		if (windowBufferWidth > 0 && windowBufferHeight > 0)
		{
			ImGui_ImplOpenGL2_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();
			if (ImGui::BeginMainMenuBar())
			{
				auto topMainMenuContent = ("Framebuffer: " + 
					std::to_string(windowBufferWidth) 
					+ " x " + 
					std::to_string(windowBufferHeight)
					+ ", Framerate: " + 
					std::to_string(ImGui::GetIO().Framerate)
					);
				ImGui::Text(topMainMenuContent.c_str());
		        ImGui::EndMainMenuBar();
			}
			if (params.showGUIBottomMenu)
			{
				static bool bottomMenuOpened = true;
				ImGui::SetNextWindowPos(ImVec2(0, windowBufferHeight - 28));
				ImGui::SetNextWindowSize(ImVec2(windowBufferWidth, 28));
				ImGui::Begin("##BottomMainMenu", &bottomMenuOpened, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoTitleBar);
				auto & CADmodelsMetrics = CADmodel.getMetrics();
				auto bottomMainMenuContent = ("Scene: " + 
					std::to_string(CADmodelsMetrics.at("scene")) + " ms, Boundary: " + 
					std::to_string(CADmodelsMetrics.at("boundary")) + " ms, Segments: " + 
					std::to_string(CADmodelsMetrics.at("segments")) + " ms, Skirt: " + 
					std::to_string(CADmodelsMetrics.at("skirt")) + " ms");
				ImGui::Text(bottomMainMenuContent.c_str());
				ImGui::End();
			}
			ImGui::Begin("Program Options");
			if (ImGui::CollapsingHeader("ImGUI Settings", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::Checkbox("Show Bottom Menu", &params.showGUIBottomMenu);
			}
			if (ImGui::CollapsingHeader("OpenGL Settings", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::Checkbox("Wireframe Mode", &params.openglWireframeModeEnabled);
				ImGui::Checkbox("Multisample", &params.openglMultisampleEnabled);
			}
			if (ImGui::CollapsingHeader("Skirt Generation Settings", ImGuiTreeNodeFlags_DefaultOpen))
			{
				auto levelsChanged = ImGui::SliderInt("Number of Levels", &params.skirtLevels, 1, 25);
				auto lengthChanged = ImGui::SliderFloat("Length", &params.skirtLength, 0.1f, 10.0f);
				auto fallingAngleChanged = ImGui::SliderFloat("Falling Angle", &params.skirtFallingAngle, 0.0f, 90.0f);
				if (levelsChanged || lengthChanged || fallingAngleChanged)
				{
					CADmodel.generateSkirt(params.skirtLevels, params.skirtLength, params.skirtFallingAngle);
				}
			}
			if (ImGui::CollapsingHeader("Drawing", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::Checkbox("Show World Axes", &params.drawWorldAxes);
				ImGui::Checkbox("Show Mesh", &params.drawingShowMesh);
				ImGui::Checkbox("Show Axis Aligned Bounding Box", &params.drawingShowAxisAlignedBoundingBox);
				ImGui::Checkbox("Show Boundary", &params.drawingShowBoundary);
				ImGui::Checkbox("Show Segments", &params.drawingShowSegments);
				ImGui::Checkbox("Show Skirt", &params.drawingShowSkirt);
			}
			ImGui::End();
			ImGui::Render();
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			auto aspectRatio = float(windowBufferWidth) / float(windowBufferHeight);
			auto projectionMatrix = glm::perspective(45.0f, aspectRatio, 0.1f, 1000.0f);
			glLoadMatrixf(glm::value_ptr(projectionMatrix));
			glMatrixMode(GL_MODELVIEW);
			auto identityMatrix = glm::mat4(1.0f);
			auto Ry = glm::rotate(identityMatrix, glm::radians(params.cameraAngleAroundY), glm::vec3(0.0f, 1.0f, 0.0f));
			auto Rx = glm::rotate(identityMatrix, glm::radians(params.cameraAngleAroundX), glm::vec3(1.0f, 0.0f, 0.0f));
			auto Tx = glm::translate(identityMatrix, glm::vec3(params.cameraTranslateX, 0.0f, 0.0f));
			auto Ty = glm::translate(identityMatrix, glm::vec3(0.0f, params.cameraTranslateY, 0.0f));
			auto Tz = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, params.cameraTranslateZ));
			auto viewMatrix = (Tx * Ty * Tz * Rx * Ry);
			//auto viewMatrix = (Tx * Ty * Tz);
			auto modelMatrix = glm::mat4(1.0f);
			auto modelviewMatrix = glm::mat4(1.0f);
			modelMatrix = identityMatrix;
			modelviewMatrix = (viewMatrix * modelMatrix);
			glLoadIdentity();
			glLoadMatrixf(glm::value_ptr(modelviewMatrix));
			if (params.drawWorldAxes) LegacyOpenGL::renderWorldAxes();
			auto S = glm::scale(identityMatrix, glm::vec3(1.0f / CADmodel.getScaleFactor()));
			auto T = glm::translate(identityMatrix, glm::vec3(-CADmodel.getCenter()));
			auto R = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
			//modelMatrix = (S * T);
			modelMatrix = (R * S * T);
			//modelMatrix = S;
			modelviewMatrix = (viewMatrix * modelMatrix);
			glLoadIdentity();
			glLoadMatrixf(glm::value_ptr(modelviewMatrix));
			if (params.drawingShowBoundary) LegacyOpenGL::renderBoundary(CADmodel);
			if (params.drawingShowSegments) LegacyOpenGL::renderSegments(CADmodel);
			if (params.drawingShowMesh) LegacyOpenGL::renderMesh(CADmodel);
			if (params.drawingShowSkirt) LegacyOpenGL::renderSkirt(CADmodel);
			if (params.drawingShowAxisAlignedBoundingBox) LegacyOpenGL::renderAxisAlignedBoundingBox(CADmodel);
			ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
		}
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	destroyWindowIcon(windowIcon);
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}