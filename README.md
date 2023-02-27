# Document Scope
The following are the steps I followed to solve the addition of a "skirt" into an STL file format. The program is interactive and a tipycal instance would look like this image
<p align="center"><img src="./OutputImages/AddSkirtToSTL_vbz8tARRTZ.png"></p>

Above image is showing a skirt generated for a simple STL file that represents a plane. The skirt has a 45 degrees angle and 10 levels of discretization

## Development Environment
I decided to use the following technologies to do the solution

1. **Microsoft Visual Studio 2019 with C++-14 Standard :** To use the Standard Template Library
2. **CMake 3.0 or higher :** To have crossplatform solution
3. **Open-Asset-Importer -Library (Assimp) :** To import | export an STL file 
4. **GLFW :** To create an OS window
5. **OpenGL Mathematics (GLM) :** To do math with vectors and matrices
6. **Dear ImGui :** To create a Graphical User Interface 

## Open STL file
Initially one should take a look at how the geometry looks. For that you can use any software that is over the internet that supports STL file format (some examples are MeshLab or Blender)

Sound easy but given that I will have to do things in code let me show you how to use the library Assimp for that purpose (based https://learnopengl.com/Model-Loading/Model)

I created 2 classes, `Mesh` and `CADModel`. A `Mesh` contains the geometry and connectivity while a `CADModel` is a composition of meshes. Here is how the interface looks for each class
<details><summary> Show Code </summary>

```Cpp
struct Vertex
{
	glm::vec3 position;
	glm::vec3 normal;
};
class Mesh
{
public:
	Mesh(const std::vector<Vertex> &, const std::vector<unsigned int> &);
	const std::vector<unsigned int> & getIndices() const;
	const std::vector<Vertex> & getVertices() const;
private:
	std::vector<unsigned int> indices;
	std::vector<Vertex> vertices;
};
class CADModel
{
public:
	void load(const std::string &);
	const Mesh & getMesh(unsigned int) const;
	unsigned int getNumberOfMeshes() const;
	glm::vec3 getCenter();
	float getScaleFactor();
private:
	Mesh processsMesh(aiMesh *);
	void processNode(aiNode *, const aiScene *);
	std::vector<Mesh> meshes;
	glm::vec3 center;
	float scaleFactor;
};
```
</details><br>

Notice how the structure `Vertex` has a **position** and a **normal**. Also, the class `CADModel` has the methods `getCenter()` and `getScaleFactor()` that I will describe in a bit

The `load` method does the work of grabbing the data from our STL file. It returns an object of the type `const aiScene *` that has the description of the model, i.e.

1. How many meshes does the model has
2. For each mesh what is the relation between the vertices and connectivity indices, a.k.a topology

and such object type should be plugin in the `processNode()` method which is the one responsible for using `processMesh()` to collect the geometrical information. The implementation for all those methods looks like
<details><summary> Show Code </summary>

```Cpp
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
void CADModel::load(const std::string & filename)
{
	Assimp::Importer importer;
	auto processingFlags = (aiProcess_JoinIdenticalVertices | aiProcess_Triangulate);
	auto scene = importer.ReadFile(filename.c_str(), processingFlags);
	if (!scene)
	{
		printf("Error loading file %s\n", filename.c_str());
		return;
	}
	meshes.clear();
	processNode(scene->mRootNode, scene);
}
```
</details><br>

And voila, that will be the code to load an STL file (in practice any other format of a polygonal mesh). It is worth noticing that

1. When fetching the normal of a vertex we are normalizing it via `glm::normalize()`. For just loading the geometry that is irrelevant but I need those normals to be unit length for a later part so better to do it as soon as possible.
2. The collection of connectivity is supposing to have **triangular face** since it is collecting 3 indices per face. This could be not entirely true but Assimp tries to *protect* such via the flag `aiProcess_Triangulate`. I could also write code there to decide what happens if a face has more than 3 indices (like a quad for example) but that is out of the scope of the task.

Everything good but how do I know the dimensions of the object? Well, a common way to solve such is to find the **axis aligned bounding box**. For such the implementation is simple, just traverse all vertices and keep the max and min in each direction
<details><summary> Show Code </summary>

```Cpp
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
	// ... previous code ...
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
}
```
</details><br>

With the `min` and `max` positions I am able to find `scaleFactor` and `center` variables. The scale factor and center is important since *most likely* the STL has its own coordinate system and in order to display in on screen I have to do some transformations (scaling and translation) to be Normalize Device Coordinates (NDC). The NDC concept and what is behind a graphics pipeline is out of the scope of the task but a good reference is https://learnopengl.com/

Since I still don't have a *graphics output* let me show you the information of previous methods via the console using `printf()`
<p align="center"><img src="./OutputImages/WindowsTerminal_CeaHphSoXv.png"></p>

## Creating Window To Display STL
I have the STL file now in memory but life is not fun if I don'tsee something on the screen. So let me show you how to render | draw the geometry I just collected into a GLFW window that uses OpenGL.

Given that the task is not to create a full renderer I will be using OpenGL Immediate Mode (a.k.a Old OpenGL or Legacy OpenGL) since it is a bit tedious to create the GPU objects (VAO,VBO) as well the shaders. Thus, no fancy lighting in the display that I will be showing neither optimization of rendering geometries.
<details><summary> Show Code </summary>

```Cpp
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

	void renderMesh(const CADModel & model)
	{
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
		glLineWidth(1.0f);
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
}
int main(int argc, char ** argv)
{
	glfwInit();
	auto window = glfwCreateWindow(800, 600, "Alejandro Guayaquil - 2023", nullptr, nullptr);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	while (!glfwWindowShouldClose(window))
	{
		int windowBufferWidth = -1;
		int windowBufferHeight = -1;
		glfwGetFramebufferSize(window, &windowBufferWidth, &windowBufferHeight);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		auto aspectRatio = float(windowBufferWidth) / float(windowBufferHeight);
		auto projectionMatrix = glm::perspective(45.0f, aspectRatio, 0.1f, 1000.0f);
		glLoadMatrixf(glm::value_ptr(projectionMatrix));
		glMatrixMode(GL_MODELVIEW);
		auto identityMatrix = glm::mat4(1.0f);
		auto Ry = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
		auto Rx = glm::rotate(identityMatrix, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
		auto Tx = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, +0.0f));
		auto Ty = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, +0.0f));
		auto Tz = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, -3.0f));
		auto viewMatrix = (Tx * Ty * Tz * Rx * Ry);
		auto modelMatrix = glm::mat4(1.0f);
		auto modelviewMatrix = glm::mat4(1.0f);
		modelMatrix = identityMatrix;
		modelviewMatrix = (viewMatrix * modelMatrix);
		glLoadIdentity();
		glLoadMatrixf(glm::value_ptr(modelviewMatrix));
		LegacyOpenGL::renderWorldAxes();
		auto S = glm::scale(identityMatrix, glm::vec3(1.0f / CADmodel.getScaleFactor()));
		modelMatrix = S;
		modelviewMatrix = (viewMatrix * modelMatrix);
		glLoadIdentity();
		glLoadMatrixf(glm::value_ptr(modelviewMatrix));
		LegacyOpenGL::renderMesh(CADmodel);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glfwDestroyWindow(window);
	glfwTerminate();
}
```
</details><br>

And thus when I run the code we can finally see a graphical output of the STL file.
<p align="center"><img src="./OutputImages/IktaTK2yeg.png"></p>

It looks about right that the STL part has a larger dimension in the -z component (the world origin coordinate axes are represented in red(x)-green(y)-blue(z)). I added also the display of the axis aligned bounding box but omitted that portion of code in previous snippet.

Notice the following
1. I am fetching the `getScaleFactor()` of the `CADModel` object since I need to rescale the geometry to *fit* in the display window.
2. I am using a common **perspective projection** with 45 degrees of field of view and keeping the aspect ratio not distorted
3. It might seem wasteful the transformation I have for the *view matrix* (`Rx`, `Ry`, `Tx`, `Ty`, `Tz`) but in the final code I added (which is not show in above snippet) a trackball camera view where you can rotate with the left button mouse, pan with the right button mouse, and zoom in | out with the scroll wheel

And just as a sanity check I decided to open other STL files to see how such geometries such look.
<p align="center"><img src="./OutputImages/34hTjVSKJv.png"></p>

Now I have a setup to start creaing the skirt of the STL file

## Finding the edge boundary

