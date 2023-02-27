# Document Scope
The following are the steps I followed to solve the addition of a "skirt" into an STL file format. The program is interactive and an instance of it would look like this image
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
Notice how the structure `Vertex` has a **position** and a **normal**. Also, the class `CADModel` has the methods `getCenter()` and `getScaleFactor()` that I will describe in a bit

The `load` method does the work of grabbing the data from our STL file. It returns an object of the type `const aiScene *` that has the description of the model, i.e.

1. How many meshes does the model has
2. For each mesh what is the relation between the vertices and connectivity indices, a.k.a topology

and such object type should be plugin in the `processNode()` method which is the one responsible for using `processMesh()` to collect the geometrical information. The implementation for all those methods looks like
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

And voila, that will be the code to load an STL file (in practice any other format of a polygonal mesh). It is worth noticing that

1. When fetching the normal of a vertex we are normalizing it via `glm::normalize()`. For just loading the geometry that is irrelevant but I need those normals to be unit length for a later part so better to do it as soon as possible.
2. The collection of connectivity is supposing to have **triangular face** since it is collecting 3 indices per face. This could be not entirely true but Assimp tries to *protect* such via the flag `aiProcess_Triangulate`. I could also write code there to decide what happens if a face has more than 3 indices (like a quad for example) but that is out of the scope of the task.

Everything good but how do I know the dimensions of the object? Well, a common way to solve such is to find the **axis aligned bounding box**. For such the implementation is simple, just traverse all vertices and keep the max and min in each direction
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
With the `min` and `max` positions I am able to find `scaleFactor` and `center` variables. The scale factor and center is important since *most likely* the STL has its own coordinate system and in order to display in on screen I have to do some transformations (scaling and translation) to be Normalize Device Coordinates (NDC). The NDC concept and what is behind a graphics pipeline is out of the scope of the task but a good reference is https://learnopengl.com/

Since I still don't have a *graphics output* let me show you the information of previous methods via the console using `printf()`
<p align="center"><img src="./OutputImages/WindowsTerminal_CeaHphSoXv.png"></p>

## Creating Window To Display STL
I have the STL file now in memory but life is not fun if I don'tsee something on the screen. So let me show you how to render | draw the geometry I just collected into a GLFW window that uses OpenGL.

Given that the task is not to create a full renderer I will be using OpenGL Immediate Mode (a.k.a Old OpenGL or Legacy OpenGL) since it is a bit tedious to create the GPU objects (VAO,VBO) as well the shaders. Thus, no fancy lighting in the display that I will be showing neither optimization of rendering geometries.
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

# Finding the edge boundary

# CompGeo_Add_Skirt

## Objective
Modify a shape, a 3D geometry, so that it can be "formed" using our Dual Sided Incremental Forming (DSIF) process. 
<p align="center"><img src=https://user-images.githubusercontent.com/91622575/174198770-e7d3b0c8-156c-4757-8b05-1a29f07e73e5.png></p>

## Context
At Machina Labs, we use our robots to push, stretch, bend and fold flat sheets of metal into various shapes. This "forming" process results in the originally flat sheet of metal, now having a 3D shape protruding out of it. Typically, the desired "part" is only a portion of the protruding shape, and which eventually gets cut out of it. 
<p align="center"><img src=https://user-images.githubusercontent.com/91622575/174198462-d68516ca-4f35-4b5a-b535-1a319c6571ba.png></p>


To be "formable" or manufacturable using our DSIF process, the geometry must resemble a protrusion off a flat surface. In other words, all its edges or the entire perimeter of the geometry must lie in a single plane. Typically however, the parts required by our clients do not satisfy this requirement! So, in order to make it manufacturable, the "part" geometry provided to us by clients must be adapted. The process of adapting the given geometry is to first, reorient the geometry to be optimal for the following steps. Next, project edges of the geometry onto a flat surface. This creates extension surfaces with perimeter ending in a single plane. The edges thus extended with the new surfaces, fulfil the criteria of being a protrusion. However, this surface will need additional modifications as explained below. Such an extension surface is what we call a "skirt" - for obvious reasons.
<p align="center"><img src=https://user-images.githubusercontent.com/91622575/174199005-22c1a9ef-db98-4096-b50f-51f5c8893031.png></p>
<p align="center"><img src=https://user-images.githubusercontent.com/91622575/174199072-42ceecb8-fc4c-4504-9f33-cdafa2a7327d.png></p>


In addition to its perimeter lying in a single plane, an ideal manufacturable surface typically has a draft, or is slanted less than 90 degrees towards the plane. Typically the transition from part surface to extended surface, over the original edge, causes a change in direction. For manufacturability, it is best for the newly added skirt surface to extend smoothly beyond the edge of the orignal part. Often this can be achieved by orienting the part optimally before extending its edges. As such, to be manufacturable, there are several other requirements and constraints - but those are beyond the scope of this assignment.

## Task
Create a program that...
- Takes in a 3D geometry in form of a .STL file. Use provided file - part.STL.
- Programmatically adds a "skirt" to the provided 3D geometry.
- Outputs a new .STL file - part_skirted.STL.

## Constraints, Assumptions and Terminology
- There are several shapes which are impossible to form, and are assumed not to be used as input. Some examples are...
  * Shapes that form a closed volume - a sphere,
  * Shapes with intersecting planes,
  * Mobius strip, etc.
- To understand the requirements below, assume the part to be in a 3D cartesian co-ordinate space, with the part entirely in negative z space (the convention we use).
- Once the part is skirted, its new perimeter, formed by extending the edges is located in the z=0 plane.
- Imagine the new skirt to be a slanting wall falling from the z=0 plane. The acute angle this wall makes with the z=0 plane will be referred to as the "wall-angle" and will be a constrained value. 
- The edges shown in images here are filleted - another requirement of manufacturable geometries. But you can ignore this requirement and create edges without fillets.

## Requirements
- The program can be a simple command line script, or optionally include a user-friendly GUI.
- Develop using your language of choice (we prefer Python).
- Constraint: The newly added surfaces will have a wall-angle no greater than 70 degrees with respect to the z=0 plane.
- The provided geometry (part.STL) may be re-oriented using your choice of CAD tools, before being used as input to your program.
- The deliverables must include documentation and all artifacts necessary to setup and run the program.
- To enable review, deliverables must include code - not just executables. And per good programming practices, must be sufficiently documented for ease of understanding.

## Optional Challenge (any or all of the features for bonus points)
- Update your program to accept any .STL file (limited by the constraints described above) and output a skirted geometry.
- Create a GUI that includes display windows to visualize the input and output 3D geometries.
- Make the visualization window interactive, allowing user to manipulate the geometry displayed in the window.
- Allow user to re-orient and update the input geometry.
- Add user controls so the constraints can be adjusted and applied to update output geometry.

## Submission
In order to submit the assignment, do the following:

1. Navigate to GitHub's project import page: [https://github.com/new/import](https://github.com/new/import)

2. In the box titled "Your old repository's clone URL", paste the homework repository's link: [https://github.com/Machina-Labs/CompGeo_Add_Skirt](https://github.com/Machina-Labs/CompGeo_Add_Skirt)

3. In the box titled "Repository Name", add a name for your local homework (ex. `Add_skirt_soln`)

4. Set privacy level to "Public", then click "Begin Import" button at bottom of the page.

5. Develop your homework solution in the cloned repository and push it to Github when you're done. Extra points for good Git hygiene.

6. Send us the link to your repository.
