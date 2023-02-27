cmake_minimum_required(VERSION 3.0)
project(ComputationGeometry)
add_executable(AddSkirtToSTL 
Main.cpp
Libraries/imgui/imgui_demo.cpp
Libraries/imgui/imgui_draw.cpp
Libraries/imgui/imgui_tables.cpp
Libraries/imgui/imgui_widgets.cpp
Libraries/imgui/imgui.cpp
Libraries/imgui/imgui.h
Libraries/imgui/imgui_internal.h
Libraries/imgui/imconfig.h
Libraries/imgui/backends/imgui_impl_glfw.cpp
Libraries/imgui/backends/imgui_impl_glfw.h
Libraries/imgui/backends/imgui_impl_opengl2.cpp
Libraries/imgui/backends/imgui_impl_opengl2.h
)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/glm)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/imgui)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/imgui/backends)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/glfw/include)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/assimp/include)
target_include_directories(AddSkirtToSTL PUBLIC Libraries/assimp/build/include)
target_link_libraries(AddSkirtToSTL 
debug 
../Libraries/glfw/build/src/Debug/glfw3 
../Libraries/assimp/build/lib/Debug/assimp-vc142-mtd
optimized
../Libraries/glfw/build/src/Release/glfw3
../Libraries/assimp/build/lib/Release/assimp-vc142-mt
opengl32
)
configure_file(
    Libraries/assimp/build/bin/Debug/assimp-vc142-mtd.dll 
    assimp-vc142-mtd.dll 
    COPYONLY
)
configure_file(
    Libraries/assimp/build/bin/Release/assimp-vc142-mt.dll 
    assimp-vc142-mt.dll 
    COPYONLY
)