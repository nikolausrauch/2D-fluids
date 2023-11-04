#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/rotate_vector.hpp>

#ifdef __APPLE__
    #include <OpenGL/gl.h>
#else

#ifdef _WIN32
    #include <windows.h>
#endif
    #include <GL/gl.h>
#endif

#include <GLFW/glfw3.h>
