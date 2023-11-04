#include "window.h"
#include "keyboard.h"
#include "mouse.h"

#include <cassert>
#include <sstream>
#include <iostream>

namespace detail
{

std::string GLFWerrorCodeToString(int code)
{
    switch (code)
    {
    case GLFW_NOT_INITIALIZED:		return "GLFW has not been initialized";
    case GLFW_NO_CURRENT_CONTEXT:	return "No context is current for this thread";
    case GLFW_INVALID_ENUM:			return "One of the arguments to the function was an invalid enum value";
    case GLFW_INVALID_VALUE:		return "One of the arguments to the function was an invalid value";
    case GLFW_OUT_OF_MEMORY:		return "A memory allocation failed";
    case GLFW_API_UNAVAILABLE:		return "GLFW could not find support for the requested client API on the system";
    case GLFW_VERSION_UNAVAILABLE:	return "The requested OpenGL or OpenGL ES version is not available";
    case GLFW_PLATFORM_ERROR:		return "A platform - specific error occurred that does not match any of the more specific categories";
    case GLFW_FORMAT_UNAVAILABLE:	return "The requested format is not supported or available";
    default:
        return std::string("Unknown error code: ") + std::to_string(code);
    }
}

void GLFWerrorCallback(int error, const char* description)
{
    std::cerr << "GLFW error (" << GLFWerrorCodeToString(error) << "): " << description << std::endl;
}

void framebufferSizeCallback(GLFWwindow* handle, int width, int height)
{
    (void) width;
    (void) height;

    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);
    (void) window;
}

void windowCloseCallback(GLFWwindow* handle)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);
    (void) window;
}

void windowFocusCallback(GLFWwindow* handle, int state)
{
    (void) state;

    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);
    (void) window;
}

void windowPositionCallback(GLFWwindow* handle, int x, int y)
{
    (void) x;
    (void) y;

    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);
    (void) window;
}

void windowRefreshCallback(GLFWwindow* handle)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);
    (void) window;
}

void windowResizeCallback(GLFWwindow* handle, int width, int height)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);

    window->onResize(width, height);
}

}



bool Window::sGLFWisInitialized = false;

Window::Window(const std::string& title, int width, int height, const ContextAttributes& context, const StyleAttributes& style, const FrameBufferAttributes& framebuffer)
    : mTitle(title), mFullScreen(false)
{
    if(!sGLFWisInitialized)
    {
        std::runtime_error("Failed to create Window. GLFW is not initialized!");
    }

    style.apply();
    context.apply();
    framebuffer.apply();

    mHandle = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    assert(mHandle);

    glfwSetWindowUserPointer(mHandle, this); // associates this Window class with the handle

    glfwSetWindowSizeCallback(mHandle, &detail::windowResizeCallback);
    glfwSetFramebufferSizeCallback(mHandle, &detail::framebufferSizeCallback);
    glfwSetWindowCloseCallback(mHandle, &detail::windowCloseCallback);
    glfwSetWindowFocusCallback(mHandle, &detail::windowFocusCallback);
    glfwSetWindowPosCallback(mHandle, &detail::windowPositionCallback);
    glfwSetWindowRefreshCallback(mHandle, &detail::windowRefreshCallback);

    mKeyboard.reset(new Keyboard(*this));
    mMouse.reset(new Mouse(*this));

    makeCurrent();
}

Window::Window(Window&& w) :
    mHandle(std::move(w.mHandle)), mTitle(std::move(w.mTitle))
{
    w.mHandle = nullptr;

    if (mHandle)
    {
        glfwSetWindowUserPointer(mHandle, this);
    }

    mKeyboard.reset(new Keyboard(*this));
    mMouse.reset(new Mouse(*this));

    makeCurrent();
}

Window& Window::operator = (Window&& w)
{
    if (mHandle)
    {
        glfwDestroyWindow(mHandle);
    }

    mHandle = std::move(w.mHandle);
    mTitle = std::move(w.mTitle);

    if (mHandle)
    {
        glfwSetWindowUserPointer(mHandle, this);
    }

    mKeyboard.reset(new Keyboard(*this));
    mMouse.reset(new Mouse(*this));

    makeCurrent();

    return *this;
}

Window::~Window()
{
    if (mHandle)
    {
        glfwDestroyWindow(mHandle);
    }
}

GLFWwindow* Window::handle() const
{
    return mHandle;
}

bool Window::closed() const
{
    return (glfwWindowShouldClose(mHandle) != 0);
}

void Window::close(bool value)
{
    glfwSetWindowShouldClose(mHandle, value ? 1 : 0);
}

void Window::title(const std::string& title)
{
    mTitle = title;
    glfwSetWindowTitle(mHandle, title.c_str());
}

void Window::position(int x, int y)
{
    glfwSetWindowPos(mHandle, x, y);
}

void Window::size(int width, int height)
{
    glfwSetWindowSize(mHandle, width, height);
}

void Window::iconify()
{
    glfwIconifyWindow(mHandle);
}

void Window::restore()
{
    glfwRestoreWindow(mHandle);
}

void Window::show()
{
    glfwShowWindow(mHandle);
}

void Window::hide()
{
    glfwHideWindow(mHandle);
}

namespace
{
    static GLFWwindow* gLastActiveContext = nullptr;
}

void Window::makeCurrent()
{
    if (gLastActiveContext != mHandle)
    {
        glfwMakeContextCurrent(mHandle);
        gLastActiveContext = mHandle;
    }
}

void Window::swapBuffers()
{
    glfwSwapBuffers(mHandle);
}

glm::ivec2 Window::position() const
{
    int x = 0;
    int y = 0;

    glfwGetWindowPos(mHandle, &x, &y);

    return glm::ivec2(x, y);
}

glm::ivec2 Window::size() const
{
    int width = 0;
    int height = 0;

    glfwGetWindowSize(mHandle, &width, &height);

    return glm::ivec2(width, height);
}

glm::ivec2 Window::framebufferSize() const
{
    int width = 0;
    int height = 0;

    glfwGetFramebufferSize(mHandle, &width, &height);

    return glm::ivec2(width, height);
}

float Window::aspectRatio() const
{
    auto winSize = size();
    return static_cast<float>(winSize.x)/static_cast<float>(winSize.y);
}

Keyboard& Window::keyboard() const
{
    return *mKeyboard;
}

Mouse& Window::mouse() const
{
    return *mMouse;
}

bool Window::focused() const
{
    return (glfwGetWindowAttrib(mHandle, GLFW_FOCUSED) != 0);
}

bool Window::iconified() const
{
    return (glfwGetWindowAttrib(mHandle, GLFW_ICONIFIED) != 0);
}

bool Window::visible() const
{
    return (glfwGetWindowAttrib(mHandle, GLFW_VISIBLE) != 0);
}

void Window::vsync(bool enable)
{
    makeCurrent();
    glfwSwapInterval(enable ? 1 : 0);
}

void Window::fullscreen(bool enable)
{
    mFullScreen = enable;

    if(enable)
    {
        mBackupSize = size();
        mBackupPos = position();

        auto monitor = glfwGetPrimaryMonitor();
        const GLFWvidmode * mode = glfwGetVideoMode(monitor);
        glfwSetWindowMonitor( mHandle, monitor, 0, 0, mode->width, mode->height, mode->refreshRate );
    }
    else
    {
        glfwSetWindowMonitor( mHandle, nullptr, mBackupPos.x, mBackupPos.y, mBackupSize.x, mBackupSize.y, 0 );
    }
}

bool Window::fullscreen() const
{
    return mFullScreen;
}


void Window::onKey(int, int, bool)
{

}

void Window::onChar(unsigned int)
{

}

void Window::onMouseButton(int, int, int, int, bool)
{

}

void Window::onMouseScroll(double)
{

}

void Window::onResize(int, int)
{

}

void Window::initializedGLFW()
{
    glfwSetErrorCallback(&detail::GLFWerrorCallback);

    if(!glfwInit())
    {
        throw std::runtime_error("failed to initialize GLFW");
    }

    sGLFWisInitialized = true;
}

void Window::shutdownGLFW()
{
    glfwTerminate();
    sGLFWisInitialized = false;
}

void Window::pollEventsGLFW()
{
    glfwPollEvents();
}

Window::ContextAttributes::ContextAttributes()
    : mAPI(eClientAPI::OPENGL),
      mVersionMajor(3), mVersionMinor(0),
      mOpenGLForwardCompatible(false), mOpenGLDebugContext(false),
      mProfile(eOpenGLProfile::ANY), mRobustness(eContextRobustness::NO_ROBUSTNESS),
      mReleaseBehaviour(eContextReleaseBehaviour::ANY)
{

}

void Window::ContextAttributes::apply() const
{
    glfwWindowHint(GLFW_CLIENT_API, static_cast<int>(mAPI));
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, mVersionMajor);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, mVersionMinor);

    glfwWindowHint(GLFW_CONTEXT_ROBUSTNESS, static_cast<int>(mRobustness));
    glfwWindowHint(GLFW_CONTEXT_RELEASE_BEHAVIOR, static_cast<int>(mReleaseBehaviour));

    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, mOpenGLForwardCompatible ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, mOpenGLDebugContext ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, static_cast<int>(mProfile));
}

Window::FrameBufferAttributes::FrameBufferAttributes()
    : mRedBits(8), mGreenBits(8), mBlueBits(8), mAlphaBits(8),
      mDepthBits(8), mStencilBits(8),
      mSamples(0), mDoubleBuffer(true), mRefreshRate(60)
{

}

void Window::FrameBufferAttributes::apply() const
{
    glfwWindowHint(GLFW_RED_BITS,       mRedBits);
    glfwWindowHint(GLFW_GREEN_BITS,     mGreenBits);
    glfwWindowHint(GLFW_BLUE_BITS,      mBlueBits);
    glfwWindowHint(GLFW_ALPHA_BITS,     mAlphaBits);
    glfwWindowHint(GLFW_DEPTH_BITS,     mDepthBits);
    glfwWindowHint(GLFW_STENCIL_BITS,   mStencilBits);

    glfwWindowHint(GLFW_SAMPLES,        mSamples);
    glfwWindowHint(GLFW_DOUBLEBUFFER,   mDoubleBuffer ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_REFRESH_RATE,   mRefreshRate);
}

Window::StyleAttributes::StyleAttributes()
    : mResizable(true), mVisible(true), mDecorated(true),
      mFocused(true), mAutoIconify(true), mFloating(true)
{

}

void Window::StyleAttributes::apply() const
{
    glfwWindowHint(GLFW_RESIZABLE,      mResizable ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_VISIBLE,        mVisible ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_DECORATED,      mDecorated ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_FOCUSED,        mFocused ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_AUTO_ICONIFY,   mAutoIconify ? GL_TRUE : GL_FALSE);
    glfwWindowHint(GLFW_FLOATING,       mFloating ? GL_TRUE : GL_FALSE);
}
