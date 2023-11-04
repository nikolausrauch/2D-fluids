#include "mouse.h"
#include "window.h"

#include <cassert>
#include <iostream>

namespace detail
{

void mouseButtonCallback(GLFWwindow* handle, int button, int action, int mods)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);

    Mouse& mouse = window->mouse();
    auto pos = mouse.position();
    mouse[static_cast<std::size_t>(button)] = (action != GLFW_RELEASE);

    window->onMouseButton(pos.x, pos.y, button, mods, (action != GLFW_RELEASE));
}

void mouseScrollCallback(GLFWwindow* handle, double xoffset, double yoffset)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);

    Mouse& mouse = window->mouse();
    mouse.scroll(xoffset);

    window->onMouseScroll(yoffset);
}

}


Mouse::Mouse(Window& sourceWindow) :
    mPosition(0, 0),
    mButtonState{false, false, false, false, false, false, false}, mScroll(0.0),
    mSourceWindow(sourceWindow), mMouseCursor(glfwCreateStandardCursor(GLFW_ARROW_CURSOR), glfwDestroyCursor)
{
    glfwSetMouseButtonCallback(sourceWindow.handle(), detail::mouseButtonCallback);
    glfwSetScrollCallback(sourceWindow.handle(), detail::mouseScrollCallback);

    std::uninitialized_fill(std::begin(mButtonState), std::end(mButtonState), false);
    mPosition = position();
}

void Mouse::position(double x, double y)
{
    glfwSetCursorPos(mSourceWindow.handle(), x, y);
    mPosition.x = x;
    mPosition.y = y;
}

const glm::dvec2& Mouse::position()
{
    glfwGetCursorPos(mSourceWindow.handle(), &mPosition.x, &mPosition.y);
    return mPosition;
}

void Mouse::cursorState(Mouse::eCursorState state)
{
    glfwSetInputMode(mSourceWindow.handle(), GLFW_CURSOR, static_cast<int>(state));
}

Mouse::eCursorState Mouse::cursorState() const
{
    // TODO: problematic if glfw ever support more states
    return static_cast<eCursorState>(glfwGetInputMode(mSourceWindow.handle(), GLFW_CURSOR));
}

void Mouse::scroll(float scroll)
{
    mScroll = scroll;
}

float Mouse::scroll() const
{
    return mScroll;
}

bool Mouse::operator[](std::size_t button) const
{
    assert(0 <= button && button < GLFW_MOUSE_BUTTON_LAST);
    return mButtonState[button];
}

bool &Mouse::operator[](std::size_t button)
{
    assert(0 <= button && button < GLFW_MOUSE_BUTTON_LAST);
    return mButtonState[button];
}

Window& Mouse::window() const
{
    return mSourceWindow;
}
