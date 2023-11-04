#include "keyboard.h"
#include "window.h"

#include <cassert>

namespace detail
{

void keyCallback(GLFWwindow* handle, int key, int scannCode, int action, int mod)
{
    (void) scannCode;

    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);

    Keyboard& kb = window->keyboard();
    kb[static_cast<std::size_t>(key)] = (action != GLFW_RELEASE);

    window->onKey(key, mod, (action != GLFW_RELEASE));
}

void charCallback(GLFWwindow* handle, unsigned int character)
{
    Window* window = static_cast<Window*>(glfwGetWindowUserPointer(handle));
    assert(window);

    window->onChar(character);
}

}

Keyboard::Keyboard(Window &sourceWindow) : mSourceWindow(sourceWindow)
{
    std::uninitialized_fill(std::begin(mKeyState), std::end(mKeyState), false);

    glfwSetKeyCallback(sourceWindow.handle(), detail::keyCallback);
    glfwSetCharCallback(sourceWindow.handle(), detail::charCallback);
}

void Keyboard::stickyKeys(bool enabled)
{
    glfwSetInputMode(mSourceWindow.handle(), GLFW_STICKY_KEYS, enabled ? GL_TRUE : GL_FALSE);
}

bool Keyboard::stickyKeys() const
{
    return glfwGetInputMode(mSourceWindow.handle(), GLFW_STICKY_KEYS);
}

Window& Keyboard::window() const
{
    return mSourceWindow;
}

bool Keyboard::operator[](std::size_t key) const
{
    assert(0 <= key && key < GLFW_KEY_LAST);
    return mKeyState[key];
}

bool& Keyboard::operator[](std::size_t button)
{
    assert(0 <= button && button < GLFW_KEY_LAST);
    return mKeyState[button];
}
