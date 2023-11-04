/**
 * @author: Nikolaus Rauch
 * @date: 12.02.2021
 */
#pragma once

#include "base.h"

#include <array>

class Window;

/**
 * @brief Keyboard class to query state (pressed keys) of keyboard
 */
class Keyboard
{
public:
    Keyboard(const Keyboard&) = delete;
    Keyboard& operator = (const Keyboard&) = delete;

    void stickyKeys(bool enabled);
    bool stickyKeys() const;

    /**
     * @brief get context window responsable for this keyboard
     * @return window reference
     */
    Window& window() const;

    /**
     * @brief operator [] to query state of key
     * @param key is a GLFW key code (e.g. GLFW_KEY_0 or GLFW_KEY_W)
     * @return true if pressed, false otherwise
     */
    bool operator[](std::size_t key) const;
    bool& operator[](std::size_t key);

private:
    Window& mSourceWindow;
    std::array<bool, GLFW_KEY_LAST> mKeyState;

    Keyboard(Window& sourceWindow);
    friend class Window;
};

