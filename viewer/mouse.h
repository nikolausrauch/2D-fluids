/**
 * @author: Nikolaus Rauch
 * @date: 12.02.2021
 */
#pragma once

#include "base.h"

#include <memory>
#include <array>

class Window;


/**
 * @brief Mouse class to query state (pressed buttons, and wheel scroll) of mouse
 */
class Mouse
{
public:
    enum eButton : int
    {
        LEFT = GLFW_MOUSE_BUTTON_1,
        RIGHT,
        MIDDLE,
        OTHER_1,
        OTHER_2,
        OTHER_3,
        OTHER_4,
        OTHER_5
    };

    enum class eCursorState : int
    {
        VISIBLE     = GLFW_CURSOR_NORMAL,
        HIDDEN      = GLFW_CURSOR_HIDDEN,
        DISABLED    = GLFW_CURSOR_DISABLED
    };

    Mouse(const Mouse&) = delete;
    Mouse& operator = (const Mouse&) = delete;

    /**
     * @brief set mouse position with respect to top-left window corner
     * @param x coordinate
     * @param y coordinate
     */
    void position(double x, double y);


    /**
     * @brief query mouse position with respect to top-left window corner
     * @param x coordinate
     * @param y coordinate
     */
    const glm::dvec2& position();

    /**
     * @brief set cursor state (visible, hidden, disabled)
     * @param state as enum eCursorState
     */
    void cursorState(eCursorState state);

    /**
     * @brief query cursor state (visible, hidden, disabled)
     * @return state as enum eCursorState
     */
    eCursorState cursorState() const;

    /**
     * @brief set croll value of mouse (should not be required to use)
     * @param float value offset during one frame (OS dependent value)
     */
    void scroll(float scroll);

    /**
     * @brief query croll value of mouse (should not be required to use)
     * @return float value offset during one frame (OS dependent value)
     */
    float scroll() const;

    /**
     * @brief operator [] to query state of button
     * @param button is a GLFW button code or enum eButton (e.g. GLFW_MOUSE_BUTTON_1 or Mouse::eButton::LEFT)
     * @return true if pressed, false otherwise
     */
    bool operator[] (std::size_t button) const;
    bool& operator[] (std::size_t button);

    /**
     * @brief get context window responsable for this keyboard
     * @return window reference
     */
    Window& window() const;

private:
    glm::dvec2 mPosition;
    std::array<bool, GLFW_MOUSE_BUTTON_LAST>  mButtonState;
    float mScroll;

    Window& mSourceWindow;
    std::unique_ptr<GLFWcursor, void(*)(GLFWcursor*)> mMouseCursor;


    Mouse(Window &sourceWindow);
    friend class Window;
};
