/**
 * @author: Nikolaus Rauch
 * @date: 12.02.2021
 *
 * Simple Window (GLFW) Wrapper Class
 */
#pragma once

#include "base.h"
#include "mouse.h"
#include "keyboard.h"

#include <memory>
#include <functional>
#include <string>


enum class eClientAPI : int
{
    OPENGL          = GLFW_OPENGL_API,
    OPENGL_ES       = GLFW_OPENGL_ES_API
};

enum class eOpenGLProfile : int
{
    ANY             = GLFW_OPENGL_ANY_PROFILE,
    COMPATIBILITY   = GLFW_OPENGL_COMPAT_PROFILE,
    CORE            = GLFW_OPENGL_CORE_PROFILE
};

enum class eContextRobustness : int
{
    NO_ROBUSTNESS           = GLFW_NO_ROBUSTNESS,
    NO_RESET_NOTIFICATION   = GLFW_NO_RESET_NOTIFICATION,
    LOSE_CONTEXT_ON_RESET   = GLFW_LOSE_CONTEXT_ON_RESET
};

enum class eContextReleaseBehaviour : int
{
    ANY     = GLFW_ANY_RELEASE_BEHAVIOR,
    FLUSH   = GLFW_RELEASE_BEHAVIOR_FLUSH,
    NONE    = GLFW_RELEASE_BEHAVIOR_NONE
};

/**
 * @brief Window class that encapsulates GLFW window and its events
 *
 * creates a OpenGL window with the fine grained control of context attributes, framebuffer attributes, and window style settings
 */
class Window
{
public:
    struct ContextAttributes
    {
        ContextAttributes();
        void apply() const;

        eClientAPI mAPI;
        int mVersionMajor;
        int mVersionMinor;

        bool mOpenGLForwardCompatible;
        bool mOpenGLDebugContext;

        eOpenGLProfile mProfile;
        eContextRobustness mRobustness;
        eContextReleaseBehaviour mReleaseBehaviour;
    };

    struct FrameBufferAttributes
    {
        FrameBufferAttributes();
        void apply() const;

        int mRedBits;
        int mGreenBits;
        int mBlueBits;
        int mAlphaBits;
        int mDepthBits;
        int mStencilBits;
        int mSamples;

        bool mDoubleBuffer;
        int mRefreshRate;
    };

    struct StyleAttributes
    {
        StyleAttributes();
        void apply() const;

        bool mResizable;
        bool mVisible;
        bool mDecorated;
        bool mFocused;
        bool mAutoIconify;
        bool mFloating;
    };

    /**
     * @brief create window with opengl context
     *
     * @param title window title
     * @param width width in pixel
     * @param height height in pixel
     * @param context OpenGL context attributes
     * @param style window style settings
     * @param framebuffer framebuffer attributes
     */
    Window(const std::string& title, int width, int height,
           const ContextAttributes& context = ContextAttributes(),
           const StyleAttributes& style = StyleAttributes(),
           const FrameBufferAttributes& framebuffer = FrameBufferAttributes());

    Window(const Window&) = delete;
    Window(Window&& w);
    Window& operator = (const Window&) = delete;
    Window& operator = (Window&& w);
    virtual ~Window();

    /**
     * @brief retrieve low lewel GLFW window handle
     * @return GLFWwindow pointer
     */
    GLFWwindow* handle() const;

    /**
     * @brief get keyboard state container
     * @return Keyboard reference
     */
    Keyboard& keyboard() const;

    /**
     * @brief get mouse state container
     * @return Mouse reference
     */
    Mouse& mouse() const;

    /**
     * @brief close window
     * @param if true window is closed in next frame
     */
    void close(bool value);

    /**
     * @brief query closed state
     * @return boolean (true window will be destroyed)
     */
    bool closed() const;

    /**
     * @brief set title of window
     * @param title title string
     */
    void title(const std::string& title);

    /**
     * @brief set window position
     * @param x coordinate in display pixels (top-left 0,0)
     * @param y coordinate in display pixels (top-left 0,0)
     */
    void position(int x, int y);

    /**
     * @brief modify window size in pixel
     * @param width width in pixels
     * @param height height in pixels
     */
    void size(int width, int height);

    /* window compositor commands */
    void iconify();
    void restore();
    void show();
    void hide();

    /**
     * @brief make this windows context current on calling thread
     */
    void makeCurrent();

    /**
     * @brief swap framebuffer contents (flushes submitted draw calls)
     */
    void swapBuffers();

    /**
     * @brief retrieve window position in display pixel coordinates
     * @return glm::ivec2 with pixel coordinates
     */
    glm::ivec2 position() const;

    /**
     * @brief retrieve window size in display pixels
     * @return glm::ivec2 with size in pixel
     */
    glm::ivec2 size() const;

    /**
     * TODO: window != framebuffer size on retina displays
     *
     * @brief retrieve frambuffer size in pixels
     * @return glm::ivec2 with size in pixel
     */
    glm::ivec2 framebufferSize() const;

    /**
     * @brief retrieve current aspect ratio of window
     * @return float value width/height
     */
    float aspectRatio() const;

    /**
     * @brief is window focused
     * @return true if focused in compositon
     */
    bool focused() const;

    /**
     * @brief window is iconified
     * @return true if iconified
     */
    bool iconified() const;

    /**
     * @brief is window visible in the compositor
     * @return true if visible
     */
    bool visible() const;

    /**
     * @brief enable vsync; depending on your driver settings (Nvidia, AMD, Intel) this may not have an effect
     * @param enable (true to activate vsync)
     */
    void vsync(bool enable);

    /**
     * @brief windowed or fullscreen mode of window
     * @param enable (true to enable fullscreen mode)
     */
    void fullscreen(bool enable);

    /**
     * @brief retrieve if fullscreen mode is activate
     * @return true if fullscreen mode active
     */
    bool fullscreen() const;

    /**
     * @brief called on key event (override to handle)
     *
     * @param key GLFW key code
     * @param mod pressed modifiers
     * @param press true if pressed
     */
    virtual void onKey(int key, int mod, bool press);

    /**
     * @brief called on character enter event (override to handle)
     *
     * @param c unsigned int char code
     */
    virtual void onChar(unsigned int c);

    /**
     * @brief called on mouse button event (override to handle)
     * @param x position in window pixel coordinates
     * @param y position in window pixel coordinates
     * @param button GLFW button code
     * @param mod pressed modifiers
     * @param press true if pressed
     */
    virtual void onMouseButton(int x, int y, int button, int mod, bool press);

    /**
     * @brief called on mouse wheel scroll (override to handle)
     * @param delta offset from previous frame
     */
    virtual void onMouseScroll(double delta);

    /**
     * @brief called on window resize event (override to handle)
     * @param width window width in pixels
     * @param height window height in pixels
     */
    virtual void onResize(int width, int height);

    /**
     * @brief initialize GLFW library (only call once per application, not per window)
     */
    static void initializedGLFW();

    /**
     * @brief shutdown GLFW library (call after all windows are closed)
     */
    static void shutdownGLFW();

    /**
     * @brief query OS events, which triggers corresponding callback functions
     */
    static void pollEventsGLFW();

private:
    GLFWwindow* mHandle;

    std::string mTitle;

    std::unique_ptr<Keyboard> mKeyboard;
    std::unique_ptr<Mouse> mMouse;

    bool mFullScreen;
    glm::vec2 mBackupSize;
    glm::vec2 mBackupPos;

    static bool sGLFWisInitialized;
};

