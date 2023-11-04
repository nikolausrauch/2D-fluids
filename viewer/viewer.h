/**
 * @author: Nikolaus Rauch (main), Marcel Ritter
 * @date: 15.10.2021
 *
 * Simple Viewer (OpenGL2 Context, ImGui, ImPlot)
 */
#pragma once

#include "window.h"

#include <string>
#include <functional>
#include <vector>

#include <imgui/imgui.h>
#include <imgui/implot.h>

namespace detail { struct WindowImGui; }


/**
 * @brief Viewer class encapsulating window, input, primitive rendering and ui (imgui) handling
 *
 * Simple Legacy OpenGL viewer that provides event handling through callback functions.
 * The events are:
 *      - init: is called after the opengl context is initialized
 *      - shutdown: is called right before the opengl context is destroyed
 *      - key: triggered if key is pressed or released
 *      - mouseButton: triggered if mouse button is pressed
 *      - update: called before rendering
 *      - draw: use to submit render calls (points, lines, quads, triangles, circles, ...)
 *      - gui: used to draw/respond to user-interface via ImGui and ImPlot
 *
 * How to use:
 *
 *     Viewer viewer;
 *     viewer.mWindow.title = "Example";
 *     viewer.mWindow.width = 1280;
 *     viewer.mWindow.height = 720;
 *     viewer.mWindow.vsync = true;
 *     viewer.mWindow.mHDPI = false;
 *
 *     // install callback functions for events
 *     viewer.onUpdate([](Window& window, double dt){ ... });
 *     viewer.onDraw([](Window& window, double dt){  ... });
 *     viewer.onGui([](Window& window, double dt){ ... });
 *     viewer.onMouseButton([](Window& window, Mouse& mouse, int button, int mod, bool press){ ... });
 *     viewer.onKey([&](Window& window, Keyboard& keyboard, int key, int mod, bool press){ ... });
 *
 *     // start viewer (OpenGL context is only created now, use drawcommands and loading functions only in callbacks)
 *     viewer.run();
 *
 */
class Viewer
{
    typedef std::function<void (Window& window, Keyboard& keyboard, int key, int mod, bool press)> KeyboardCallback;
    typedef std::function<void (Window& window, Mouse& mouse, int button, int mod, bool press)> MouseCallback;
    typedef std::function<void (Window& window, double dt)> UpdateCallback;
    typedef std::function<void (Window& window, double dt)> DrawCallback;
    typedef std::function<void (Window& window, double dt)> GuiCallback;
    typedef std::function<void (void)> SetupCallback;

public:
    Viewer();
    ~Viewer();

    /**
     * @brief initialize GL context and starts main loop
     *        repeatedly calling installed callback functions (e.g. onKey, onDraw, etc.)
     *        runs until window is closed
     */
    void run();

    /**
     * @brief retrieve current frame per second estimate
     * @return float representing current number of frames rendererd during one second
     */
    double fps() const;

    /**
     * @brief retrieve elapsed time from previous frame
     * @return elapsed time in seconds
     */
    double dt() const;

    /**
     * @brief get the world space position from screen space coordinates
     * @param position in pixel with respect to top-left window corner
     * @return reprojected position
     */
    glm::vec2 worldSpacePosition(const glm::dvec2& windowPos);


    /**
     * @brief install callback function for GL context required initialization
     * @param function reference: void (void)
     */
    void onInit(const SetupCallback& initCB);

    /**
     * @brief install callback function to cleanup GL related resources
     * @param function reference: void (void)
     */
    void onShutdown(const SetupCallback& shutdownCB);

    /**
     * @brief install callback function triggered on key event
     *
     * @param function reference: void (Window& window, Keyboard& keyboard, int key, int mod, bool press)
     *        - window: context responsable for event
     *        - keyboard: state container of keyboard
     *        - key: key code related to event (e.g. GLFW_KEY_0 or GLFW_KEY_W)
     *        - mod: modifier (ctrl, alt, etc.)
     *        - press: true if pressed, false otherwise
     */
    void onKey(const KeyboardCallback& keyCB);

    /**
     * @brief install callback function triggered on mouse button event
     *
     * @param function reference: void (Window& window, Mouse& mouse, int button, int mod, bool press):
     *        - window: context responsable for event
     *        - mouse: state container of mouse
     *        - button: button code related to event (e.g. GLFW_MOUSE_BUTTON_1 or Mouse::eButton::LEFT)
     *        - mod: modifier (ctrl, alt, etc.)
     *        - press: true if pressed, false otherwise
     */
    void onMouseButton(const MouseCallback& mouseCB);

    /**
     * @brief callback initiated once per frame (for animation, physic steps, etc.)
     *
     * @param function reference: void (Window& window, double dt)
     *        - window: context responsable for event
     *        - dt: elapsed time from previous frame
     */
    void onUpdate(const UpdateCallback& updateCB);

    /**
     * @brief install callback to draw primitives with (drawPoints, drawLines, etc.)
     * @param function reference: void (Window& window, double dt)
     *        - window: context responsable for event
     *        - dt: elapsed time from previous frame
     */
    void onDraw(const DrawCallback& drawCB);

    /**
     * @brief install callback to draw/respond to user-interface (ImGui/ImPlot)
     *        use ImGui and ImPlot commands to draw your user interface (see ImGui and ImPlot examples)
     *
     * @param function reference: void (Window& window, double dt)
     *        - window: context responsable for event
     *        - dt: elapsed time from previous frame
     */
    void onGui(const GuiCallback& guiCB);

    /**
     * @brief load texture resource from filepath (only call after GL context has been created)
     * @param filepath to supported image file (.png, .jpeg, etc.)
     * @return resource identifier
     */
    unsigned int loadTexture(const std::string& path);

    /**
     * @brief load texture resource from provided data
     * @param width of texture
     * @param height of texture
     * @param data rgba unsigned char values
     * @return resource identifier
     */
    unsigned int createTexture(unsigned int width, unsigned int height, const unsigned char* data = nullptr);

    /**
     * @brief resize texture resource with provided data
     * @param width of texture
     * @param height of texture
     * @param data rgba unsigned char values
     * @return resource identifier
     */
    void resizeTexture(unsigned int id, unsigned int width, unsigned int height, const unsigned char* data = nullptr);

    /**
     * @brief update texture resource
     * @param id of texture resource
     * @param width of texture
     * @param height of texture
     * @param data rgba unsigned char values
     * @return resource identifier
     */
    void updateTexture(unsigned int id, unsigned int width, unsigned int height, const unsigned char* data);

    /**
     * @brief draw points from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and position and color reference which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& coord, glm::vec4& color)
     *         - data_type: type in container
     *         - coord: position in world space of point
     *         - color: color of point
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawPoints(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw lines from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and color references which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& start, glm::vec2& end, glm::vec4& color)
     *         - data_type: type in container
     *         - start: start position in world space of line
     *         - end: end position in world space of line
     *         - color: color of line
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawLines(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw line
     *
     * @param p0: first vertex position in worldspace
     * @param p1: second vertex position in worldspace
     * @param color: color of line
     */
    void drawLine(const glm::vec2& p0, const glm::vec2& p1, const glm::vec4& color);

    /**
     * @brief draw circles from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and color references which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& start, glm::vec2& end, glm::vec4& color)
     *         - data_type: type in container
     *         - position: position in world space of circle
     *         - radius: radius of circle
     *         - color: color of circle
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawCircles(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw triangles from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and color references which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec4& color)
     *         - data_type: type in container
     *         - p0: first vertex position in worldspace
     *         - p1: second vertex position in worldspace
     *         - p2: thrid vertex position in worldspace
     *         - color: color of triangle
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawTriangles(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw textured triangles from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and color references which need to be set
     *
     * @param textureID texture identifier from loadTexture
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& uv0, glm::vec2& uv1, glm::vec2& uv2)
     *         - data_type: type in container
     *         - p0: first vertex position in worldspace
     *         - p1: second vertex position in worldspace
     *         - p2: thrid vertex position in worldspace
     *         - uv0: texture coordinate for first vertex
     *         - uv1: texture coordinate for second vertex
     *         - uv2: texture coordinate for third vertex
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawTrianglesTextured(unsigned int textureID, InputIt first, InputIt last, UnaryFunction f);


    /**
     * @brief draw quads from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and texture coordinates references which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, glm::vec4& color)
     *         - data_type: type in container
     *         - p0: first vertex position in worldspace
     *         - p1: second vertex position in worldspace
     *         - p2: thrid vertex position in worldspace
     *         - uv0: texture coordinate for first vertex
     *         - uv1: texture coordinate for second vertex
     *         - uv2: texture coordinate for third vertex
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawQuads2D(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw textured quads from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and positions and texture coordinates references which need to be set
     *
     * @param textureID texture identifier from loadTexture
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, unsigned int textureID)
     *         - data_type: type in container
     *         - p0: first vertex position in worldspace
     *         - p1: second vertex position in worldspace
     *         - p2: thrid vertex position in worldspace
     *         - textureID: texture resource in use
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawQuads2DTextured(unsigned int textureID, InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw connected lines from generic container content
     *
     * @tparam InputIt Iterator of generic container (needs to support operators *, !=; ++)
     * @tparam UnaryFunction function taking as first argument a const reference to the data type in the container
     *         and vector of vertices and color references which need to be set
     *
     * @param first iterator start
     * @param last iterator stop
     * @param f function reference: bool (const auto& data_type, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, unsigned int texture ud color)
     *         - data_type: type in container
     *         - p0: first vertex position in worldspace
     *         - p1: second vertex position in worldspace
     *         - p2: thrid vertex position in worldspace
     *         - uv0: texture coordinate for first vertex
     *         - uv1: texture coordinate for second vertex
     *         - uv2: texture coordinate for third vertex
     *         - return false if draw command should be skipped
     */
    template<class InputIt, class UnaryFunction>
    void drawOutline(InputIt first, InputIt last, UnaryFunction f);

    /**
     * @brief draw a single point at a world space position with a specified color
     * @param pos position in worldspace of point
     * @param color color of point
     */
    void drawPoint(const glm::vec2& pos, const glm::vec4& color);

    /**
     * @brief draw a single circle at a given world space position with a specified color
     *
     * @param position in world space coordinates
     * @param radius of circle
     * @param color of circle
     */
    void drawCircle(const glm::vec2 &position, float radius, glm::vec4& color);

    /**
     * @brief draw a single circle (only outline) at a given world space position with a specified color
     *
     * @param position in world space coordinates
     * @param radius of circle
     * @param color of circle
     */
    void drawCircleOutline(const glm::vec2 &position, float radius, glm::vec4& color);

    /**
     * @brief draw a single quad at a world space position and relative vertex coordinates with a specified color
     *
     * @param p0 first vertex position in world space
     * @param p1 second vertex position in world space
     * @param p2 third vertex position in world space
     * @param p3 fourth vertex position in world space
     * @param color quad color
     */
    void drawQuad(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3, const glm::vec4& color);

    /**
     * @brief draw a single texture quad at a world space position and relative vertex coordinates with a specified texture
     *
     * @param textureID texture resource to use
     * @param p0 first vertex position in world space
     * @param p1 second vertex position in world space
     * @param p2 third vertex position in world space
     * @param p3 fourth vertex position in world space
     */
    void drawQuad2DTextured(unsigned int textureID, const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3);

    /**
     * @brief draw visual boundary to indicate collission planes
     *
     * @param pos position in world space
     * @param normal direction of boundary
     * @param size visual length of boundary
     * @param pseudoShadingSize length of shading lines
     */
    void drawBoundary(const glm::vec2& pos, const glm::vec2& normal, float size, float pseudoShadingSize = -1.0f);


private:
    float renderScale() const;



public:
    /* window settings */
    struct
    {
        std::string title;
        int width;
        int height;
        bool vsync;
        bool mHDPI;
    } mWindow;

    /* render settings */
    struct
    {
        float pointRadius;
        float lineWidth;
        unsigned int circleVertices;

        bool wireframe;

        glm::vec3 bgColor;
    } mRender;

    /* camera controls */
    struct
    {
        glm::vec2 position;
        glm::vec2 size;
        float zoom;
    } mCamera;

    struct
    {
        bool keyboardCaptured;
        bool mouseCaptured;
    } mUI;

private:
    struct
    {
        double accumTime;
        double fps;
        unsigned int frames;

        double dt;
    } mFPSCounter;

    std::vector<GLuint> mTextureStorage;

private:
    MouseCallback mOnMouseButton;
    KeyboardCallback mOnKey;
    UpdateCallback mOnUpdate;
    DrawCallback mOnDraw;
    GuiCallback mOnGui;
    SetupCallback mOnInit;
    SetupCallback mOnShutdown;

    friend struct detail::WindowImGui;
};

#include "viewer.inl"

