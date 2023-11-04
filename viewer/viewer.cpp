#include "viewer.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <stb_image/stb_image.h>

namespace detail
{

/* Window Wrapped in ImGui functionality */
struct WindowImGui : public Window
{
    WindowImGui(Viewer& viewer,
                const std::string& title, int width, int height, bool hdpi,
                const ContextAttributes& context = ContextAttributes(),
                const StyleAttributes& style = StyleAttributes(),
                const FrameBufferAttributes& framebuffer = FrameBufferAttributes())
        : Window(title, width, height, context, style, framebuffer), mViewer(viewer)
    {
        initImGui(*this, hdpi);
    }

    ~WindowImGui()
    {
        shutdownImGui();
    }

    void onKey(int key, int mod, bool press) override
    {
        ImGuiIO& io = ImGui::GetIO();
        io.KeysDown[key] = press;

        // TODO: imgui modifiers missing
        if(mViewer.mOnKey) mViewer.mOnKey(*this, keyboard(), key, mod, press);
    }

    void onChar(unsigned int c) override
    {
        ImGuiIO& io = ImGui::GetIO();
        io.AddInputCharacter(c);
    }

    void onMouseButton(int x, int y, int button, int mod, bool pressed) override
    {
        (void) x; (void) y;

        std::size_t idx = static_cast<std::size_t>(button);
        mImGuiData.mMouseJustPressed[idx] = pressed || mImGuiData.mMouseJustPressed[idx];

        if(mViewer.mOnMouseButton) mViewer.mOnMouseButton(*this, mouse(), button, mod, pressed);
    }

    void onMouseScroll(double delta) override
    {
        ImGuiIO& io = ImGui::GetIO();
        io.MouseWheel += static_cast<float>(delta);
    }

    void onResize(int width, int height) override
    {
        glViewport(0, 0, width, height);

        /* rescale camera view */
        glm::vec2 scaleChange{ width / static_cast<float>(mViewer.mWindow.width), height / static_cast<float>(mViewer.mWindow.height)};
        mViewer.mCamera.size *= scaleChange;

        mViewer.mWindow.width = width;
        mViewer.mWindow.height = height;
    }

    void initImGui(Window &window, float hdpi = false)
    {
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGui::StyleColorsClassic();

        /**** GLFW ImGui Bindings ****/
        ImGuiIO& io = ImGui::GetIO();
        io.BackendFlags |= ImGuiBackendFlags_HasMouseCursors;
        io.BackendPlatformName = "imgui_impl_glfw";
        io.ClipboardUserData = window.handle();

        /* Keyboard mapping. ImGui will use those indices to peek into the io.KeysDown[] array.*/
        io.KeyMap[ImGuiKey_Tab]         = GLFW_KEY_TAB;
        io.KeyMap[ImGuiKey_LeftArrow]   = GLFW_KEY_LEFT;
        io.KeyMap[ImGuiKey_RightArrow]  = GLFW_KEY_RIGHT;
        io.KeyMap[ImGuiKey_UpArrow]     = GLFW_KEY_UP;
        io.KeyMap[ImGuiKey_DownArrow]   = GLFW_KEY_DOWN;
        io.KeyMap[ImGuiKey_PageUp]      = GLFW_KEY_PAGE_UP;
        io.KeyMap[ImGuiKey_PageDown]    = GLFW_KEY_PAGE_DOWN;
        io.KeyMap[ImGuiKey_Home]        = GLFW_KEY_HOME;
        io.KeyMap[ImGuiKey_End]         = GLFW_KEY_END;
        io.KeyMap[ImGuiKey_Insert]      = GLFW_KEY_INSERT;
        io.KeyMap[ImGuiKey_Delete]      = GLFW_KEY_DELETE;
        io.KeyMap[ImGuiKey_Backspace]   = GLFW_KEY_BACKSPACE;
        io.KeyMap[ImGuiKey_Space]       = GLFW_KEY_SPACE;
        io.KeyMap[ImGuiKey_Enter]       = GLFW_KEY_ENTER;
        io.KeyMap[ImGuiKey_Escape]      = GLFW_KEY_ESCAPE;
        io.KeyMap[ImGuiKey_KeyPadEnter] = GLFW_KEY_KP_ENTER;
        io.KeyMap[ImGuiKey_A]           = GLFW_KEY_A;
        io.KeyMap[ImGuiKey_C]           = GLFW_KEY_C;
        io.KeyMap[ImGuiKey_V]           = GLFW_KEY_V;
        io.KeyMap[ImGuiKey_X]           = GLFW_KEY_X;
        io.KeyMap[ImGuiKey_Y]           = GLFW_KEY_Y;
        io.KeyMap[ImGuiKey_Z]           = GLFW_KEY_Z;


        /*** font loading ***/
        if(hdpi)
        {
            // io.FontGlobalScale *= 2.0f; // blurry font, therfore this workaround
            ImFontConfig config;
            config.SizePixels = 24.0f;
            config.OversampleH = config.OversampleV = 1;
            config.PixelSnapH = true;
            io.Fonts->AddFontDefault(&config);
        }

        /**** OpenGL ImGui Setup ****/
        io.BackendRendererName = "imgui_impl_opengl2";

        unsigned char* pixels;
        int width, height;
        io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);

        /* Upload texture to graphics system */
        glGenTextures(1, &mImGuiData.mFontTexture);
        glBindTexture(GL_TEXTURE_2D, mImGuiData.mFontTexture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
        io.Fonts->TexID = (ImTextureID) (intptr_t) mImGuiData.mFontTexture;
    }

    void shutdownImGui()
    {
        if(mImGuiData.mFontTexture)
        {
            ImGuiIO& io = ImGui::GetIO();
            glDeleteTextures(1, &mImGuiData.mFontTexture);
            io.Fonts->TexID = nullptr;
            mImGuiData.mFontTexture = 0;
        }

        ImPlot::DestroyContext();
        ImGui::DestroyContext();
    };

    void newFrameImGui(Window &window, double dt)
    {
        ImGuiIO& io = ImGui::GetIO();
        assert(io.Fonts->IsBuilt());

        io.DisplaySize = ImVec2( window.size().x, window.size().y );
        io.DisplayFramebufferScale = ImVec2( window.framebufferSize().x / io.DisplaySize.x, window.framebufferSize().y / io.DisplaySize.y );
        io.DeltaTime = static_cast<float>(dt);

        /* update mouse */
        auto& mouse = window.mouse();
        for(std::size_t i = 0; i < IM_ARRAYSIZE(io.MouseDown); i++)
        {
            io.MouseDown[i] = mImGuiData.mMouseJustPressed[i] || mouse[i];
            mImGuiData.mMouseJustPressed[i] = false;
        }

        if(window.focused())
        {
            io.MousePos = ImVec2(static_cast<float>(mouse.position().x), static_cast<float>(mouse.position().y));
        }

        ImGui::NewFrame();
    }

    void renderImGui()
    {
        ImGui::Render();

        ImDrawData* draw_data = ImGui::GetDrawData();

        int fb_width = static_cast<int>(draw_data->DisplaySize.x * draw_data->FramebufferScale.x);
        int fb_height = static_cast<int>(draw_data->DisplaySize.y * draw_data->FramebufferScale.y);
        if (fb_width == 0 || fb_height == 0)
        {
            return;
        }

        /* Backup GL state */
        GLint last_texture; glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
        GLint last_polygon_mode[2]; glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
        GLint last_viewport[4]; glGetIntegerv(GL_VIEWPORT, last_viewport);
        GLint last_scissor_box[4]; glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
        glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_TRANSFORM_BIT);


        /* Setup render state: alpha-blending enabled, no face culling, no depth testing, scissor enabled, vertex/texcoord/color pointers, polygon fill. */
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_COLOR_MATERIAL);
        glEnable(GL_SCISSOR_TEST);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glEnable(GL_TEXTURE_2D);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        /* Our visible imgui space lies from draw_data->DisplayPos (top left) to draw_data->DisplayPos+data_data->DisplaySize (bottom right). DisplayPos is (0,0) for single viewport apps. */
        glViewport(0, 0, static_cast<GLsizei>(fb_width), static_cast<GLsizei>(fb_height));
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(draw_data->DisplayPos.x, draw_data->DisplayPos.x + draw_data->DisplaySize.x, draw_data->DisplayPos.y + draw_data->DisplaySize.y, draw_data->DisplayPos.y, -1.0f, +1.0f);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();


        /* Will project scissor/clipping rectangles into framebuffer space */
        ImVec2 clip_off = draw_data->DisplayPos;         // (0,0) unless using multi-viewports
        ImVec2 clip_scale = draw_data->FramebufferScale; // (1,1) unless using retina display which are often (2,2)


        /* Render command lists */
        for (int n = 0; n < draw_data->CmdListsCount; n++)
        {
            const ImDrawList* cmd_list = draw_data->CmdLists[n];
            const ImDrawVert* vtx_buffer = cmd_list->VtxBuffer.Data;
            const ImDrawIdx* idx_buffer = cmd_list->IdxBuffer.Data;
            glVertexPointer(2, GL_FLOAT, sizeof(ImDrawVert), (const GLvoid*)((const char*)vtx_buffer + IM_OFFSETOF(ImDrawVert, pos)));
            glTexCoordPointer(2, GL_FLOAT, sizeof(ImDrawVert), (const GLvoid*)((const char*)vtx_buffer + IM_OFFSETOF(ImDrawVert, uv)));
            glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(ImDrawVert), (const GLvoid*)((const char*)vtx_buffer + IM_OFFSETOF(ImDrawVert, col)));

            for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++)
            {
                const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];

                /* Project scissor/clipping rectangles into framebuffer space */
                ImVec4 clip_rect;
                clip_rect.x = (pcmd->ClipRect.x - clip_off.x) * clip_scale.x;
                clip_rect.y = (pcmd->ClipRect.y - clip_off.y) * clip_scale.y;
                clip_rect.z = (pcmd->ClipRect.z - clip_off.x) * clip_scale.x;
                clip_rect.w = (pcmd->ClipRect.w - clip_off.y) * clip_scale.y;

                if (clip_rect.x < fb_width && clip_rect.y < fb_height && clip_rect.z >= 0.0f && clip_rect.w >= 0.0f)
                {
                    /* Apply scissor/clipping rectangle */
                    glScissor((int)clip_rect.x, (int)(fb_height - clip_rect.w), (int)(clip_rect.z - clip_rect.x), (int)(clip_rect.w - clip_rect.y));

                    /* Bind texture, Draw */
                    glBindTexture(GL_TEXTURE_2D, (GLuint)(intptr_t)pcmd->TextureId);
                    glDrawElements(GL_TRIANGLES, (GLsizei)pcmd->ElemCount, sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, idx_buffer);
                }
                idx_buffer += pcmd->ElemCount;
            }
        }


        /* Restore modified GL state */
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
        glBindTexture(GL_TEXTURE_2D, (GLuint)last_texture);
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glPopAttrib();
        glPolygonMode(GL_FRONT, (GLenum)last_polygon_mode[0]); glPolygonMode(GL_BACK, (GLenum)last_polygon_mode[1]);
        glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
        glScissor(last_scissor_box[0], last_scissor_box[1], (GLsizei)last_scissor_box[2], (GLsizei)last_scissor_box[3]);
    }

private:
    struct
    {
        GLuint mFontTexture;
        std::array<bool, GLFW_MOUSE_BUTTON_LAST> mMouseJustPressed; // detect one frame clicking
    } mImGuiData;

    Viewer& mViewer;
};

}


Viewer::Viewer()
    : mWindow{ "GL2 Viewer", 1280, 720, false, false }
    , mRender{ 0.5, 1.0, 16, true, {1.0, 1.0, 1.0} }
    , mCamera{ {0.0, 0.0}, {mWindow.width, mWindow.height}, 1.0f }
    , mUI{ false, false }
    , mFPSCounter{ 0.0, 0.0, 0, 0.0 }
{
    Window::initializedGLFW();
}

Viewer::~Viewer()
{
    Window::shutdownGLFW();
}

void Viewer::run()
{
    /* Window Creation */
    Window::ContextAttributes context;
    context.mAPI = eClientAPI::OPENGL;
    context.mProfile = eOpenGLProfile::ANY;
    context.mVersionMajor = 2;
    context.mVersionMinor = 0;

    detail::WindowImGui window(*this, mWindow.title,
                                      renderScale() * mWindow.width,
                                      renderScale() * mWindow.height,
                                      mWindow.mHDPI, context);

    /* mac display fix */
    auto windowSize = window.size();
    mWindow.width = windowSize.x;
    mWindow.height = windowSize.y;

    window.vsync(mWindow.vsync);

    /* Timer init */
    double previousTime = 0.0;
    glfwSetTime(previousTime);
    mFPSCounter.fps = 60;
    mFPSCounter.frames = 0;

    /* Render init */
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    /* User init */
    if(mOnInit) mOnInit();

    /* Main Loop */
    while(!window.closed())
    {
        Window::pollEventsGLFW();

        glClearColor(mRender.bgColor.r, mRender.bgColor.g, mRender.bgColor.b, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        /* Timer and fps counter*/
        double currentTime = glfwGetTime();
        mFPSCounter.dt = currentTime - previousTime;
        previousTime = currentTime;

        mFPSCounter.accumTime += mFPSCounter.dt;
        mFPSCounter.frames++;
        if(mFPSCounter.accumTime >= 1.0)
        {
            mFPSCounter.fps = static_cast<double>(mFPSCounter.frames) * 0.5 + mFPSCounter.fps * 0.5;
            mFPSCounter.accumTime -= 1.0f;
            mFPSCounter.frames = 0;
        }

        /* user update */
        if(mOnUpdate) mOnUpdate(window, mFPSCounter.dt);

        /* setup render state */
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        auto size = mCamera.size * mCamera.zoom;
        glOrtho(-size.x/2, size.x/2, -size.y/2, size.y/2, 0, 1);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(-mCamera.position.x, -mCamera.position.y, 0.0f);

        /* user render */
        if(mOnDraw) mOnDraw(window, mFPSCounter.dt);

        /* ImGui Rendering */
        window.newFrameImGui(window, mFPSCounter.dt);
        if(mOnGui) mOnGui(window, mFPSCounter.dt);
        mUI.keyboardCaptured = ImGui::GetIO().WantCaptureKeyboard;
        mUI.mouseCaptured = ImGui::GetIO().WantCaptureMouse;
        window.renderImGui();

        window.swapBuffers();
    }

    /* User shutdown */
    if(mOnShutdown) mOnShutdown();

    /* clean-up opengl textures */
    for(auto id : mTextureStorage)
    {
        glDeleteTextures(1, &id);
    }
    mTextureStorage.clear();
}

void Viewer::onInit(const Viewer::SetupCallback &initCB)
{
    mOnInit = initCB;
}

void Viewer::onShutdown(const Viewer::SetupCallback &shutdownCB)
{
    mOnShutdown = shutdownCB;
}

void Viewer::onKey(const Viewer::KeyboardCallback &keyCB)
{
    mOnKey = keyCB;
}

void Viewer::onMouseButton(const Viewer::MouseCallback &mouseCB)
{
    mOnMouseButton = mouseCB;
}

void Viewer::onUpdate(const Viewer::UpdateCallback &updateCB)
{
    mOnUpdate = updateCB;
}

void Viewer::onDraw(const Viewer::DrawCallback &drawCB)
{
    mOnDraw = drawCB;
}

void Viewer::onGui(const Viewer::GuiCallback &guiCB)
{
    mOnGui = guiCB;
}

unsigned int Viewer::loadTexture(const std::string& path)
{
    int width = 0, height = 0, components = 0;

    /* flip image to match opengl's texture coordinates */
    stbi_set_flip_vertically_on_load(true);

    /* load image */
    unsigned char* data = stbi_load(path.c_str(), &width, &height, &components, 4);
    if(data == nullptr)
    {
        std::cerr << "[Texture] couldn't load image file " << path << std::endl;
        std::cerr.flush();
        throw std::runtime_error("[Texture] couldn't load image file " + path);
    }

    GLuint id = 0;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    stbi_image_free(data);
    glBindTexture(GL_TEXTURE_2D, 0);

    mTextureStorage.push_back(id);

    return id;
}

unsigned int Viewer::createTexture(unsigned int width, unsigned int height, const unsigned char* data)
{
    GLuint id = 0;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

    glBindTexture(GL_TEXTURE_2D, 0);

    mTextureStorage.push_back(id);

    return id;
}

void Viewer::resizeTexture(unsigned int id, unsigned int width, unsigned int height, const unsigned char *data)
{
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, 0);
}

void Viewer::updateTexture(unsigned int id, unsigned int width, unsigned int height, const unsigned char *data)
{
    glBindTexture(GL_TEXTURE_2D, id);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, 0);
}

void Viewer::drawLine(const glm::vec2 &p0, const glm::vec2 &p1, const glm::vec4 &color)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glBegin(GL_LINES);
        glColor4f(color.r, color.g, color.b, color.a);
        glVertex2f(p0.x, p0.y);
        glVertex2f(p1.x, p1.y);
    glEnd();
}

void Viewer::drawPoint(const glm::vec2 &pos, const glm::vec4 &color)
{
    float renderscale = renderScale();
    glPointSize(renderscale * 2.0f * mRender.pointRadius);

    glBegin(GL_POINTS);
        glColor4f(color.r, color.g, color.b, color.a);
        glVertex2f(pos.x, pos.y);
    glEnd();
}

void Viewer::drawCircle(const glm::vec2 &position, float radius, glm::vec4 &color)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_TRIANGLE_FAN);
        glColor4f(color.r, color.g, color.b, color.a);

        glVertex2f(position.x, position.y);
        for(unsigned int i = 0; i <= mRender.circleVertices; i++)
        {
            glVertex2f(position.x + radius * glm::cos(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices),
                       position.y + radius * glm::sin(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices));
        }
    glEnd();

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Viewer::drawCircleOutline(const glm::vec2 &position, float radius, glm::vec4 &color)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);

    glBegin(GL_LINE_STRIP);
    {
        glColor4f(color.r, color.g, color.b, color.a);

        for(unsigned int i = 0; i <= mRender.circleVertices; i++)
        {
            glVertex2f(position.x + radius * glm::cos(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices),
                       position.y + radius * glm::sin(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices));
        }
    }
    glEnd();
}

void Viewer::drawQuad(const glm::vec2 &p0, const glm::vec2 &p1, const glm::vec2 &p2, const glm::vec2 &p3, const glm::vec4 &color)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_QUADS);
        glColor4f(color.r, color.g, color.b, color.a);
        glVertex2f(p0.x, p0.y);
        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
        glVertex2f(p3.x, p3.y);
    glEnd();

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Viewer::drawQuad2DTextured(unsigned int textureID, const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glShadeModel(GL_SMOOTH);

    glBegin(GL_QUADS);
        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex2f(p0.x, p0.y);
        glTexCoord2f(1.0f, 0.0f); glVertex2f(p1.x, p1.y);
        glTexCoord2f(1.0f, 1.0f); glVertex2f(p2.x, p2.y);
        glTexCoord2f(0.0f, 1.0f); glVertex2f(p3.x, p3.y);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Viewer::drawBoundary(const glm::vec2 &pos, const glm::vec2 &normal, float size, float pseudoShadingSize)
{
    float renderscale = renderScale();

    glLineWidth(renderscale * mRender.lineWidth);
    glColor3f(0.5, 0.5, 0.5);

    glm::vec2 line = glm::normalize(glm::vec2{normal.y, -normal.x});
    glm::vec2 shadedNormal = glm::rotate(normal, glm::pi<float>() / 4.0f);

    glm::vec2 start = (pos - (size / 2.0f) * line);
    glm::vec2 end =  (pos + (size / 2.0f) * line);

    glBegin(GL_LINES);
    glVertex2f(start.x, start.y);
    glVertex2f(end.x, end.y);

    if(pseudoShadingSize > 0.0f)
    {
        for(int i = 1; i < static_cast<int>(size / pseudoShadingSize); i++)
        {
            glVertex2f(start.x + (i*pseudoShadingSize) * line.x, start.y + (i*pseudoShadingSize) * line.y);

            glVertex2f((-pseudoShadingSize*shadedNormal).x + start.x + (i*pseudoShadingSize) * line.x,
                       (-pseudoShadingSize*shadedNormal).y + start.y + (i*pseudoShadingSize) * line.y);
        }
    }
    glEnd();
}

float Viewer::renderScale() const
{
    return (mWindow.mHDPI) ? 2.0f : 1.0f;
}

double Viewer::fps() const
{
    return mFPSCounter.fps;
}

double Viewer::dt() const
{
    return mFPSCounter.dt;
}

glm::vec2 Viewer::worldSpacePosition(const glm::dvec2 &windowPos)
{
    auto size = mCamera.size * mCamera.zoom;
    auto pos = glm::unProject(glm::vec3(windowPos.x, mWindow.height - windowPos.y, 0.0), glm::translate(-glm::vec3(mCamera.position, 0.0)),
                   glm::ortho(-size.x/2, size.x/2, -size.y/2, size.y/2, 0.0f, 1.0f), glm::vec4(0, 0, mWindow.width, mWindow.height));

    return glm::vec2(pos);
}
