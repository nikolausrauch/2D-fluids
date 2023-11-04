template<class InputIt, class UnaryFunction>
void Viewer::drawPoints(InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glPointSize(renderscale * 2.0f * mRender.pointRadius);
    glm::vec2 coord;
    glm::vec4 color;

    glBegin(GL_POINTS);
    for (; first != last; ++first)
    {
        if(!f(*first, coord, color)) continue;

        glColor4f(color.r, color.g, color.b, color.a);
        glVertex2f(coord.x, coord.y);
    }
    glEnd();
}

template<class InputIt, class UnaryFunction>
void Viewer::drawLines(InputIt first, InputIt last, UnaryFunction f)
{    
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 start;
    glm::vec2 end;
    glm::vec4 color;

    glBegin(GL_LINES);
    for (; first != last; ++first)
    {
        if(!f(*first, start, end, color)) continue;

        glColor4f(color.r, color.g, color.b, color.a);
        glVertex2f(start.x, start.y);
        glVertex2f(end.x, end.y);
    }
    glEnd();
}

template<class InputIt, class UnaryFunction>
void Viewer::drawCircles(InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 position = {0.0, 0.0};
    glm::vec4 color = {0.0, 0.0, 0.0, 1.0};
    float radius = 1.0f;

    std::vector<glm::vec2> vertices(mRender.circleVertices);
    for(unsigned int i = 0; i <= mRender.circleVertices; i++)
    {
        vertices[i] =
        {
            radius * glm::cos(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices),
            radius * glm::sin(2.0*glm::pi<float>() * static_cast<float>(i) / mRender.circleVertices)
        };
    }


    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    for (; first != last; ++first)
    {
        if(!f(*first, position, radius, color)) continue;

        glBegin(GL_TRIANGLE_FAN);
        glColor4f(color.r, color.g, color.b, color.a);

        glVertex2f(position.x, position.y);

        for(unsigned int i = 0; i <= mRender.circleVertices; i++)
        {
            glVertex2f(position.x + radius * vertices[i].x,
                       position.y + radius * vertices[i].y);
        }
        glEnd();
    }

    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

template<class InputIt, class UnaryFunction>
void Viewer::drawTriangles(InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 p0;
    glm::vec2 p1;
    glm::vec2 p2;
    glm::vec4 color;

    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    glBegin(GL_TRIANGLES);
    for (; first != last; ++first)
    {
        if(!f(*first, p0, p1, p2, color)) continue;

        glColor4f(color.r, color.g, color.b, color.a);

        glVertex2f(p0.x, p0.y);
        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
    }
    glEnd();

    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

template<class InputIt, class UnaryFunction>
void Viewer::drawTrianglesTextured(unsigned int textureID, InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 p0;
    glm::vec2 p1;
    glm::vec2 p2;
    glm::vec2 uv0;
    glm::vec2 uv1;
    glm::vec2 uv2;

    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glBegin(GL_TRIANGLES);
    for (; first != last; ++first)
    {
        if(!f(*first, p0, p1, p2, uv0, uv1, uv2)) continue;

        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

        glTexCoord2f(uv0.x, uv0.y); glVertex2f(p0.x, p0.y);
        glTexCoord2f(uv1.x, uv1.y); glVertex2f(p1.x, p1.y);
        glTexCoord2f(uv2.x, uv2.y); glVertex2f(p2.x, p2.y);
    }
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

    if(mRender.wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

template<class InputIt, class UnaryFunction>
void Viewer::drawQuads2D(InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 p0;
    glm::vec2 p1;
    glm::vec2 p2;
    glm::vec2 p3;
    glm::vec4 color;

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_QUADS);
    for (; first != last; ++first)
    {
        if (!f(*first, p0, p1, p2, p3, color)) continue;

        glColor4f(color.r, color.g, color.b, color.a);

        glVertex2f(p0.x, p0.y);
        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
        glVertex2f(p3.x, p3.y);

    }
    glEnd();

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

template<class InputIt, class UnaryFunction>
void Viewer::drawQuads2DTextured(unsigned int textureID, InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    glm::vec2 p0;
    glm::vec2 p1;
    glm::vec2 p2;
    glm::vec2 p3;

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glBegin(GL_QUADS);
    for (; first != last; ++first)
    {
        if (!f(*first, p0, p1, p2, p3, textureID)) continue;

        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

        glTexCoord2f(0.0f, 0.0f); glVertex2f(p0.x, p0.y);
        glTexCoord2f(1.0f, 0.0f); glVertex2f(p1.x, p1.y);
        glTexCoord2f(1.0f, 1.0f); glVertex2f(p2.x, p2.y);
        glTexCoord2f(0.0f, 1.0f); glVertex2f(p3.x, p3.y);

    }
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (mRender.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

template<class InputIt, class UnaryFunction>
void Viewer::drawOutline(InputIt first, InputIt last, UnaryFunction f)
{
    float renderscale = renderScale();
    glLineWidth(renderscale * mRender.lineWidth);
    std::vector<glm::vec2> vertices;
    glm::vec4 color;

    for (; first != last; ++first)
    {
        vertices.clear();
        if(!f(*first, vertices, color)) continue;

        glBegin(GL_LINE_LOOP);

        glColor4f(color.r, color.g, color.b, color.a);

        for(auto& v : vertices)
        {
            glVertex2f(v.x, v.y);
        }

        glEnd();
    }
}
