#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/norm.hpp>

namespace kernel
{

/* kernel radius w.r.t. particle radius */
static inline constexpr float radiusMultiplier = 4.0f;


struct std
{
    std(float h) : h(h), h_2(h*h), h_3(h_2*h), h_4(h_3*h) {}

    float operator() (float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = 1.0f - d_2 / h_2;
        return 4.0f / (glm::pi<float>() * h_2) * x * x * x;
    }

    float firstDerivative(float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = 1.0f - d_2 / h_2;
        return -24.0f * d / (glm::pi<float>() * h_4) * x * x;
    }

    float secondDerivative(float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = d_2 / h_2;
        return 24.0f / (glm::pi<float>() * h_4) * (1.0 - x) * (5 * x - 1);
    }

    float secondDerivative(const glm::vec2& x_i, const glm::vec2& x_j) const
    {
        float d = glm::length(x_i - x_j);
        return secondDerivative(d, d*d);
    }

    glm::vec2 gradient(const glm::vec2& x_i, const glm::vec2& x_j) const
    {
        auto dir = x_j - x_i;
        auto d = glm::length(dir);

        return (d != 0.0f) ? - firstDerivative(d, d*d) * (dir / d) : glm::vec2{0.0f, 0.0f};
    }

private:
    float h;
    float h_2;
    float h_3;
    float h_4;
};

struct spiky
{
    spiky(float h) : h(h), h_2(h*h), h_3(h_2*h), h_4(h_3*h) {}

    float operator() (float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = 1.0f - d / h;
        return 10.0f / (glm::pi<float>() * h_2) * x * x * x;
    }

    float firstDerivative(float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = 1.0f - d/h;
        return -30.0f / (glm::pi<float>() * h_3) * x * x;
    }

    float secondDerivative(float d, float d_2) const
    {
        if(d >= h) return 0.0f;

        float x = 1.0f - d/h;
        return 60.0f / (glm::pi<float>() * h_4) * x;
    }

    glm::vec2 gradient(const glm::vec2& x_i, const glm::vec2& x_j) const
    {
        auto dir = x_j - x_i;
        auto d = glm::length(dir);

        return (d != 0.0f) ? - firstDerivative(d, d*d) * (dir / d) : glm::vec2{0.0f, 0.0f};
    }

private:
    float h;
    float h_2;
    float h_3;
    float h_4;
};

}
