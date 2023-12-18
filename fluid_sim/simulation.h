#pragma once

#include <glm/vec2.hpp>

#include <vector>

struct Scene
{
    enum class eType
    {
        FLUID_BODY,
        BOUNDARY
    };

    struct Box
    {
        glm::vec2 min;
        glm::vec2 max;
        eType type;
    };

    struct Circle
    {
        glm::vec2 center;
        float radius;
        eType type;
    };


    /* construct scene */
    void clear();
    void add(const Box& box);
    void add(const Circle& box);

    /* retrieve scene geometry */
    const std::vector<Box>& boxes() const;
    const std::vector<Circle>& circles() const;


private:
    std::vector<Box> mBoxes;
    std::vector<Circle> mCircles;
};


struct Simulation
{
    virtual void create(const Scene& desc) = 0;
    virtual void clear() = 0;
    virtual void update() = 0;
};
