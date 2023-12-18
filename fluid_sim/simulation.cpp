#include "simulation.h"

void Scene::clear()
{
    mBoxes.clear();
    mCircles.clear();
}

void Scene::add(const Box& box)
{
    mBoxes.emplace_back(box);
}

void Scene::add(const Circle& circle)
{
    mCircles.emplace_back(circle);
}

const std::vector<Scene::Box>& Scene::boxes() const
{
    return mBoxes;
}

const std::vector<Scene::Circle>& Scene::circles() const
{
    return mCircles;
}
