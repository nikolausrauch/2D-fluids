#include "particlehelper.h"

#include <algorithm>
#include <ctime>

#include <glm/gtx/norm.hpp>

#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))

namespace helper
{

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2 &min, const glm::vec2 &max, float jitterNoise)
{
    std::srand(std::time(nullptr));

    std::vector<glm::vec2> positions;
    positions.reserve(256);

    for(float x = min.x; x <= max.x; x += spacing)
    {
        for(float y = min.y; y <= max.y; y += spacing)
        {

            glm::vec2 pos = {x, y};

            pos.x += ((RAND_0_1 <= 0.5) ? -1.0 : 1.0) * jitterNoise * spacing * RAND_0_1;
            pos.y += ((RAND_0_1 <= 0.5) ? -1.0 : 1.0) * jitterNoise * spacing * RAND_0_1;

            positions.emplace_back(pos);
        }
    }

    return positions;
}

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& center, float radius, float jitterNoise)
{
    auto samples = randomPositions(spacing, center - glm::vec2{radius, radius}, center + glm::vec2{radius, radius}, jitterNoise);

    samples.erase(std::remove_if(samples.begin(), samples.end(),
                  [&](const glm::vec2& p){ return glm::length2(p - center) > radius*radius; }),
                  samples.end());

    return samples;
}


}
