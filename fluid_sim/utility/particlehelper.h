#pragma once

#include <glm/vec2.hpp>

#include <vector>

namespace helper
{

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& min, const glm::vec2& max, float jitterNoise = 0.0f);
std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& center, float radius, float jitterNoise = 0.0f);

}
