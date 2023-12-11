#pragma once

#include <algorithm>
#include <vector>
#include <glm/glm.hpp>

template<typename T, const glm::vec2&(*Pos)(const T&)>
struct HashGrid
{
    HashGrid(unsigned int max_entries, unsigned int size = 512*512, unsigned int max_neighbors = 32)
    {
        mHashgrid.resize( size );
        std::for_each(mHashgrid.begin(), mHashgrid.end(), [&](auto& n){ n.reserve( max_neighbors ); });

        mNeighbors.resize(max_entries);
        std::for_each(mNeighbors.begin(), mNeighbors.end(), [&](auto& n){ n.reserve( max_neighbors ); });
    }

    void clear()
    {
        #pragma omp parallel for schedule(static)
        for(auto& cell : mHashgrid)
        {
            cell.clear();
        }
    }

    void fillGrid(const std::vector<T>& data, float radius)
    {
        for(unsigned int i = 0; i < data.size(); i++)
        {
            const auto& p = data[i];
            glm::ivec2 bucket = glm::floor(Pos(p) / radius);

            unsigned int hash = ( (73856093 * bucket.x) ^ (19349663 * bucket.y) );
            mHashgrid[ hash % mHashgrid.size() ].emplace_back(i);
        }
    }

    void updateNeighbour(const T &p, std::vector<T>& data, float radius)
    {
        updateNeighbour(std::distance(data.data(), &p), data, radius);
    }

    void updateNeighbour(unsigned int i, std::vector<T>& data, float radius)
    {
        const auto& p = data[i];
        glm::ivec2 bucket = glm::floor(Pos(p) / radius);

        mNeighbors[i].clear();
        for(int x = -1; x <= 1; x++)
        {
            for(int y = -1; y <= 1; y++)
            {
                glm::ivec2 search(bucket.x + x,bucket.y + y);
                unsigned int hash = ( (73856093 * search.x) ^ (19349663 * search.y) );
                for(unsigned int j : mHashgrid[ hash % mHashgrid.size() ])
                {
                    if(i == j) continue;

                    const auto& pj = data[j];
                    if(glm::length(Pos(p) - Pos(pj)) < radius)
                    {
                        mNeighbors[i].emplace_back(j);
                    }
                }
            }
        }
    }

    void updateGhostNeighbour(unsigned int i, const glm::vec2& pos, std::vector<T>& particles, float radius)
    {
        glm::ivec2 bucket = glm::floor(pos / radius);

        mNeighbors[i].clear();
        for(int x = -1; x <= 1; x++)
        {
            for(int y = -1; y <= 1; y++)
            {
                glm::ivec2 search(bucket.x + x,bucket.y + y);
                unsigned int hash = ( (73856093 * search.x) ^ (19349663 * search.y) );
                for(unsigned int j : mHashgrid[ hash % mHashgrid.size() ])
                {
                    const auto& pj = particles[j];
                    if(glm::length(pos - pj.position) < radius)
                    {
                        mNeighbors[i].emplace_back(j);
                    }
                }
            }
        }
    }

    void fillNeighbors(std::vector<T>& data, float radius)
    {
        #pragma omp parallel for schedule(static)
        for(unsigned int i = 0; i < data.size(); i++)
        {
            updateNeighbour(i, data, radius);
        }
    }

    const std::vector<unsigned int>& neighbors(unsigned int p_index)
    {
        return mNeighbors[p_index];
    }

    const std::vector<unsigned int>& neighbors(const T& p, const std::vector<T>& data)
    {
        return neighbors(std::distance(data.data(), &p));
    }

protected:
    std::vector< std::vector<unsigned int> > mNeighbors;
    std::vector< std::vector<unsigned int> > mHashgrid;
};
