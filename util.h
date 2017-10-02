#pragma once

#include <cstddef>
#include <set>
#include <vector>
#include <numeric>
#include <algorithm>

using std::size_t;

template<typename T1, typename T2>
size_t index(const std::vector<T1> &v, const T2 &a)
{
    return std::distance(v.begin(), std::find(v.begin(), v.end(), a));
}

template<typename T>
std::vector<size_t> order(const std::vector<T> &v)
{
    std::vector<size_t> z(v.size());
    std::iota(z.begin(), z.end(), size_t(0));
    std::sort(z.begin(), z.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });
    return z;
}

template<typename T>
std::vector<T> unique(std::vector<T> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

template<typename T>
std::vector<T> stable_unique(std::vector<T> v)
{
    std::set<T> seen;

    auto last = v.begin();
    for (auto itr = v.begin(); itr != v.end(); ++itr) {
        if (seen.insert(*itr).second) {
            if (last != itr)
                *last = *itr;
            ++last;
        }
    }

    v.erase(last, v.end());

    return v;
}

template<typename T1, typename T2>
std::vector<T1> subset(const std::vector<T1> &v, const std::vector<T2> &idx)
{
    std::vector<T1> z;
    z.reserve(idx.size());
    for (auto i : idx)
        z.push_back(v[i]);
    return z;
}
