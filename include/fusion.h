#pragma once

#include "bin.h"

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

typedef robin_hood::unordered_map<std::pair<size_t, size_t>,
                                  uint32_t, pair_hash>
    UMAP_PAIR_K;

// define bin_read_map
extern UMAP_PAIR_K bin_read_map;

// detect fusion
void Detect_fusion(h_Bins &, h_Reads &, uint32_t, uint32_t, int);