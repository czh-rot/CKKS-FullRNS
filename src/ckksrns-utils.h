//
// Created by EYx on 2023/5/4.
//

#ifndef FRNS_BOOTSTRAP_CKKSRNS_UTILS_H
#define FRNS_BOOTSTRAP_CKKSRNS_UTILS_H

#include "Common.h"

/**
 * Computes parameters to ensure the encoding and decoding computations take exactly the
 * specified number of levels. More specifically, it returns a vector than contains
 * layers (the number of layers to collapse in one level), rows (how many such levels),
 * rem (the number of layers remaining to be collapsed in one level)
 *
 * @param logSlots the base 2 logarithm of the number of slots.
 * @param budget the allocated level budget for the computation.
 */
std::vector<uint32_t> SelectLayers(uint32_t logSlots, uint32_t budget = 4);

/**
 * Computes all parameters needed for the homomorphic encoding and decoding in the bootstrapping
 * operation and returns them as a vector. The returned vector's data can be accessed using
 * enum'ed indices from CKKS_BOOT_PARAMS that are defined below.
 *
 * @param slots number of slots.
 * @param levelBudget the allocated level budget for the computation.
 * @param dim1 the value for the inner dimension in the baby-step giant-step strategy
 * @return vector with parameters for the homomorphic encoding and decoding in bootstrapping
 */
std::vector<int32_t> GetCollapsedFFTParams(uint32_t slots, uint32_t levelBudget = 4, uint32_t dim1 = 0);

uint32_t ReduceRotation(int32_t index, uint32_t slots);

namespace CKKS_BOOT_PARAMS {
/**
   * Enums representing indices for the vector returned by GetCollapsedFFTParams()
   */
    enum {
        LEVEL_BUDGET,  // the level budget
        LAYERS_COLL,   // the number of layers to collapse in one level
        LAYERS_REM,  // the number of layers remaining to be collapsed in one level to have exactly the number of levels specified in the level budget
        NUM_ROTATIONS,      // the number of rotations in one level
        BABY_STEP,          // the baby step in the baby-step giant-step strategy
        GIANT_STEP,         // the giant step in the baby-step giant-step strategy
        NUM_ROTATIONS_REM,  // the number of rotations in the remaining level
        BABY_STEP_REM,      // the baby step in the baby-step giant-step strategy for the remaining level
        GIANT_STEP_REM,     // the giant step in the baby-step giant-step strategy for the remaining level
        TOTAL_ELEMENTS      // total number of elements in the vector
    };
}  // namespace CKKS_BOOT_PARAMS

#endif //FRNS_BOOTSTRAP_CKKSRNS_UTILS_H
