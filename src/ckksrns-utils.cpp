//
// Created by EYx on 2023/5/4.
//

#include "ckksrns-utils.h"
#include "Common.h"

std::vector<uint32_t> SelectLayers(uint32_t logSlots, uint32_t budget) {
    uint32_t layers = ceil(static_cast<double>(logSlots) / budget);
    uint32_t rows   = logSlots / layers;
    uint32_t rem    = logSlots % layers;

    uint32_t dim = rows;
    if (rem != 0) {
        dim = rows + 1;
    }

    // the above choice ensures dim <= budget
    if (dim < budget) {
        layers -= 1;
        rows = logSlots / layers;
        rem  = logSlots - rows * layers;
        dim  = rows;

        if (rem != 0) {
            dim = rows + 1;
        }

        // the above choice endures dim >=budget
        while (dim != budget) {
            rows -= 1;
            rem = logSlots - rows * layers;
            dim = rows;
            if (rem != 0) {
                dim = rows + 1;
            }
        }
    }

    return {layers, rows, rem};
}

std::vector<int32_t> GetCollapsedFFTParams(uint32_t slots, uint32_t levelBudget, uint32_t dim1) {
    // Need to compute how many layers are collapsed in each of the level from the budget.
    // If there is no exact division between the maximum number of possible levels (log(slots)) and the
    // level budget, the last level will contain the remaining layers collapsed.
    int32_t layersCollapse;
    int32_t remCollapse;

    std::vector<uint32_t> dims = SelectLayers(std::log2(slots), levelBudget);
    layersCollapse             = dims[0];
    remCollapse                = dims[2];

    int32_t flagRem = 0;
    if (remCollapse == 0) {
        flagRem = 0;
    }
    else {
        flagRem = 1;
    }

    uint32_t numRotations    = (1 << (layersCollapse + 1)) - 1;
    uint32_t numRotationsRem = (1 << (remCollapse + 1)) - 1;

    // Computing the baby-step b and the giant-step g for the collapsed layers for decoding.
    int32_t b, g;
    if (dim1 == 0 || dim1 > numRotations) {
        if (numRotations > 7) {
            g = (1 << (int32_t(layersCollapse / 2) + 2));
        }
        else {
            g = (1 << (int32_t(layersCollapse / 2) + 1));
        }
    }
    else {
        g = dim1;
    }

    b            = (numRotations + 1) / g;
    int32_t bRem = 0;
    int32_t gRem = 0;

    if (flagRem) {
        if (numRotationsRem > 7) {
            gRem = (1 << (int32_t(remCollapse / 2) + 2));
        }
        else {
            gRem = (1 << (int32_t(remCollapse / 2) + 1));
        }
        bRem = (numRotationsRem + 1) / gRem;
    }

    // If this return statement changes then CKKS_BOOT_PARAMS should be altered as well
    return {int32_t(levelBudget),     layersCollapse, remCollapse, int32_t(numRotations), b, g,
            int32_t(numRotationsRem), bRem,           gRem};
}

//DongYilin add
uint32_t ReduceRotation(int32_t index, uint32_t slots) {
    int32_t islots = int32_t(slots);

    // if slots is a power of 2
    if ((slots & (slots - 1)) == 0) {
        int32_t n = std::log2(slots);
        if (index >= 0) {
            return index - ((index >> n) << n);
        }
        return index + islots + ((int32_t(std::fabs(index)) >> n) << n);
    }
    return (islots + index % islots) % islots;
}