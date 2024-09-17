//
// Created by EYx on 2023/5/5.
//

#ifndef FRNS_BOOTSTRAP_CONSTANTS_H
#define FRNS_BOOTSTRAP_CONSTANTS_H

#define NORMAL_CIPHER_SIZE 2

enum LargeScalingFactorConstants {
    MAX_BITS_IN_WORD = 62,
    MAX_LOG_STEP     = 60,
};


/**
 * @brief  BASE_NUM_LEVELS_TO_DROP is the most common value for levels/towers to drop (do not make it a default argument
 * as default arguments work differently for virtual functions)
 */
// TODO (dsuponit): remove BASE_NUM_LEVELS_TO_DROP
enum {
    BASE_NUM_LEVELS_TO_DROP = 1,
};

/**
 * @brief Lists all modes for RLWE schemes, such as BGV and BFV
 */
enum SecretKeyDist {
    GAUSSIAN        = 0,
    UNIFORM_TERNARY = 1,
    SPARSE_TERNARY  = 2,
};

enum ScalingTechnique {
    FIXEDMANUAL = 0,
    FIXEDAUTO,
    FLEXIBLEAUTO,
    FLEXIBLEAUTOEXT,
    NORESCALE,
    INVALID_RS_TECHNIQUE,  // TODO (dsuponit): make this the first value
};

#endif //FRNS_BOOTSTRAP_CONSTANTS_H
