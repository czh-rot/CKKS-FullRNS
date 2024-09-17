//
// Created by EYx on 2023/4/25.
//

#ifndef FRNS_BOOTSTRAP_MYTESTROTATION_H
#define FRNS_BOOTSTRAP_MYTESTROTATION_H
#include "Context.h"
#include "SecretKey.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"
#include "TimeUtils.h"
void testAutomorphismTransform();
void testPrecomputeAutoMap();

void testHRotate(long logN, long L, long logp, long dnum, long logSlots);

void testFastRotate_demo(long logN, long L, long logp, long dnum, long rotSlots, long logSlots);
void testEvalFastRotation(long logN, long L, long logp, long dnum, long rotSlots, long logSlots);

void testEvalFastRotationExt(long logN, long L, long logp, long dnum, long rotSlots, long logSlots);

void testPartialSum(long logN, long L, long logp, long dnum, long logSlots);

void testConjugate_demo(long logN, long L, long logp, long dnum, long logSlots);
#endif //FRNS_BOOTSTRAP_MYTESTROTATION_H
