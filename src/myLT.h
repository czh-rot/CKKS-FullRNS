//
// Created by EYx on 2023/5/4.
//

#ifndef INC_0503_C2S_DEMO_MYLT_H
#define INC_0503_C2S_DEMO_MYLT_H

void EvalLinearTransform(Ciphertext &result,
                         Plaintext *A, uint32_t A_len,
                         const Ciphertext &ct,
                         Scheme scheme);

void EvalSlotsToCoeffs(Ciphertext &result,
                       Plaintext** A,uint32_t A_len,
                       Ciphertext ctxt,
                       Scheme scheme);

void EvalCoeffsToSlots(Ciphertext &result,
                       Plaintext** A,uint32_t A_len,
                       Ciphertext ctxt,
                       Scheme scheme);

#endif //INC_0503_C2S_DEMO_MYLT_H
