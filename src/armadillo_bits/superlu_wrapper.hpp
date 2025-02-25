// This file is part of the Armadillo C++ library
// It is provided as a convenience for users of Armadillo
// The content below is a copy of the relevant part from Armadillo's source code
// to provide the SuperLU wrapper functions that Armadillo needs

#ifndef ARMA_SUPERLU_WRAPPER
#define ARMA_SUPERLU_WRAPPER

#include "../superlu_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

// Function prototypes for SuperLU wrappers
void wrapper_sgssv(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperLUStat_t*, int*);
void wrapper_dgssv(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperLUStat_t*, int*);
void wrapper_cgssv(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperLUStat_t*, int*);
void wrapper_zgssv(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperMatrix*, superlu::SuperLUStat_t*, int*);

void wrapper_sgssvx(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, int*, char*, float*, float*, superlu::SuperMatrix*, superlu::SuperMatrix*, void*, int, superlu::SuperMatrix*, superlu::SuperMatrix*, float*, float*, float*, float*, superlu::GlobalLU_t*, superlu::mem_usage_t*, superlu::SuperLUStat_t*, int*);
void wrapper_dgssvx(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, int*, char*, double*, double*, superlu::SuperMatrix*, superlu::SuperMatrix*, void*, int, superlu::SuperMatrix*, superlu::SuperMatrix*, double*, double*, double*, double*, superlu::GlobalLU_t*, superlu::mem_usage_t*, superlu::SuperLUStat_t*, int*);
void wrapper_cgssvx(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, int*, char*, float*, float*, superlu::SuperMatrix*, superlu::SuperMatrix*, void*, int, superlu::SuperMatrix*, superlu::SuperMatrix*, float*, float*, float*, float*, superlu::GlobalLU_t*, superlu::mem_usage_t*, superlu::SuperLUStat_t*, int*);
void wrapper_zgssvx(superlu::superlu_options_t*, superlu::SuperMatrix*, int*, int*, int*, char*, double*, double*, superlu::SuperMatrix*, superlu::SuperMatrix*, void*, int, superlu::SuperMatrix*, superlu::SuperMatrix*, double*, double*, double*, double*, superlu::GlobalLU_t*, superlu::mem_usage_t*, superlu::SuperLUStat_t*, int*);

void wrapper_StatInit(superlu::SuperLUStat_t*);
void wrapper_StatFree(superlu::SuperLUStat_t*);
void wrapper_set_default_options(superlu::superlu_options_t*);

void* wrapper_superlu_malloc(size_t);
void wrapper_superlu_free(void*);

void wrapper_Destroy_SuperNode_Matrix(superlu::SuperMatrix*);
void wrapper_Destroy_CompCol_Matrix(superlu::SuperMatrix*);
void wrapper_Destroy_CompCol_Permuted(superlu::SuperMatrix*);
void wrapper_Destroy_Dense_Matrix(superlu::SuperMatrix*);

// Implementations of the wrapper functions
void wrapper_sgssv(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, superlu::SuperMatrix* L, superlu::SuperMatrix* U, superlu::SuperMatrix* B, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_sgssv((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, (SuperMatrix*)L, (SuperMatrix*)U, (SuperMatrix*)B, (SuperLUStat_t*)stat, info);
}

void wrapper_dgssv(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, superlu::SuperMatrix* L, superlu::SuperMatrix* U, superlu::SuperMatrix* B, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_dgssv((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, (SuperMatrix*)L, (SuperMatrix*)U, (SuperMatrix*)B, (SuperLUStat_t*)stat, info);
}

void wrapper_cgssv(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, superlu::SuperMatrix* L, superlu::SuperMatrix* U, superlu::SuperMatrix* B, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_cgssv((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, (SuperMatrix*)L, (SuperMatrix*)U, (SuperMatrix*)B, (SuperLUStat_t*)stat, info);
}

void wrapper_zgssv(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, superlu::SuperMatrix* L, superlu::SuperMatrix* U, superlu::SuperMatrix* B, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_zgssv((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, (SuperMatrix*)L, (SuperMatrix*)U, (SuperMatrix*)B, (SuperLUStat_t*)stat, info);
}

void wrapper_sgssvx(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C, superlu::SuperMatrix* L, superlu::SuperMatrix* U, void* work, int lwork, superlu::SuperMatrix* B, superlu::SuperMatrix* X, float* rpg, float* rcond, float* ferr, float* berr, superlu::GlobalLU_t* Glu, superlu::mem_usage_t* mem_usage, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_sgssvx((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, etree, equed, R, C, (SuperMatrix*)L, (SuperMatrix*)U, work, lwork, (SuperMatrix*)B, (SuperMatrix*)X, rpg, rcond, ferr, berr, (GlobalLU_t*)Glu, (mem_usage_t*)mem_usage, (SuperLUStat_t*)stat, info);
}

void wrapper_dgssvx(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C, superlu::SuperMatrix* L, superlu::SuperMatrix* U, void* work, int lwork, superlu::SuperMatrix* B, superlu::SuperMatrix* X, double* rpg, double* rcond, double* ferr, double* berr, superlu::GlobalLU_t* Glu, superlu::mem_usage_t* mem_usage, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_dgssvx((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, etree, equed, R, C, (SuperMatrix*)L, (SuperMatrix*)U, work, lwork, (SuperMatrix*)B, (SuperMatrix*)X, rpg, rcond, ferr, berr, (GlobalLU_t*)Glu, (mem_usage_t*)mem_usage, (SuperLUStat_t*)stat, info);
}

void wrapper_cgssvx(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C, superlu::SuperMatrix* L, superlu::SuperMatrix* U, void* work, int lwork, superlu::SuperMatrix* B, superlu::SuperMatrix* X, float* rpg, float* rcond, float* ferr, float* berr, superlu::GlobalLU_t* Glu, superlu::mem_usage_t* mem_usage, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_cgssvx((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, etree, equed, R, C, (SuperMatrix*)L, (SuperMatrix*)U, work, lwork, (SuperMatrix*)B, (SuperMatrix*)X, rpg, rcond, ferr, berr, (GlobalLU_t*)Glu, (mem_usage_t*)mem_usage, (SuperLUStat_t*)stat, info);
}

void wrapper_zgssvx(superlu::superlu_options_t* options, superlu::SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C, superlu::SuperMatrix* L, superlu::SuperMatrix* U, void* work, int lwork, superlu::SuperMatrix* B, superlu::SuperMatrix* X, double* rpg, double* rcond, double* ferr, double* berr, superlu::GlobalLU_t* Glu, superlu::mem_usage_t* mem_usage, superlu::SuperLUStat_t* stat, int* info) {
    wrapper_zgssvx((superlu_options_t*)options, (SuperMatrix*)A, perm_c, perm_r, etree, equed, R, C, (SuperMatrix*)L, (SuperMatrix*)U, work, lwork, (SuperMatrix*)B, (SuperMatrix*)X, rpg, rcond, ferr, berr, (GlobalLU_t*)Glu, (mem_usage_t*)mem_usage, (SuperLUStat_t*)stat, info);
}

void wrapper_StatInit(superlu::SuperLUStat_t* stat) {
    wrapper_StatInit((SuperLUStat_t*)stat);
}

void wrapper_StatFree(superlu::SuperLUStat_t* stat) {
    wrapper_StatFree((SuperLUStat_t*)stat);
}

void wrapper_set_default_options(superlu::superlu_options_t* options) {
    wrapper_set_default_options((superlu_options_t*)options);
}

void* wrapper_superlu_malloc(size_t size) {
    return wrapper_superlu_malloc(size);
}

void wrapper_superlu_free(void* ptr) {
    wrapper_superlu_free(ptr);
}

void wrapper_Destroy_SuperNode_Matrix(superlu::SuperMatrix* A) {
    wrapper_Destroy_SuperNode_Matrix((SuperMatrix*)A);
}

void wrapper_Destroy_CompCol_Matrix(superlu::SuperMatrix* A) {
    wrapper_Destroy_CompCol_Matrix((SuperMatrix*)A);
}

void wrapper_Destroy_CompCol_Permuted(superlu::SuperMatrix* A) {
    wrapper_Destroy_CompCol_Permuted((SuperMatrix*)A);
}

void wrapper_Destroy_Dense_Matrix(superlu::SuperMatrix* A) {
    wrapper_Destroy_Dense_Matrix((SuperMatrix*)A);
}

#ifdef __cplusplus
}
#endif

#endif 