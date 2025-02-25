#ifndef SUPERLU_WRAPPER_H
#define SUPERLU_WRAPPER_H

#include <stdlib.h>
#include <slu_ddefs.h>

#ifdef __cplusplus
extern "C" {
#endif

// Function wrappers for SuperLU functions (double precision only)
void wrapper_dgssv(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info);
void wrapper_dgssvx(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C, SuperMatrix* L, SuperMatrix* U, void* work, int lwork, SuperMatrix* B, SuperMatrix* X, double* rpg, double* rcond, double* ferr, double* berr, GlobalLU_t* Glu, mem_usage_t* mem_usage, SuperLUStat_t* stat, int* info);

void wrapper_StatInit(SuperLUStat_t* stat);
void wrapper_StatFree(SuperLUStat_t* stat);
void wrapper_set_default_options(superlu_options_t* options);
void* wrapper_superlu_malloc(size_t size);
void wrapper_superlu_free(void* ptr);
void wrapper_Destroy_SuperNode_Matrix(SuperMatrix* A);
void wrapper_Destroy_CompCol_Matrix(SuperMatrix* A);
void wrapper_Destroy_CompCol_Permuted(SuperMatrix* A);
void wrapper_Destroy_Dense_Matrix(SuperMatrix* A);

#ifdef __cplusplus
}
#endif

#endif // SUPERLU_WRAPPER_H 