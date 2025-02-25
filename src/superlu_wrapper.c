#include <stdlib.h>
#include <slu_ddefs.h>

// Function wrappers for SuperLU functions (double precision only)
void wrapper_dgssv(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info) {
    dgssv(options, A, perm_c, perm_r, L, U, B, stat, info);
}

void wrapper_dgssvx(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, double* R, double* C, SuperMatrix* L, SuperMatrix* U, void* work, int lwork, SuperMatrix* B, SuperMatrix* X, double* rpg, double* rcond, double* ferr, double* berr, GlobalLU_t* Glu, mem_usage_t* mem_usage, SuperLUStat_t* stat, int* info) {
    dgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work, lwork, B, X, rpg, rcond, ferr, berr, Glu, mem_usage, stat, info);
}

void wrapper_StatInit(SuperLUStat_t* stat) {
    StatInit(stat);
}

void wrapper_StatFree(SuperLUStat_t* stat) {
    StatFree(stat);
}

void wrapper_set_default_options(superlu_options_t* options) {
    set_default_options(options);
}

void* wrapper_superlu_malloc(size_t size) {
    return superlu_malloc(size);
}

void wrapper_superlu_free(void* ptr) {
    superlu_free(ptr);
}

void wrapper_Destroy_SuperNode_Matrix(SuperMatrix* A) {
    Destroy_SuperNode_Matrix(A);
}

void wrapper_Destroy_CompCol_Matrix(SuperMatrix* A) {
    Destroy_CompCol_Matrix(A);
}

void wrapper_Destroy_CompCol_Permuted(SuperMatrix* A) {
    Destroy_CompCol_Permuted(A);
}

void wrapper_Destroy_Dense_Matrix(SuperMatrix* A) {
    Destroy_Dense_Matrix(A);
} 