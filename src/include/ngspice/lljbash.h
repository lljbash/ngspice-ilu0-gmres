#pragma once

#include <stdint.h>

typedef struct {
    int size;
    int max_nnz;
    int* row_ptr;
    int* col_idx;
    double* value;
} LLJBASH_CsrMatrix;

typedef struct {
    double tolerance;
    int krylov_subspace_dimension;
    int max_iterations;
} LLJBASH_GmresParameters;

typedef struct {
    double total_init_time;
    double total_fgmr_time;
    double total_mv_time;
    double total_precon_time;
} LLJBASH_GmresStat;

typedef struct MatrixFrame SMPmatrix;
typedef struct CKTcircuit CKTcircuit;
typedef struct MatrixElement* ElementPtr;

typedef struct {
    void* ilu;
    void* gmres;
    SMPmatrix* smp;
    LLJBASH_CsrMatrix csr;
    int need_setup;
    ElementPtr* element_mapping;
    int first_ilu;
    double* intermediate;

    struct {
        int dc_finished;
        double direct_dc_time;
        double direct_transient_time;
        double ilu_setup_time;
        double total_precon_time;
        double total_gmres_time;
    } stat;
} LLJBASH_Solver;
extern LLJBASH_Solver* lljbash_this;

extern struct LLJBASH_SolverFunctions {
    void (*Init)(LLJBASH_Solver*);
    void (*Free)(LLJBASH_Solver*);
    void (*ImportMatrix)(LLJBASH_Solver*, SMPmatrix*);
    void (*InitPreconditoner)(LLJBASH_Solver*);
    int (*GmresSolve)(LLJBASH_Solver*, double*);
    int (*NIiter)(LLJBASH_Solver*, CKTcircuit*, int);
    int (*CKTop)(LLJBASH_Solver*, CKTcircuit*, long, long, int);
} lljbash_solver;

extern struct LLJBASH_Functions {
    void* dlhandler;

    void (*SetupCsrMatrix)(LLJBASH_CsrMatrix*, int, int);
    void (*DestroyCsrMatrix)(LLJBASH_CsrMatrix*);
    void (*CopyCsrMatrix)(LLJBASH_CsrMatrix*, const LLJBASH_CsrMatrix*);
    const LLJBASH_CsrMatrix* CSR_MATRIX_DEFAULT;

    void* (*IluSolverCreate)(int);
    void (*IluSolverDestroy)(void*);
    LLJBASH_CsrMatrix* (*IluSolverGetMatrix)(void*);
    int (*IluSolverSetup)(void*);
    int (*IluSolverFactorize)(void*);

    void* (*GmresCreate)(void);
    void (*GmresDestroy)(void*);
    LLJBASH_GmresParameters* (*GmresGetParameters)(void*);
    LLJBASH_GmresStat* (*GmresGetStat)(void*);
    void (*GmresSetPreconditioner)(void*, void*);
    int (*GmresSolve)(void*, const LLJBASH_CsrMatrix*, const double*, double*, int*);
} lljbash;

#define LLJBASH_CPU_FREQUENCY 2.4e9

extern uint64_t LLJBASH_cycles;

inline uint64_t LLJBASH_Rdtsc(void) {
    uint32_t lo, hi;
    __asm__ __volatile__ (
            "xorl %%eax, %%eax\n"
            "cpuid\n"
            "rdtsc\n"
            : "=a" (lo), "=d" (hi)
            :
            : "%ebx", "%ecx");
    return (uint64_t) hi << 32 | lo;
}

inline void LLJBASH_Tic(void) {
    LLJBASH_cycles = LLJBASH_Rdtsc();
}

inline double LLJBASH_Toc(void) {
    return (double)(LLJBASH_Rdtsc() - LLJBASH_cycles) / LLJBASH_CPU_FREQUENCY;
}

#undef CPU_FREQUENCY
