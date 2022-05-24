#pragma once

typedef struct {
    int size;
    long max_nnz;
    long* row_ptr;
    int* col_idx;
    double* value;
} LLJBASH_CsrMatrix;

typedef struct {
    double tolerance;
    int krylov_subspace_dimension;
    int max_iterations;
} LLJBASH_GmresParameters;

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
} LLJBASH_Solver;

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

    void (*SetupCsrMatrix)(LLJBASH_CsrMatrix*, int, long);
    void (*DestroyCsrMatrix)(LLJBASH_CsrMatrix*);
    void (*CopyCsrMatrix)(LLJBASH_CsrMatrix*, const LLJBASH_CsrMatrix*);
    const LLJBASH_CsrMatrix* CSR_MATRIX_DEFAULT;

    void* (*IluSolverCreate)(int);
    void (*IluSolverDestroy)(void*);
    LLJBASH_CsrMatrix* (*IluSolverGetMatrix)(void*);
    int (*IluSolverFactorize)(void*, int);

    void* (*GmresCreate)(void);
    void (*GmresDestroy)(void*);
    LLJBASH_GmresParameters* (*GmresGetParameters)(void*);
    void (*GmresSetPreconditioner)(void*, void*);
    int (*GmresSolve)(void*, const LLJBASH_CsrMatrix*, const double*, double*, int*);
} lljbash;
