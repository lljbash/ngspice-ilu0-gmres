#pragma once

typedef struct {
    int size;
    long max_nnz;
    long* col_ptr;
    int* row_idx;
    double* value;
} LLJBASH_CscMatrix;

typedef struct {
    double tolerance;
    int krylov_subspace_dimension;
    int max_iterations;
} LLJBASH_GmresParameters;

typedef struct MatrixFrame SMPmatrix;
typedef struct CKTcircuit CKTcircuit;

typedef struct {
    void* ilu;
    void* gmres;
} LLJBASH_Solver;

extern struct LLJBASH_SolverFunctions {
    void (*Init)(LLJBASH_Solver*);
    void (*Free)(LLJBASH_Solver*);
    void (*ImportMatrix)(LLJBASH_Solver*, const SMPmatrix*);
    void (*InitPreconditoner)(LLJBASH_Solver*);
    void (*GmresSolve)(LLJBASH_Solver*, double*);
    int (*NIiter)(LLJBASH_Solver*, CKTcircuit*, int);
} lljbash_solver;

extern struct LLJBASH_Functions {
    void* dlhandler;

    void (*SetupCscMatrix)(LLJBASH_CscMatrix*, int, long);
    void (*DestroyCscMatrix)(LLJBASH_CscMatrix*);

    void* (*IluSolverCreate)(int);
    void (*IluSolverDestroy)(void*);
    LLJBASH_CscMatrix* (*IluSolverGetMatrix)(void*);
    int (*IluSolverFactorize)(void*, int);

    void* (*GmresCreate)(void);
    void (*GmresDestroy)(void*);
    LLJBASH_GmresParameters* (*GmresGetParameters)(void*);
    void (*GmresSetPreconditioner)(void*, void*);
    int (*GmresSolve)(void*, const LLJBASH_CscMatrix*, const double*, double*, int*);
} lljbash;
