#include "ngspice/lljbash.h"
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <dlfcn.h>
#include "ngspice/spmatrix.h"
#include "../sparse/spdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/memory.h"
#include "ngspice/cpextern.h"
#ifdef XSPICE
#include "ngspice/enh.h"
#endif

void LLJBASH_InitSolver(LLJBASH_Solver* solver) {
    solver->ilu = lljbash.IluSolverCreate(0);
    solver->gmres = lljbash.GmresCreate();
    lljbash.GmresGetParameters(solver->gmres)->tolerance = 1e-4;
    lljbash.GmresGetParameters(solver->gmres)->krylov_subspace_dimension = 60;
    lljbash.GmresGetParameters(solver->gmres)->max_iterations = 2500;
    lljbash.GmresSetPreconditioner(solver->gmres, solver->ilu);
    solver->smp = NULL;
    solver->csr = *lljbash.CSR_MATRIX_DEFAULT;
    solver->need_setup = TRUE;
    solver->element_mapping = NULL;
    solver->first_ilu = TRUE;
    solver->intermediate = NULL;
    memset(&solver->stat, 0, sizeof(solver->stat));
}

void LLJBASH_FreeSolver(LLJBASH_Solver* solver) {
    lljbash.IluSolverDestroy(solver->ilu);
    lljbash.GmresDestroy(solver->gmres);
    FREE(solver->element_mapping);
    FREE(solver->intermediate);
}

void LLJBASH_SetupMatrix(LLJBASH_Solver* solver, SMPmatrix* smp) {
    if (!smp->RowsLinked) {
        SMPpreOrder(smp);
        spcLinkRows(smp);
    }
    solver->smp = smp;

    LLJBASH_CsrMatrix* mat = &solver->csr;
    lljbash.SetupCsrMatrix(mat, smp->Size, smp->Elements + smp->Size);
    solver->element_mapping = TMALLOC(ElementPtr, mat->max_nnz);
    int nnz = 0;
    for (int row = 0; row < smp->Size; ++row) {
        mat->row_ptr[row] = nnz;
        bool has_diag = FALSE;
        for (ElementPtr e = smp->FirstInRow[row+1]; e; e = e->NextInRow) {
            if (!e->Fillin) {
                int col = e->Col - 1;
                if (!has_diag && col >= row) {
                    if (col > row) {
                        mat->col_idx[nnz] = row;
                        mat->value[nnz] = 0.0;
                        solver->element_mapping[nnz] = NULL;
                        ++nnz;
                        printf("no diag %d\n", row);
                    }
                    else {
                        if (e->Real == 0.0) {
                            printf("diag %d is 0\n", row);
                        }
                    }
                    has_diag = TRUE;
                }
                mat->col_idx[nnz] = col;
                mat->value[nnz] = e->Real;
                solver->element_mapping[nnz] = e;
                ++nnz;
            }
        }
        if (!has_diag) {
            mat->col_idx[nnz] = row;
            mat->value[nnz] = 0.0;
            solver->element_mapping[nnz] = NULL;
            ++nnz;
            printf("no diag %d\n", row);
        }
    }
    mat->row_ptr[smp->Size] = nnz;

    LLJBASH_CsrMatrix* ilumat = lljbash.IluSolverGetMatrix(solver->ilu);
    lljbash.CopyCsrMatrix(ilumat, mat);

    if (!solver->intermediate) {
        solver->intermediate = TMALLOC(double, mat->size * 2);
    }
    double* sol = solver->intermediate + mat->size;
    /*srand(114514);*/
    for (int i = 0; i < mat->size; ++i) {
        /*sol[i] = ((double) rand() / ((double) RAND_MAX + 1.0));*/
        sol[i] = 0.0;
    }

    solver->need_setup = FALSE;

    printf("n = %d, nnz = %d\n", mat->size, nnz);
}

void LLJBASH_ImportMatrix(LLJBASH_Solver* solver, SMPmatrix* smp) {
    if (solver->need_setup) {
        LLJBASH_SetupMatrix(solver, smp);
    }
    else {
        LLJBASH_CsrMatrix* mat = &solver->csr;
        for (int i = 0; i < mat->row_ptr[mat->size]; ++i) {
            mat->value[i] = solver->element_mapping[i] ? solver->element_mapping[i]->Real : 0.0;
        }
    }
}

void LLJBASH_InitPreconditoner(LLJBASH_Solver* solver) {
    LLJBASH_CsrMatrix* ilumat = lljbash.IluSolverGetMatrix(solver->ilu);
    memcpy(ilumat->value, solver->csr.value, sizeof(double[ilumat->row_ptr[ilumat->size]]));
    if (solver->first_ilu) {
        LLJBASH_Tic();
        lljbash.IluSolverSetup(solver->ilu);
        solver->first_ilu = FALSE;
        solver->stat.ilu_setup_time += LLJBASH_Toc();
    }
    LLJBASH_Tic();
    lljbash.IluSolverFactorize(solver->ilu);
    solver->stat.total_precon_time += LLJBASH_Toc();
}

int LLJBASH_GmresSolve(LLJBASH_Solver* solver, double* x) {
    LLJBASH_CsrMatrix* mat = &solver->csr;
    int* p_ext_order = &solver->smp->IntToExtRowMap[mat->size];
    double* rhs = solver->intermediate;
    double* sol = solver->intermediate + mat->size;
    for (int i = mat->size; i > 0; i--) {
        rhs[i-1] = x[*(p_ext_order--)-1];
    }
    /*memcpy(sol, rhs, sizeof(double[mat->size]));*/
    int its;
    int success = lljbash.GmresSolve(solver->gmres, mat, rhs, sol, &its);
    if (!success) {
        /*for (int i = 0; i < mat->size; ++i) {*/
            /*sol[i] = ((double) rand() / ((double) RAND_MAX + 1.0));*/
        /*}*/
        /*success = lljbash.GmresSolve(solver->gmres, mat, rhs, sol, NULL);*/
        /*if (!success) {*/
            fputs("\nGMRES failed!\n", stderr);
            return E_INTERN;
        /*}*/
    }
    printf("its: %d\n", its);
    p_ext_order = &solver->smp->IntToExtColMap[mat->size];
    for (int i = mat->size; i > 0; i--) {
        x[*(p_ext_order--)-1] = sol[i-1];
    }
    return OK;
}

int LLJBASH_NIiter(LLJBASH_Solver* solver, CKTcircuit* ckt, int maxIter) {
    double *OldCKTstate0 = NULL;
    int error, i;

    int iterno = 0;
    int ipass = 0;

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;

    if ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)) {
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        error = CKTload(ckt);
        if (error)
            return(error);
        return(OK);
    }

#ifdef WANT_SENSE2
    if (ckt->CKTsenInfo) {
        error = NIsenReinit(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt);
        if (error) {
#ifdef STEPDEBUG
            printf("re-init returned error \n");
#endif
            return(error);
        }
    }

    /* OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1); */

    for (;;) {

        ckt->CKTnoncon = 0;

#ifdef NEWPRED
        if (!(ckt->CKTmode & MODEINITPRED))
#endif
        {

            error = CKTload(ckt);
            /* printf("loaded, noncon is %d\n", ckt->CKTnoncon); */
            /* fflush(stdout); */
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                printf("load returned error \n");
#endif
                FREE(OldCKTstate0);
                return (error);
            }

            puts("import");
            lljbash_solver.ImportMatrix(solver, ckt->CKTmatrix);
            puts("ILU0");
            lljbash_solver.InitPreconditoner(solver);

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

            /*SMPprint(ckt->CKTmatrix, NULL);*/
            /*SMPprintRHS(ckt->CKTmatrix, NULL, ckt->CKTrhs, ckt->CKTirhs);*/
            puts("GMRES");
            LLJBASH_Tic();
            error = lljbash_solver.GmresSolve(solver, &ckt->CKTrhs[1]);
            solver->stat.total_gmres_time += LLJBASH_Toc();
            if (error) {
                exit(-1);
            }
            /*SMPprintRHS(ckt->CKTmatrix, NULL, ckt->CKTrhs, ckt->CKTirhs);*/
            /*static char buf[256];*/
            /*fgets(buf, 256, stdin);*/
#ifdef STEPDEBUG
            /*XXXX*/
            if (ckt->CKTrhs[0] != 0.0)
                printf("NIiter: CKTrhs[0] = %g\n", ckt->CKTrhs[0]);
            if (ckt->CKTrhsSpare[0] != 0.0)
                printf("NIiter: CKTrhsSpare[0] = %g\n", ckt->CKTrhsSpare[0]);
            if (ckt->CKTrhsOld[0] != 0.0)
                printf("NIiter: CKTrhsOld[0] = %g\n", ckt->CKTrhsOld[0]);
            /*XXXX*/
#endif
            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");
#ifdef STEPDEBUG
                    fprintf(stderr, "too many iterations without convergence: %d iter's (max iter == %d)\n",
                    iterno, maxIter);
#endif
                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }

            if ((ckt->CKTnoncon == 0) && (iterno != 1))
                ckt->CKTnoncon = NIconvTest(ckt);
            else
                ckt->CKTnoncon = 1;

#ifdef STEPDEBUG
            printf("noncon is %d\n", ckt->CKTnoncon);
#endif
        }

        if ((ckt->CKTnodeDamping != 0) && (ckt->CKTnoncon != 0) &&
            ((ckt->CKTmode & MODETRANOP) || (ckt->CKTmode & MODEDCOP)) &&
            (iterno > 1))
        {
            CKTnode *node;
            double diff, maxdiff = 0;
            for (node = ckt->CKTnodes->next; node; node = node->next)
                if (node->type == SP_VOLTAGE) {
                    diff = fabs(ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number]);
                    if (maxdiff < diff)
                        maxdiff = diff;
                }

            if (maxdiff > 10) {
                double damp_factor = 10 / maxdiff;
                if (damp_factor < 0.1)
                    damp_factor = 0.1;
                for (node = ckt->CKTnodes->next; node; node = node->next) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    ckt->CKTrhs[node->number] =
                        ckt->CKTrhsOld[node->number] + (damp_factor * diff);
                }
                for (i = 0; i < ckt->CKTnumStates; i++) {
                    diff = ckt->CKTstate0[i] - OldCKTstate0[i];
                    ckt->CKTstate0[i] = OldCKTstate0[i] + (damp_factor * diff);
                }
            }
        }

        if (ckt->CKTmode & MODEINITFLOAT) {
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
                if (ipass)
                    ckt->CKTnoncon = ipass;
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return(OK);
            }
        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
            ipass = 1;
        } else if (ckt->CKTmode & MODEINITSMSIG) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITTRAN) {
            if (iterno <= 1)
                ckt->CKTniState |= NISHOULDREORDER;
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITPRED) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else {
            ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
            printf("bad initf state \n");
#endif
            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        /* build up the lvnim1 array from the lvn array */
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        /* printf("after loading, after solving\n"); */
        /* CKTdump(ckt); */
    }
    /*NOTREACHED*/

}

int dynamic_gmin(CKTcircuit *, long int, long int, int);
int spice3_gmin(CKTcircuit *, long int, long int, int);
int new_gmin(CKTcircuit*, long int, long int, int);
int gillespie_src(CKTcircuit *, long int, long int, int);
int spice3_src(CKTcircuit *, long int, long int, int);

int LLJBASH_CKTop(LLJBASH_Solver* solver, CKTcircuit *ckt, long firstmode, long continuemode, int iterlim) {
    int converged;

#ifdef HAS_PROGREP
    SetAnalyse("op", 0);
#endif

    ckt->CKTmode = firstmode;

    if (!ckt->CKTnoOpIter) {
#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call =
            (ckt->CKTnumGminSteps <= 0) && (ckt->CKTnumSrcSteps <= 0);
#endif
        converged = LLJBASH_NIiter (solver, ckt, iterlim);
        if (converged == 0)
            return converged;   /* successfull */
    } else {
        converged = 1;          /* the 'go directly to gmin stepping' option */
    }


    /* no convergence on the first try, so we do something else */
    /* first, check if we should try gmin stepping */

    if (ckt->CKTnumGminSteps >= 1) {
        if (ckt->CKTnumGminSteps == 1) {
            /* only the old gmin */
            if (cp_getvar("dyngmin", CP_BOOL, NULL, 0)) {
                converged = dynamic_gmin(ckt, firstmode, continuemode, iterlim);
            }
            /* first the old gmin, then the new gmin */
            else {
                converged = dynamic_gmin(ckt, firstmode, continuemode, iterlim);
                if(converged != 0) {
                    converged = new_gmin(ckt, firstmode, continuemode, iterlim);
                }

            }
        }
        else {
            converged = spice3_gmin(ckt, firstmode, continuemode, iterlim);
        }
        if (converged == 0) /* If gmin-stepping worked... move out */
            return converged;
    }

    /* ... otherwise try stepping sources ...
     * now, we'll try source stepping - we scale the sources
     * to 0, converge, then start stepping them up until they
     * are at their normal values
     */

    if (ckt->CKTnumSrcSteps >= 1) {
        if (ckt->CKTnumSrcSteps == 1)
            converged = gillespie_src(ckt, firstmode, continuemode, iterlim);
        else
            converged = spice3_src(ckt, firstmode, continuemode, iterlim);
        if (converged == 0) /* If gmin-stepping worked... move out */
            return converged;
    }

    /* If command 'optran' is not given, the function
       returns immediately with the previous 'converged' */
    converged = OPtran(ckt, converged);

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif

    return converged;
}

LLJBASH_Solver lljbash_global_solver;
LLJBASH_Solver* lljbash_this = &lljbash_global_solver;
struct LLJBASH_SolverFunctions lljbash_solver;
struct LLJBASH_Functions lljbash;

uint64_t LLJBASH_cycles;

void __attribute__((constructor)) LLJBASH_LoadFunctions() {
    lljbash.dlhandler = dlopen("/home/lilj/par_ilu0_gmres/lib/libpar-ilu0-gmres_c.so", RTLD_NOW);
#define LLJBASH_LOAD_SYMBOL(SYMBOL) \
    lljbash.SYMBOL = dlsym(lljbash.dlhandler, "LLJBASH_" #SYMBOL)
    LLJBASH_LOAD_SYMBOL(SetupCsrMatrix);
    LLJBASH_LOAD_SYMBOL(DestroyCsrMatrix);
    LLJBASH_LOAD_SYMBOL(CopyCsrMatrix);
    LLJBASH_LOAD_SYMBOL(CSR_MATRIX_DEFAULT);
    LLJBASH_LOAD_SYMBOL(IluSolverCreate);
    LLJBASH_LOAD_SYMBOL(IluSolverDestroy);
    LLJBASH_LOAD_SYMBOL(IluSolverGetMatrix);
    LLJBASH_LOAD_SYMBOL(IluSolverSetup);
    LLJBASH_LOAD_SYMBOL(IluSolverFactorize);
    LLJBASH_LOAD_SYMBOL(GmresCreate);
    LLJBASH_LOAD_SYMBOL(GmresDestroy);
    LLJBASH_LOAD_SYMBOL(GmresGetParameters);
    LLJBASH_LOAD_SYMBOL(GmresGetStat);
    LLJBASH_LOAD_SYMBOL(GmresSetPreconditioner);
    LLJBASH_LOAD_SYMBOL(GmresSolve);
#undef LLJBASH_LOAD_SYMBOL
    lljbash_solver.Init = LLJBASH_InitSolver;
    lljbash_solver.Free = LLJBASH_FreeSolver;
    lljbash_solver.ImportMatrix = LLJBASH_ImportMatrix;
    lljbash_solver.InitPreconditoner = LLJBASH_InitPreconditoner;
    lljbash_solver.GmresSolve = LLJBASH_GmresSolve;
    lljbash_solver.NIiter = LLJBASH_NIiter;
    lljbash_solver.CKTop = LLJBASH_CKTop;
}

void __attribute__((destructor)) LLJBASH_UnloadFunctions() {
    dlclose(lljbash.dlhandler);
}
