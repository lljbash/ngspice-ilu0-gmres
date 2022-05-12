#include "ngspice/lljbash.h"
#include <dlfcn.h>
#include "ngspice/spmatrix.h"
#include "../sparse/spdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/memory.h"

void LLJBASH_InitSolver(LLJBASH_Solver* solver) {
    solver->ilu = lljbash.IluSolverCreate(8);
    solver->gmres = lljbash.GmresCreate();
    lljbash.GmresSetPreconditioner(solver->gmres, solver->ilu);
}

void LLJBASH_FreeSolver(LLJBASH_Solver* solver) {
    lljbash.IluSolverDestroy(solver->ilu);
    lljbash.GmresDestroy(solver->gmres);
}

void LLJBASH_ImportMatrix(LLJBASH_Solver* solver, const SMPmatrix* smp) {
    LLJBASH_CscMatrix* mat = lljbash.IluSolverGetMatrix(solver->ilu);
    lljbash.SetupCscMatrix(mat, smp->Size, smp->Elements);
    long nnz = 0;
    for (int col = 0; col < smp->Size; ++col) {
        mat->col_ptr[col] = nnz;
        for (ElementPtr e = smp->FirstInCol[col+1]; e; e = e->NextInCol, ++nnz) {
            mat->row_idx[nnz] = e->Row - 1;
            mat->value[nnz] = e->Real;
        }
    }
    mat->col_ptr[smp->Size] = nnz;
}

void LLJBASH_InitPreconditoner(LLJBASH_Solver* solver) {
    lljbash.IluSolverFactorize(solver->ilu, TRUE);
}

void LLJBASH_GmresSolve(LLJBASH_Solver* solver, double* x) {
    lljbash.GmresSolve(solver->gmres, lljbash.IluSolverGetMatrix(solver->ilu), x, x, NULL);
}

int LLJBASH_NIiter(LLJBASH_Solver* solver, CKTcircuit* ckt, int maxIter) {
    double startTime, *OldCKTstate0 = NULL;
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

            lljbash_solver.ImportMatrix(solver, ckt->CKTmatrix);
            lljbash_solver.InitPreconditoner(solver);

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            if (!OldCKTstate0)
                OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

            startTime = SPfrontEnd->IFseconds();
            /*SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);*/
            lljbash_solver.GmresSolve(solver, ckt->CKTrhs);
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
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

struct LLJBASH_SolverFunctions lljbash_solver;
struct LLJBASH_Functions lljbash;

void __attribute__((constructor)) LLJBASH_LoadFunctions() {
    lljbash.dlhandler = dlopen("/home/lilj/par_ilu0_gmres/lib/libpar-ilu0-gmres_c.so", RTLD_NOW);
#define LLJBASH_LOAD_SYMBOL(SYMBOL) \
    lljbash.SYMBOL = dlsym(lljbash.dlhandler, "LLJBASH_" #SYMBOL)
    LLJBASH_LOAD_SYMBOL(SetupCscMatrix);
    LLJBASH_LOAD_SYMBOL(DestroyCscMatrix);
    LLJBASH_LOAD_SYMBOL(IluSolverCreate);
    LLJBASH_LOAD_SYMBOL(IluSolverDestroy);
    LLJBASH_LOAD_SYMBOL(IluSolverGetMatrix);
    LLJBASH_LOAD_SYMBOL(IluSolverFactorize);
    LLJBASH_LOAD_SYMBOL(GmresCreate);
    LLJBASH_LOAD_SYMBOL(GmresDestroy);
    LLJBASH_LOAD_SYMBOL(GmresGetParameters);
    LLJBASH_LOAD_SYMBOL(GmresSetPreconditioner);
    LLJBASH_LOAD_SYMBOL(GmresSolve);
#undef LLJBASH_LOAD_SYMBOL
    lljbash_solver.Init = LLJBASH_InitSolver;
    lljbash_solver.Free = LLJBASH_FreeSolver;
    lljbash_solver.ImportMatrix = LLJBASH_ImportMatrix;
    lljbash_solver.InitPreconditoner = LLJBASH_InitPreconditoner;
    lljbash_solver.GmresSolve = LLJBASH_GmresSolve;
    lljbash_solver.NIiter = LLJBASH_NIiter;
}

void __attribute__((destructor)) LLJBASH_UnloadFunctions() {
    dlclose(lljbash.dlhandler);
}
