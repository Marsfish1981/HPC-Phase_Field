#ifndef _LP_F_H_
#define _LP_F_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 1)
#define FORT_LINESOLVE linesolve1dsim
#define FORT_ADOTX     adotx1dsim
#define FORT_FLUX      flux1dsim
#endif

#if (BL_SPACEDIM == 2)
#define FORT_GSRB      gsrb2dsim
#define FORT_ADOTX     adotx2dsim
#define FORT_FLUX      flux2dsim
#endif

#if (BL_SPACEDIM == 3)
#define FORT_GSRB      gsrb3dsim
#define FORT_ADOTX     adotx3dsim
#define FORT_FLUX      flux3dsim
#endif

#else

#if (BL_SPACEDIM == 1)
#if defined(BL_FORT_USE_UPPERCASE)
#define FORT_LINESOLVE LINESOLVE1DSIM
#define FORT_ADOTX     ADOTX1DSIM
#define FORT_FLUX      FLUX1DSIM
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_LINESOLVE linesolve1dsim
#define FORT_ADOTX     adotx1dsim
#define FORT_FLUX      flux1dsim
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_LINESOLVE linesolve1dsim_
#define FORT_ADOTX     adotx1dsim_
#define FORT_FLUX      flux1dsim_
#endif
#endif

#if (BL_SPACEDIM == 2)
#if defined(BL_FORT_USE_UPPERCASE)
#define FORT_GSRB      GSRB2DSIM
#define FORT_ADOTX     ADOTX2DSIM
#define FORT_FLUX      FLUX2DSIM
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_GSRB      gsrb2dsim
#define FORT_ADOTX     adotx2dsim
#define FORT_FLUX      flux2dsim
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_GSRB      gsrb2dsim_
#define FORT_ADOTX     adotx2dsim_
#define FORT_FLUX      flux2dsim_
#endif
#endif

#if (BL_SPACEDIM == 3)
#if   defined(BL_FORT_USE_UPPERCASE)
#define FORT_GSRB      GSRB3DSIM
#define FORT_ADOTX     ADOTX3DSIM
#define FORT_FLUX      FLUX3DSIM
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_GSRB      gsrb3dsim
#define FORT_ADOTX     adotx3dsim
#define FORT_FLUX      flux3dsim
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_GSRB      gsrb3dsim_
#define FORT_ADOTX     adotx3dsim_
#define FORT_FLUX      flux3dsim_
#endif

#endif

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 1)
    void FORT_LINESOLVE (
        Real* phi       , ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        const Real* rhs , ARLIM_P(rhs_lo),  ARLIM_P(rhs_hi),
        const Real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const Real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const int* lo, const int* hi, const int *nc,
        const Real *h
        );
    
    void FORT_FLUX(
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const int *xlo, const int *xhi, const int *nc,
        const Real *h,
        Real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi)
        );
#endif
#if (BL_SPACEDIM == 2)
    void FORT_GSRB (
        Real* phi       , ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        const Real* rhs , ARLIM_P(rhs_lo),  ARLIM_P(rhs_hi),
        const Real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const Real* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
        const Real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const Real* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const Real *h, const  int* redblack
        );
    
    void FORT_FLUX(
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
	const int *nc,
        const Real *h,
        Real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi),
        Real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi)
        );
#endif
#if (BL_SPACEDIM == 3)
    void FORT_GSRB (
        Real* phi       , ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        const Real* rhs , ARLIM_P(rhs_lo),  ARLIM_P(rhs_hi),
        const Real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const Real* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
        const Real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const Real* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
        const Real* den4, ARLIM_P(den4),    ARLIM_P(den4_hi),
        const int* m4   , ARLIM_P(m4),      ARLIM_P(m4_hi),
        const Real* den5, ARLIM_P(den5),    ARLIM_P(den5_hi),
        const int* m5   , ARLIM_P(m5),      ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const Real *h, const  int* redblack
        );
    
    void FORT_FLUX(
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
        const int *zlo, const int *zhi, 
	const int *nc,
        const Real *h,
        Real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi),
        Real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi),
        Real* zflux, ARLIM_P(zflux_lo), ARLIM_P(zflux_hi)
        );
#endif

    void FORT_ADOTX(
        Real *y,       ARLIM_P(y_lo), ARLIM_P(y_hi),
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );
}
#endif

#endif /*_LP_F_H_*/
