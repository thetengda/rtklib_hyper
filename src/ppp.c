/*------------------------------------------------------------------------------
 * ppp.c : precise point positioning
 *
 *          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
 *
 * options : -DIERS_MODEL  use IERS tide model
 *           -DOUTSTAT_AMB output ambiguity parameters to solution status
 *
 * references :
 *    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
 *    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
 *        2003, November 2003
 *    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
 *        Space Technology Library, 2004
 *    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
 *        May 2009
 *    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
 *        Code Biases, URA
 *    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
 *        celestial reference frames, Geophys. Res. Let., 1997
 *    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
 *         Conventions (2010), 2010
 *    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
 *        GPS Solutions, 13:1-12, 2009
 *    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
 *        modeling, InsideGNSS, September, 2010
 *    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
 *        Space Research, 2010
 *    [11] IGS MGEX (http://igs.org/mgex)
 *
 * version : $Revision:$ $Date:$
 * history : 2010/07/20 1.0  new
 *                           added api:
 *                               tidedisp()
 *           2010/12/11 1.1  enable exclusion of eclipsing satellite
 *           2012/02/01 1.2  add gps-glonass h/w bias correction
 *                           move windupcorr() to rtkcmn.c
 *           2013/03/11 1.3  add otl and pole tides corrections
 *                           involve iers model with -DIERS_MODEL
 *                           change initial variances
 *                           suppress acos domain error
 *           2013/09/01 1.4  pole tide model by iers 2010
 *                           add mode of ionosphere model off
 *           2014/05/23 1.5  add output of trop gradient in solution status
 *           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
 *                           fix bug on m2 computation in tide_pole()
 *           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
 *           2015/05/10 1.8  add function to detect slip by MW-LC jump
 *                           fix ppp solutin problem with large clock variance
 *           2015/06/08 1.9  add precise satellite yaw-models
 *                           cope with day-boundary problem of satellite clock
 *           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
 *                           pppoutsolsat() -> pppoutstat()
 *           2015/11/13 1.11 add L5-receiver-dcb estimation
 *                           merge post-residual validation by rnx2rtkp_test
 *                           support support option opt->pppopt=-GAP_RESION=nnnn
 *           2016/01/22 1.12 delete support for yaw-model bug
 *                           add support for ura of ephemeris
 *           2018/10/10 1.13 support api change of satexclude()
 *           2020/11/30 1.14 use sat2freq() to get carrier frequency
 *                           use E1-E5b for Galileo iono-free LC
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

// #define RCVCLK_BDS23 1
// #include <assert.h>
// #include <memory.h>
#define MAX_STD_FIX 0.15 /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4   /* min satellite number for solution */
#define THRES_REJECT 4.0 /* reject threshold of posfit-res (sigma) */

// #define VAR_POS SQR     (60.0)    /* init variance receiver position (m^2) */
// #define VAR_VEL SQR     (10.0)    /* init variance of receiver vel ((m/s)^2) */
// #define VAR_ACC SQR     (10.0)    /* init variance of receiver acc ((m/ss)^2) */
// #define VAR_CLK SQR     (60.0)    /* init variance receiver clock (m^2) */
// #define VAR_ZTD SQR     (0.6 )     /* init variance ztd (m^2) */
// #define VAR_GRA SQR     (0.01)    /* init variance gradient (m^2) */
// #define VAR_DCB SQR     (30.0)    /* init variance dcb (m^2) */
// #define VAR_BIAS SQR    (60.0)   /* init variance phase-bias (m^2) */
// #define VAR_IONO SQR    (60.0)   /* init variance iono-delay */
// #define VAR_GLO_IFB SQR (0.6 ) /* variance of glonass ifb */

#define ERR_SAAS 0.3   /* saastamoinen model error std (m) */
#define ERR_BRDCI 0.5  /* broadcast iono model error factor */
#define ERR_CBIAS 0.3  /* code bias error std (m) */
#define REL_HUMI 0.7   /* relative humidity for saastamoinen model */
#define GAP_RESION 120 /* default gap to reset ionos parameters (ep) */

// #define EFACT_GPS_L5 10.0 /* error factor of GPS/QZS L5 */

/* number and index of states */
#define NF(opt) ((opt)->ionoopt == IONOOPT_IFLC ? 1 : (opt)->nf)
#define NP(opt) ((opt)->dynamics ? 9 : 3)
#define NC(opt) (NSYS + 1)
#define NT(opt) ((opt)->tropopt < TROPOPT_EST ? 0 : ((opt)->tropopt == TROPOPT_EST ? 1 : 3))
#define NI(opt) ((opt)->ionoopt == IONOOPT_EST ? MAXSAT : 0)
#define ND(opt) ((opt)->nf >= 3 ? 1 : 0)
#define NR(opt) (NP(opt) + NC(opt) + NT(opt) + NI(opt) + ND(opt))
#define NB(opt) (NF(opt) * MAXSAT)
#define NX(opt) (NR(opt) + NB(opt))
#define IC(s, opt) (NP(opt) + (s))
#define IT(opt) (NP(opt) + NC(opt))
#define II(s, opt) (NP(opt) + NC(opt) + NT(opt) + (s) -1)
#define ID(opt) (NP(opt) + NC(opt) + NT(opt) + NI(opt))
#define IB(s, f, opt) (NR(opt) + MAXSAT * (f) + (s) -1)

#define MAXNVOBS (MAXSAT + 1) /* NSAT vion + 1 vtrp */
#define ENABLE_VIONO 0
#define ENABLE_VTROP 0

//pppar 
// #define MIN_ARC_GAP     300.0       /* min arc gap (s) */
// #define LOG_PI          1.14472988584940017 /* log(pi) */
// #define SQRT2           1.41421356237309510 /* sqrt(2) */

// #define CONST_AMB       0.001       /* constraint to fixed ambiguity */
// #define THRES_RES       0.3         /* threashold of residuals test (m) */

// #define SWAP_I(x,y)     do {int _z=x; x=y; y=_z;} while (0)
// #define SWAP_D(x,y)     do {double _z=x; x=y; y=_z;} while (0)

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i) {
    if (rtk->sol.stat == SOLQ_FIX)
        return SQRT(rtk->Pa[i + i * rtk->nx]);
    return SQRT(rtk->P[i + i * rtk->nx]);
}
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t *rtk, char *buff) {
    ssat_t *ssat;
    double tow, pos[3], vel[3], acc[3], *x;
    int i, j, k, week;
    char id[32], *p = buff;

    if (!rtk->sol.stat)
        return 0;

    trace(3, "pppoutstat:\n");

    tow = time2gpst(rtk->sol.time, &week);

    x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;

    /* receiver position */
    p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow, rtk->sol.stat, x[0], x[1], x[2],
                 STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->x + 3, vel);
        ecef2enu(pos, rtk->x + 6, acc);
        p += sprintf(p,
                     "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
                     "%.4f,%.5f,%.5f,%.5f\n",
                     week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0);
    }
    /* receiver clocks */
    i = IC(0, &rtk->opt);
    p += sprintf(p, "$CLK,%d,%.3f,%d,%d", week, tow, rtk->sol.stat, NSYS+1);
    for(j=0;j<NSYS+1;++j){
        p += sprintf(p,",%.3f",x[i+j] * 1E9 / CLIGHT); // clocks(ns) 
    }
    for(j=0;j<NSYS+1;++j){
        p += sprintf(p,",%.3f", STD(rtk, i+j) * 1E9 / CLIGHT); // std(ns)
    }
    p += sprintf(p, "\n");

    /* tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat, 1, x[i], STD(rtk, i));
    }
    if (rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow, rtk->sol.stat, 1, x[i + 1], x[i + 2],
                     STD(rtk, i + 1), STD(rtk, i + 2));
    }
    /* ionosphere parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST) {
        for (i = 0; i < MAXSAT; i++) {
            ssat = rtk->ssat + i;
            if (!ssat->vs)
                continue;
            j = II(i + 1, &rtk->opt);
            if (rtk->x[j] == 0.0)
                continue;
            satno2id(i + 1, id);
            p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow, rtk->sol.stat, id,
                         rtk->ssat[i].azel[0] * R2D, rtk->ssat[i].azel[1] * R2D, x[j], STD(rtk, j));
        }
    }
#ifdef OUTSTAT_AMB
    /* ambiguity parameters */
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < NF(&rtk->opt); j++) {
            k = IB(i + 1, j, &rtk->opt);
            if (rtk->x[k] == 0.0)
                continue;
            satno2id(i + 1, id);
            p +=
                sprintf(p, "$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat, id, j + 1, x[k], STD(rtk, k));
        }
#endif
    return (int) (p - buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs) {
    double rsun[3], esun[3], r, ang, erpv[5] = {0}, cosa;
    int i, j;
    const char *type;

    trace(3, "testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
    normv3(rsun, esun);

    for (i = 0; i < n; i++) {
        type = nav->pcvs[obs[i].sat - 1].type;

        if ((r = norm(rs + i * 6, 3)) <= 0.0)
            continue;

        /* only block IIA */
        if (*type && !strstr(type, "BLOCK IIA"))
            continue;

        /* sun-earth-satellite angle */
        cosa = dot(rs + i * 6, esun, 3) / r;
        cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
        ang  = acos(cosa);

        /* test eclipse */
        if (ang < PI / 2.0 || r * sin(ang) > RE_WGS84)
            continue;

        trace(3, "eclipsing sat excluded %s sat=%2d\n", time_str(obs[0].time, 0), obs[i].sat);

        for (j = 0; j < 3; j++)
            rs[j + i * 6] = 0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu) {
    if (fabs(beta) < 1E-12 && fabs(mu) < 1E-12)
        return PI;
    return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu, double *yaw) {
    *yaw = yaw_nominal(beta, mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt, const double *rs, double *exs, double *eys) {
    double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
    double yaw, cosy, siny, erpv[5] = {0};
    int i;

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

    /* beta and orbit angle */
    matcpy(ri, rs, 6, 1);
    ri[3] -= OMGE * ri[1];
    ri[4] += OMGE * ri[0];
    cross3(ri, ri + 3, n);
    cross3(rsun, n, p);
    if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) || !normv3(p, ep))
        return 0;
    beta = PI / 2.0 - acos(dot(esun, en, 3));
    E    = acos(dot(es, ep, 3));
    mu   = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
    if (mu < -PI / 2.0)
        mu += 2.0 * PI;
    else if (mu >= PI / 2.0)
        mu -= 2.0 * PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat, type, opt, beta, mu, &yaw))
        return 0;

    /* satellite fixed x,y-vector */
    cross3(en, es, ex);
    cosy = cos(yaw);
    siny = sin(yaw);
    for (i = 0; i < 3; i++) {
        exs[i] = -siny * en[i] + cosy * ex[i];
        eys[i] = -cosy * en[i] - siny * ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt, const double *rs, const double *rr,
                     double *phw) {
    double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
    int i;

    if (opt <= 0)
        return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time, sat, type, opt, rs, exs, eys))
        return 0;

    /* unit vector satellite to receiver */
    for (i = 0; i < 3; i++)
        r[i] = rr[i] - rs[i];
    if (!normv3(r, ek))
        return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);
    xyz2enu(pos, E);
    exr[0] = E[1];
    exr[1] = E[4];
    exr[2] = E[7]; /* x = north */
    eyr[0] = -E[0];
    eyr[1] = -E[3];
    eyr[2] = -E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek, eys, eks);
    cross3(ek, eyr, ekr);
    for (i = 0; i < 3; i++) {
        ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];
        dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];
    }
    cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
    if (cosp < -1.0)
        cosp = -1.0;
    else if (cosp > 1.0)
        cosp = 1.0;
    ph = acos(cosp) / 2.0 / PI;
    cross3(ds, dr, drs);
    if (dot(ek, drs, 3) < 0.0)
        ph = -ph;

    *phw = ph + floor(*phw - ph + 0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int idx, int type, const prcopt_t *opt) {
    double fact = 1.0, sinel = sin(el);

    if (type == 1)
        fact *= opt->eratio[idx == 0 ? 0 : 1]; // todo adjuct varance by frequency index 
    if(sys == SYS_GPS) fact*=EFACT_GPS;
    else if(sys == SYS_GLO) fact*=EFACT_GLO;
    else if(sys == SYS_GAL) fact*=EFACT_GAL;
    else if(sys == SYS_QZS) fact*=EFACT_QZS;
    else if(sys == SYS_SBS) fact*=EFACT_SBS;
    else if(sys == SYS_CMP) fact*=EFACT_CMP;
    else if(sys == SYS_IRN) fact*=EFACT_IRN;
    else if(sys == SYS_LEO) fact*=EFACT_LEO;
    else fact*=EFACT_GPS;

    // if (sys == SYS_GPS || sys == SYS_QZS) {
    //     if (idx == 2)
    //         fact *= EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    // }
    if (opt->ionoopt == IONOOPT_IFLC)
        fact *= 3.0;
    return SQR(fact * opt->err[1]) + SQR(fact * opt->err[2] / sinel);
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i) {
    int j;
    rtk->x[i] = xi;
    for (j = 0; j < rtk->nx; j++) {
        rtk->P[i + j * rtk->nx] = rtk->P[j + i * rtk->nx] = i == j ? var : 0.0;
    }
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav) {
    double freq1, freq2;

    freq1 = sat2freq(obs->sat, obs->code[0], nav);
    freq2 = sat2freq(obs->sat, obs->code[1], nav);
    if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0)
        return 0.0;
    return (obs->L[0] / freq1 - obs->L[1] / freq2) * CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav) {
    double freq1, freq2;

    freq1 = sat2freq(obs->sat, obs->code[0], nav);
    freq2 = sat2freq(obs->sat, obs->code[1], nav);

    if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0 || obs->P[0] == 0.0 || obs->P[1] == 0.0)
        return 0.0;
    return (obs->L[0] - obs->L[1]) * CLIGHT / (freq1 - freq2) -
           (freq1 * obs->P[0] + freq2 * obs->P[1]) / (freq1 + freq2); // m, (N2-N1)*CLIGHT/(freq1-fre12)
}

static int model_bds2_cbias(int sat, const double *azel, double *mp) {

    static const double IGSOCOEF[3][10] = {/* m */    
        {-0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35}, // B1
        {-0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32}, // B3
        {-0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33}, // B2
    };
    static const double MEOCOEF[3][10] = {/* m */
        {-0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05}, // B1
        {-0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47}, // B3
        {-0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69}, // B2
    };

    int i, ielev, sattype;
    double a;

    if (bds2or3(sat) != 2) {
        for (i = 0; i < 3; ++i)
            mp[i] = 0.0;
        return 1;
    }

    ielev   = azel[1] * R2D / 10;
    sattype = sat2orbtype(sat);
    if (ielev < 0) {
        for (i = 0; i < 3; ++i) {
            if (sattype == ORB_MEO) {
                mp[i] = MEOCOEF[i][0];
            } else if (sattype == ORB_IGSO) {
                mp[i] = IGSOCOEF[i][0];
            } else {
                mp[i] = 0.0;
            }
        }
    } else if (ielev >= 9) {
        for (i = 0; i < 3; ++i) {
            if (sattype == ORB_MEO) {
                mp[i] = MEOCOEF[i][9];
            } else if (sattype ==ORB_IGSO) {
                mp[i] = IGSOCOEF[i][9];
            } else {
                mp[i] = 0.0;
            }
        }
    } else {
        a = (azel[1] * R2D - ielev * 10) / 10;
        for (i = 0; i < 3; ++i) {
            if (sattype == ORB_IGSO) {
                mp[i] = (1.0 - a) * MEOCOEF[i][ielev] + a * MEOCOEF[i][ielev + 1];
            } else if (sattype == ORB_IGSO) {
                mp[i] = (1.0 - a) * IGSOCOEF[i][ielev] + a * IGSOCOEF[i][ielev + 1];
            } else {
                mp[i] = 0.0;
            }
        }
    }
    return 1;
}

/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel, const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P, double *Lc, double *Pc) {
    double freq[NFREQ] = {0}, C1, C2;
    int i, sys = satsys(obs->sat, NULL);
    double cbias[NFREQ] = {0};

    if (sys == SYS_CMP)
        model_bds2_cbias(obs->sat, azel, cbias);

    for (i = 0; i < NFREQ; i++) {
        L[i] = P[i] = 0.0;
        freq[i]     = sat2freq(obs->sat, obs->code[i], nav);
        if (freq[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0)
            continue;
        if (testsnr(0, 0, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask))
            continue;

        /* antenna phase center and phase windup correction */
        L[i] = obs->L[i] * CLIGHT / freq[i] - dants[i] - dantr[i] - phw * CLIGHT / freq[i];
        P[i] = obs->P[i] - dants[i] - dantr[i] + cbias[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        P[i] -= nav->satosb[obs->sat - 1][0][obs->code[i] - 1];
        L[i] -= nav->satosb[obs->sat - 1][1][obs->code[i] - 1];
    }
    /* iono-free LC */
    *Lc = *Pc = 0.0;
    if (freq[0] == 0.0 || freq[1] == 0.0)
        return;
    C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
    C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

    if (L[0] != 0.0 && L[1] != 0.0)
        *Lc = C1 * L[0] + C2 * L[1];
    if (P[0] != 0.0 && P[1] != 0.0)
        *Pc = C1 * P[0] + C2 * P[1];
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n) {
    int i, j;

    trace(3, "detslp_ll: n=%d\n", n);

    for (i = 0; i < n && i < MAXOBS; i++)
        for (j = 0; j < rtk->opt.nf; j++) {
            if (obs[i].L[j] == 0.0 || !(obs[i].LLI[j] & 3))
                continue;

            trace(3, "detslp_ll: slip detected sat=%2d f=%d\n", obs[i].sat, j + 1);

            rtk->ssat[obs[i].sat - 1].slip[j] = 1;
        }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    double g0, g1;
    int i, j;

    trace(3, "detslp_gf: n=%d\n", n);

    for (i = 0; i < n && i < MAXOBS; i++) {

        if ((g1 = gfmeas(obs + i, nav)) == 0.0)
            continue;

        g0                              = rtk->ssat[obs[i].sat - 1].gf[0];
        rtk->ssat[obs[i].sat - 1].gf[0] = g1;

        trace(4, "detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n", obs[i].sat, g0, g1);

        if (g0 != 0.0 && fabs(g1 - g0) > rtk->opt.thresslip[0]) {
            trace(3, "detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n", obs[i].sat, g0, g1);

            for (j = 0; j < rtk->opt.nf; j++)
                rtk->ssat[obs[i].sat - 1].slip[j] |= 1;
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    double w0, w1;
    int i, j;

    trace(3, "detslp_mw: n=%d\n", n);

    for (i = 0; i < n && i < MAXOBS; i++) {
        if ((w1 = mwmeas(obs + i, nav)) == 0.0)
            continue;

        w0                              = rtk->ssat[obs[i].sat - 1].mw[0];
        rtk->ssat[obs[i].sat - 1].mw[0] = w1;

        trace(4, "detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f\n", obs[i].sat, w0, w1);

        if (w0 != 0.0 && fabs(w1 - w0) > rtk->opt.thresslip[1]) {
            trace(3, "detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n", obs[i].sat, w0, w1);

            for (j = 0; j < rtk->opt.nf; j++)
                rtk->ssat[obs[i].sat - 1].slip[j] |= 1;
        }
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk) {
    double *F, *P, *FP, *x, *xp, pos[3], Q[9] = {0}, Qv[9];
    int i, j, *ix, nx;

    trace(3, "udpos_ppp:\n");

    /* fixed mode */
    if (rtk->opt.mode == PMODE_PPP_FIXED) {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->opt.ru[i], 1E-8, i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x, 3) <= 0.0) {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.std[0]), i);
        if (rtk->opt.dynamics) {
            for (i = 3; i < 6; i++)
                initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.std[1]), i);
            for (i = 6; i < 9; i++)
                initx(rtk, 1E-6, SQR(rtk->opt.std[2]), i);
        }
    }
    /* static ppp mode */
    if (rtk->opt.mode == PMODE_PPP_STATIC) {
        for (i = 0; i < 3; i++) {
            rtk->P[i * (1 + rtk->nx)] += SQR(rtk->opt.prn[0]) * fabs(rtk->tt);
        }
        return;
    }
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i = 0; i < 3; i++) {
            initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.std[0]), i);
        }
        return;
    }
    /* generate valid state index */
    ix = imat(rtk->nx, 1);
    for (i = nx = 0; i < rtk->nx; i++) {
        if (rtk->x[i] != 0.0 && rtk->P[i + i * rtk->nx] > 0.0)
            ix[nx++] = i;
    }
    if (nx < 9) {
        free(ix);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F  = eye(nx);
    P  = mat(nx, nx);
    FP = mat(nx, nx);
    x  = mat(nx, 1);
    xp = mat(nx, 1);

    for (i = 0; i < 6; i++) {
        F[i + (i + 3) * nx] = rtk->tt;
    }
    for (i = 0; i < 3; i++) {
        F[i + (i + 6) * nx] = SQR(rtk->tt) / 2.0;
    }
    for (i = 0; i < nx; i++) {
        x[i] = rtk->x[ix[i]];
        for (j = 0; j < nx; j++) {
            P[i + j * nx] = rtk->P[ix[i] + ix[j] * rtk->nx];
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN", nx, 1, nx, 1.0, F, x, 0.0, xp);
    matmul("NN", nx, nx, nx, 1.0, F, P, 0.0, FP);
    matmul("NT", nx, nx, nx, 1.0, FP, F, 0.0, P);

    for (i = 0; i < nx; i++) {
        rtk->x[ix[i]] = xp[i];
        for (j = 0; j < nx; j++) {
            rtk->P[ix[i] + ix[j] * rtk->nx] = P[i + j * nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0] = Q[4] = SQR(rtk->opt.prn[1]) * fabs(rtk->tt);
    Q[8]        = SQR(rtk->opt.prn[2]) * fabs(rtk->tt);
    ecef2pos(rtk->x, pos);
    covecef(pos, Q, Qv);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            rtk->P[i + 6 + (j + 6) * rtk->nx] += Qv[i + j * 3];
        }
    free(ix);
    free(F);
    free(P);
    free(FP);
    free(x);
    free(xp);
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk) {
    double dtr;
    int i,j;

    trace(3, "udclk_ppp:\n");

    /* initialize every epoch for clock (white noise) */
    for (i = 0; i < NSYS + 1; i++) {
        // if (rtk->opt.sateph == EPHOPT_PREC) {
            /* time of prec ephemeris is based gpst */
            /* negelect receiver inter-system bias  */
            // dtr = rtk->sol.dtr[0];
            // fixed: dtr are derived from spp for each system
            dtr = rtk->sol.dtr[i];
            // if(dtr==0){ // spp do not get dtr(i)
            //     for(j=0;j<NSYS+1;++j){
            //         if(rtk->sol.dtr[j]!=0)
            //             dtr = rtk->sol.dtr[j]; // use dtr(G) instead
            //     }
            // }
        // } else {
        //     dtr = i == 0 ? rtk->sol.dtr[0] : rtk->sol.dtr[0] + rtk->sol.dtr[i];
        // }
        initx(rtk, CLIGHT * dtr, SQR(rtk->opt.std[3]), IC(i, &rtk->opt));
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk) {
    double pos[3], azel[] = {0.0, PI / 2.0}, ztd, var;
    int i = IT(&rtk->opt), j;

    trace(3, "udtrop_ppp:\n");

    if (rtk->x[i] == 0.0) {
        ecef2pos(rtk->sol.rr, pos);
        tropcorr(rtk->sol.time, NULL, pos, azel, rtk->opt.tropopt, &ztd, &var);
        initx(rtk, ztd, var, i);

        if (rtk->opt.tropopt >= TROPOPT_ESTG) {
            for (j = i + 1; j < i + 3; j++)
                initx(rtk, 1E-6, SQR(rtk->opt.std[5]), j);
        }
    } else {
        rtk->P[i + i * rtk->nx] += SQR(rtk->opt.prn[4]) * fabs(rtk->tt);

        if (rtk->opt.tropopt >= TROPOPT_ESTG) {
            for (j = i + 1; j < i + 3; j++) {
                rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[5]) * fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    double freq1, freq2, ion, sinel, pos[3], *azel;
    char *p;
    int i, j;

    trace(3, "udiono_ppp:\n");

    for (i = 0; i < MAXSAT; i++) {
        j = II(i + 1, &rtk->opt);
        if (rtk->x[j] != 0.0 && (int) rtk->ssat[i].outc[0] > rtk->opt.maxgapiono) {
            rtk->x[j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        j = II(obs[i].sat, &rtk->opt);
        if (rtk->x[j] == 0.0) {
            freq1 = sat2freq(obs[i].sat, obs[i].code[0], nav);
            freq2 = sat2freq(obs[i].sat, obs[i].code[1], nav);
            if (obs[i].P[0] == 0.0 || obs[i].P[1] == 0.0 || freq1 == 0.0 || freq2 == 0.0) {
                continue;
            }
            ion = (obs[i].P[0] - obs[i].P[1]) / (SQR(FREQ1 / freq1) - SQR(FREQ1 / freq2));
            ecef2pos(rtk->sol.rr, pos);
            azel = rtk->ssat[obs[i].sat - 1].azel;
            ion /= ionmapf(pos, azel);
            initx(rtk, ion, SQR(rtk->opt.std[6]), j);
        } else {
            sinel = sin(MAX(rtk->ssat[obs[i].sat - 1].azel[1], 5.0 * D2R));
            rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[6] / sinel) * fabs(rtk->tt);
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk) {
    int i = ID(&rtk->opt);

    trace(3, "uddcb_ppp:\n");

    if (rtk->x[i] == 0.0) {
        initx(rtk, 1E-6, SQR(rtk->opt.std[7]), i);
    }else {
        rtk->P[i + i * rtk->nx] += SQR(rtk->opt.prn[7]) * fabs(rtk->tt);
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    double L[NFREQ], P[NFREQ], Lc, Pc, bias[MAXOBS], offset = 0.0, pos[3] = {0};
    double freq1, freq2, f2, ion, dantr[NFREQ] = {0}, dants[NFREQ] = {0};
    int i, j, k, f, sat, slip[MAXOBS] = {0}, clk_jump = 0;

    trace(3, "udbias  : n=%d\n", n);

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump = ROUND(time2gpst(obs[0].time, NULL) * 10) % 864000 == 0;
    }
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < rtk->opt.nf; j++) {
            rtk->ssat[i].slip[j] = 0;
        }
    /* detect cycle slip by LLI */
    detslp_ll(rtk, obs, n);

    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk, obs, n, nav);

    /* detect slip by Melbourne-Wubbena linear combination jump */
    detslp_mw(rtk, obs, n, nav);

    ecef2pos(rtk->sol.rr, pos);

    for (f = 0; f < NF(&rtk->opt); f++) {

        /* reset phase-bias if expire obs outage counter */
        for (i = 0; i < MAXSAT; i++) {
            if (++rtk->ssat[i].outc[f] > (uint32_t) rtk->opt.maxout || rtk->opt.modear == ARMODE_INST || clk_jump) {
                initx(rtk, 0.0, 0.0, IB(i + 1, f, &rtk->opt));
            }
        }
        for (i = k = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            j   = IB(sat, f, &rtk->opt);
            corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants, 0.0, L, P, &Lc, &Pc);

            bias[i] = 0.0;

            if (rtk->opt.ionoopt == IONOOPT_IFLC) {
                bias[i] = Lc - Pc;
                slip[i] = rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1];
            } else if (L[f] != 0.0 && P[f] != 0.0) {
                freq1   = sat2freq(sat, obs[i].code[0], nav);
                freq2   = sat2freq(sat, obs[i].code[f], nav);
                slip[i] = rtk->ssat[sat - 1].slip[f];
                if (obs[i].P[0] == 0.0 || obs[i].P[f] == 0.0 || freq1 == 0.0 || freq2 == 0.0) {
                    continue;
                }
                if (f != 0) {
                    ion = (obs[i].P[0] - obs[i].P[f]) / (1.0 - SQR(freq1 / freq2));
                } else {
                    // f == 0, using f1-f2
                    f2 = sat2freq(sat, obs[i].code[1], nav);
                    if (obs[i].P[1] == 0.0 || f2 == 0.0) {
                        continue;
                    }
                    ion = (obs[i].P[0] - obs[i].P[1]) / (1.0 - SQR(freq1 / f2));
                }
                bias[i] = L[f] - P[f] + 2.0 * ion * SQR(freq1 / freq2);
            }
            if (rtk->x[j] == 0.0 || slip[i] || bias[i] == 0.0)
                continue;

            offset += bias[i] - rtk->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */
        if (k >= 2 && fabs(offset / k) > 0.0005 * CLIGHT) {
            for (i = 0; i < MAXSAT; i++) {
                j = IB(i + 1, f, &rtk->opt);
                if (rtk->x[j] != 0.0)
                    rtk->x[j] += offset / k;
            }
            trace(2, "phase-code jump corrected: %s n=%2d dt=%12.9fs\n", time_str(rtk->sol.time, 0), k,
                  offset / k / CLIGHT);
        }
        for (i = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            j   = IB(sat, f, &rtk->opt);

            rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[8]) * fabs(rtk->tt);

            if (bias[i] == 0.0 || (rtk->x[j] != 0.0 && !slip[i]))
                continue;

            /* reinitialize phase-bias if detecting cycle slip */
            initx(rtk, bias[i], SQR(rtk->opt.std[8]), IB(sat, f, &rtk->opt));

            /* reset fix flags */
            for (k = 0; k < MAXSAT; k++)
                rtk->ambc[sat - 1].flags[k] = 0;

            trace(5, "udbias_ppp: sat=%2d bias=%.3f\n", sat, bias[i]);
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    trace(3, "udstate_ppp: n=%d\n", n);

    /* temporal update of position */
    udpos_ppp(rtk);

    /* temporal update of clock */
    udclk_ppp(rtk);

    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
        udtrop_ppp(rtk);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST) {
        udiono_ppp(rtk, obs, n, nav);
    }
    /* temporal update of L5-receiver-dcb parameters */
    if (rtk->opt.nf >= 3) {
        uddcb_ppp(rtk);
    }
    /* temporal update of phase-bias */
    udbias_ppp(rtk, obs, n, nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv, double *dant) {
    double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
    int i;

    for (i = 0; i < 3; i++) {
        ru[i] = rr[i] - rs[i];
        rz[i] = -rs[i];
    }
    if (!normv3(ru, eu) || !normv3(rz, ez))
        return;

    cosa  = dot(eu, ez, 3);
    cosa  = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
    nadir = acos(cosa);

    antmodel_s(pcv, nadir, dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos, const double *azel, const double *x, double *dtdx,
                              double *var) {
    const double zazel[] = {0.0, PI / 2.0};
    double zhd, m_h, m_w, cotz, grad_n, grad_e;

    /* zenith hydrostatic delay */
    zhd = tropmodel(time, pos, zazel, 0.0);

    /* mapping function */
    m_h = tropmapf(time, pos, azel, &m_w);

    if (azel[1] > 0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz   = 1.0 / tan(azel[1]);
        grad_n = m_w * cotz * cos(azel[0]);
        grad_e = m_w * cotz * sin(azel[0]);
        m_w += grad_n * x[1] + grad_e * x[2];
        dtdx[1] = grad_n * (x[0] - zhd);
        dtdx[2] = grad_e * (x[0] - zhd);
    }
    dtdx[0] = m_w;
    *var    = SQR(0.01);
    return m_h * zhd + m_w * (x[0] - zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, const double *x,
                      double *dtdx, const nav_t *nav, double *dtrp, double *var) {
    double trp[3] = {0};

    if (opt->tropopt == TROPOPT_SAAS) {
        *dtrp = tropmodel(time, pos, azel, REL_HUMI);
        *var  = SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
        matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);
        *dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
        return 1;
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, int sat,
                      const double *x, const nav_t *nav, double *dion, double *var) {
    if (opt->ionoopt == IONOOPT_TEC) {
        return iontec(time, nav, pos, azel, 1, dion, var);
    }
    if (opt->ionoopt == IONOOPT_BRDC) {
        *dion = ionmodel(time, nav->ion_gps, pos, azel);
        *var  = SQR(*dion * ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt == IONOOPT_EST) {
        *dion = x[II(sat, opt)];
        *var  = 0.0;
        return 1;
    }
    if (opt->ionoopt == IONOOPT_IFLC) {
        *dion = *var = 0.0;
        return 1;
    }
    return 0;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs, const double *dts, const double *var_rs,
                   const int *svh, const double *dr, int *exc, const nav_t *nav, const double *x, rtk_t *rtk, double *v,
                   double *H, double *R, double *azel) {
    prcopt_t *opt = &rtk->opt;
    double y, r, cdtr, bias, C = 0.0, wavelen=0.0, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
    double var[MAXOBS * 2 * NFREQ + MAXNVOBS], dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb, freq[NFREQ] = {0};
    double dantr[NFREQ] = {0}, dants[NFREQ] = {0}, satpco[3*NFREQ]={0}, rsat[3]={0};
    double ve[MAXOBS * 2 * NFREQ + MAXNVOBS] = {0}, vmax = 0;
    char str[32];
    int ne = 0, obsi[MAXOBS * 2 * NFREQ + MAXNVOBS] = {0}, frqi[MAXOBS * 2 * NFREQ + MAXNVOBS], maxobs, maxfrq, rej;
    int i, j, k, sat, sys, nv = 0, nx = rtk->nx, stat = 1,isys=-1,ifreq =-1,itype =-1;

    time2str(obs[0].time, str, 2);

    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < opt->nf; j++)
            rtk->ssat[i].vsat[j] = 0;

    for (i = 0; i < 3; i++)
        rr[i] = x[i] + dr[i];
    ecef2pos(rr, pos);

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        rtk->ssat[sat - 1].vs = 0;

        for(j=0;j<3;++j){
            rsat[j]= rs[i*6+j]; // has not corrected pco
        }

        if ((r = geodist(rsat, rr, e)) <= 0.0 || satazel(pos, e, azel + i * 2) < opt->elmin) {
            exc[i] = 1;
            continue;
        }
        if (!(sys = satsys(sat, NULL)) || (isys = sys2idx(sys)) < 0 || isys >= NSYS || /* !rtk->ssat[sat - 1].vs || */ // fixed: ppp will reset spp vs
            satexclude(obs[i].sat, var_rs[i], svh[i], opt) || exc[i]) {
            exc[i] = 1;
            continue;
        }
        /* tropospheric and ionospheric model */
        if (!model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart) ||
            !model_iono(obs[i].time, pos, azel + i * 2, opt, sat, x, nav, &dion, &vari)) { // L1 vertical iono delay 
            continue;
        }

        /* satellite and receiver antenna model */
        if (opt->sateph == EPHOPT_PREC) {       // sat pco
            if (opt->ionoopt == IONOOPT_IFLC) { // IF12 pco
                satantoff(obs[i].time, rsat, sat, nav, 1, satpco);
            } else { // UC pco
                satantoff(obs[i].time, rsat, sat, nav, 0, satpco);
            }
        }
        if (opt->posopt[0]) // sat pcv 
            satantpcv(rsat, rr, nav->pcvs + sat - 1, dants);
        antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr); // rcv pco pcv

        /* phase windup model */
        if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rsat, rr,
                       &rtk->ssat[sat - 1].phw)) {
            continue;
        }
        /* corrected phase and code measurements */
        corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j = 0; j < 2 * NF(opt); j++) {

            ifreq = j/2; 
            itype = j%2; // 0: phase, 1: code

            dcb = bias = 0.0;

            if (opt->ionoopt == IONOOPT_IFLC) {
                if ((y = itype == 0 ? Lc : Pc) == 0.0)
                    continue;
                for(k=0;k<3;++k){
                    rsat[k] = rs[i*6+k] + satpco[k];
                }
            } else {
                if ((y = itype == 0 ? L[ifreq] : P[ifreq]) == 0.0)
                    continue;
                if ((freq[ifreq] = sat2freq(sat, obs[i].code[ifreq], nav)) == 0.0)
                    continue;
                C = SQR(FREQ1 / freq[ifreq]) * ionmapf(pos, azel + i * 2) * (itype == 0 ? -1.0 : 1.0); // iono coef
                for(k=0;k<3;++k){
                    rsat[k] = rs[i*6+k] + satpco[k + ifreq*3];
                }
            }
            // rsat: phase center
            // update sat-rcv vector/distance for pco
            if ((r = geodist(rsat, rr, e)) <= 0.0) {
                continue;
            }

            for (k = 0; k < nx; k++)
                H[k + nx * nv] = k < 3 ? -e[k] : 0.0;

            /* receiver clock */
            k = isys;

            if(bds2or3(sat) == 2){
                H[IC(NC(opt) - 1, opt) + nx * nv] = 1.0;
                cdtr                     = x[IC(NC(opt) - 1, opt)];
            }else {
                H[IC(k, opt) + nx * nv] = 1.0;
                cdtr                    = x[IC(k, opt)];
            }

            if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
                for (k = 0; k < (opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
                    H[IT(opt) + k + nx * nv] = dtdx[k];
                }
            }
            if (opt->ionoopt == IONOOPT_EST) {
                if (rtk->x[II(sat, opt)] == 0.0)
                    continue;
                H[II(sat, opt) + nx * nv] = C;
            }
            if (ifreq == 2 && itype == 1) { /* L5-receiver-dcb */
                dcb += rtk->x[ID(opt)];
                H[ID(opt) + nx * nv] = 1.0;
            }
            if (itype == 0) { /* phase bias */
                if ((bias = x[IB(sat, ifreq, opt)]) == 0.0)
                    continue;
                H[IB(sat, ifreq, opt) + nx * nv] = 1.0;
            }
            /* residual */
            v[nv] = y - (r + cdtr - CLIGHT * dts[i * 2] + dtrp + C * dion + dcb + bias);

            if (itype == 0)
                rtk->ssat[sat - 1].resc[ifreq] = v[nv];
            else
                rtk->ssat[sat - 1].resp[ifreq] = v[nv];

            /* variance */
            var[nv] = varerr(obs[i].sat, sys, azel[1 + i * 2], ifreq, itype, opt) + vart + SQR(C) * vari + var_rs[i];
            if (sys == SYS_GLO && itype == 1)
                var[nv] += SQR(rtk->opt.std[9]);
            // adjust variance by orbtype
            switch (sat2orbtype(sat)) {
                case ORB_LEO:
                    var[nv] *= 1;
                    break;
                case ORB_MEO:
                    var[nv] *= 1;
                    break;
                case ORB_IGSO:
                    var[nv] *= 1;
                    break; 
                case ORB_GEO:
                    var[nv] *= 100;
                    break; 
            }

            trace(4, "%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, sat, itype ? "P" : "L", ifreq + 1, v[nv],
                  sqrt(var[nv]), azel[1 + i * 2] * R2D);

            /* reject satellite by pre-fit residuals */
            if (!post && opt->maxinno > 0.0 && fabs(v[nv]) > opt->maxinno) { //todo does it make sense?
                trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n", post, str, sat,
                      itype ? "P" : "L", ifreq + 1, v[nv], azel[1 + i * 2] * R2D);
                exc[i] = 1;
                rtk->ssat[sat - 1].rejc[itype]++;
                continue;
            }
            /* record large post-fit residuals */
            if (post && fabs(v[nv]) > sqrt(var[nv]) * THRES_REJECT) {
                obsi[ne] = i;
                frqi[ne] = j;
                ve[ne]   = v[nv];
                ne++;
            }
            if (itype == 0)
                rtk->ssat[sat - 1].vsat[ifreq] = 1;
            nv++;
            rtk->ssat[sat - 1].vs = 1;
        }
        // virtual observation for dion: Ion = I + noise
        if (ENABLE_VIONO) {
            if (rtk->x[II(sat, opt)] == 0.0)
                continue;
            //!! just example: ionex as vobs
            //!! warn: if VTEC is same for all sat, the H matrix is ill-conditioned
            // if (iontec(obs[i].time, nav, pos, azel + i * 2, 1, &vobs_ion, &var[nv])) { // get L1 L1 vertical iono delay vobs
            //     H[II(sat, opt) + nx * nv] = ionmapf(pos, azel + i * 2);
            //     v[nv]                     = vobs_ion - rtk->x[II(sat, opt)];
            //     nv++;
            // }
            v[nv]                     = 0.0 - rtk->x[II(sat, opt)];
            H[II(sat, opt) + nx * nv] = 1.0;
            var[nv]                   = 10000.0;
            nv++;
        }
    }
    // virtual observation for dtrp: Trp = T + noise
    if (ENABLE_VTROP) {
        v[nv]                = 0.0 - rtk->x[IT(opt)];
        H[IT(opt) + nx * nv] = 1.0;
        var[nv]              = 10000.0;
        nv++;
    }
    /* reject satellite with large and max post-fit residual */
    if (post && ne > 0) {
        vmax   = ve[0];
        maxobs = obsi[0];
        maxfrq = frqi[0];
        rej    = 0;
        for (j = 1; j < ne; j++) {
            if (fabs(vmax) >= fabs(ve[j]))
                continue;
            vmax   = ve[j];
            maxobs = obsi[j];
            maxfrq = frqi[j];
            rej    = j;
        }
        sat = obs[maxobs].sat;
        trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n", post, str, sat, maxfrq % 2 ? "P" : "L",
              maxfrq / 2 + 1, vmax, azel[1 + maxobs * 2] * R2D);
        exc[maxobs] = 1;
        rtk->ssat[sat - 1].rejc[maxfrq % 2]++;
        stat    = 0;
        ve[rej] = 0;
    }
    for (i = 0; i < nv; i++)
        for (j = 0; j < nv; j++) {
            R[i + j * nv] = i == j ? var[i] : 0.0;
        }
    return post ? stat : nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt) {
    return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat) {
    const prcopt_t *opt = &rtk->opt;
    int i, j;

    /* test # of valid satellites */
    rtk->sol.ns = 0;
    for (i = 0; i < n && i < MAXOBS; i++) {
        for (j = 0; j < opt->nf; j++) {
            if (!rtk->ssat[obs[i].sat - 1].vsat[j])
                continue;
            rtk->ssat[obs[i].sat - 1].lock[j]++;
            rtk->ssat[obs[i].sat - 1].outc[j] = 0;
            if (j == 0)
                rtk->sol.ns++;
        }
    }
    rtk->sol.stat = rtk->sol.ns < MIN_NSAT_SOL ? SOLQ_NONE : stat;

    if (rtk->sol.stat == SOLQ_FIX) {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float) rtk->Pa[i + i * rtk->na];
        }
        rtk->sol.qr[3] = (float) rtk->Pa[1];
        rtk->sol.qr[4] = (float) rtk->Pa[1 + 2 * rtk->na];
        rtk->sol.qr[5] = (float) rtk->Pa[2];
    } else {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float) rtk->P[i + i * rtk->nx];
        }
        rtk->sol.qr[3] = (float) rtk->P[1];
        rtk->sol.qr[4] = (float) rtk->P[2 + rtk->nx];
        rtk->sol.qr[5] = (float) rtk->P[2];
    }
    for (i = 0; i < NC(&rtk->opt); ++i)
        rtk->sol.dtr[i] = rtk->x[IC(i, opt)];
    for (i = 0; i < n && i < MAXOBS; i++)
        for (j = 0; j < opt->nf; j++) {
            rtk->ssat[obs[i].sat - 1].snr[j] = obs[i].SNR[j];
        }
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < opt->nf; j++) {
            if (rtk->ssat[i].slip[j] & 3)
                rtk->ssat[i].slipc[j]++;
            if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX)
                rtk->ssat[i].fix[j] = 1;
        }
}
/* test hold ambiguity -------------------------------------------------------*/
static int test_hold_amb(rtk_t *rtk) {
    int i, j, stat = 0;

    /* no fix-and-hold mode */
    if (rtk->opt.modear != ARMODE_FIXHOLD)
        return 0;

    /* reset # of continuous fixed if new ambiguity introduced */
    for (i = 0; i < MAXSAT; i++) {
        if (rtk->ssat[i].fix[0] != 2 && rtk->ssat[i].fix[1] != 2)
            continue;
        for (j = 0; j < MAXSAT; j++) {
            if (rtk->ssat[j].fix[0] != 2 && rtk->ssat[j].fix[1] != 2)
                continue;
            if (!rtk->ambc[j].flags[i] || !rtk->ambc[i].flags[j])
                stat = 1;
            rtk->ambc[j].flags[i] = rtk->ambc[i].flags[j] = 1;
        }
    }
    if (stat) {
        rtk->nfix = 0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix >= rtk->opt.minfix;
}

// static double LCfreq(const int *coef, const double* freq){
//     double f =0;
//     int i;
//     for(i=0;i<NFREQ;++i){
//         if(coef[i]==0) continue;
//         if(freq[i]==0.0) return 0.0;
//         f += coef[i]*freq[i];
//     }
//     return f; 
// }

// static double L_LCval(const int *coef, const double *v, const double* freq) {
//     double lc =0;
//     int i;
//     for(i=0;i<NFREQ;++i){
//         if(coef[i]==0) continue;
//         if(freq[i]==0.0 || v[i]==0.0) return 0.0;
//         lc += coef[i]*freq[i]*v[i]*CLIGHT/freq[i];
//     }
//     return lc;
// }

// static double P_LCval(const int *coef, const double *v, const double* freq){
//     double pc=0;
//     int i;
//         for(i=0;i<NFREQ;++i){
//         if(coef[i]==0) continue;
//         if(freq[i]==0.0 || v[i]==0.0) return 0.0;
//         pc += coef[i]*freq[i]*v[i];
//     }
//     return pc;
// }

// static double LCvar(const int *coef,const double* freq,double sig){
//     int i;
//     double a=0,b=0;
//     for(i=0;i<NFREQ;++i){
//         a += SQR(coef[i]* freq[i]);
//         b += coef[i]* freq[i];
//     }
//     return a/SQR(b) * sig;
// }

// static void average_LC(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, const double *azel) {
// #define NLC 3
//     int i, j, sat, isys;
//     double LC[4], LCv[4],LCf[4], sig, freq[NFREQ],PLC,LLC,P[NFREQ],L[NFREQ];
//     const int coef[NLC][NFREQ] = {{1, -1, 0}, {0, 1, -1}, {1, -6, 5}};
//     ambc_t *amb;
//     ssat_t * ssat;

    
//     for (i = 0; i < n; ++i) {
//         sat = obs[i].sat;
//         ssat = rtk->ssat + sat - 1;

//         if ((isys = sys2idx(satsys(sat, NULL))) < 0)
//             continue;
//         if (azel[1 + 2 * i] < rtk->opt.elmin)
//             continue;

//         for (j = 0; j < NFREQ; ++j) {
//             freq[j] = sat2freq(sat, obs[i].code[j], nav);
//             if (freq[j] == 0) {
//                 continue;
//             }
//             // OSB correction
//             P[j] = obs[i].P[j] - nav->satosb[obs->sat - 1][0][obs[i].code[j] - 1];
//             L[j] = obs[i].L[j] - nav->satosb[obs->sat - 1][1][obs[i].code[j] - 1]/(CLIGHT/freq[j]);
//         }
//         sig = sqrt(SQR(rtk->opt.err[1]) + SQR(rtk->opt.err[2] / sin(azel[1 + 2 * i])));
//         for (j = 0; j < NLC; ++j) {
//             PLC=P_LCval(coef[j], P, FreqVals[isys]);
//             LLC=L_LCval(coef[j], L, FreqVals[isys]);
//             LCf[j] =LCfreq(coef[j],freq);
//             if(PLC==0||LLC==0 || LCf[j]==0){
//                 LC[j]=0;
//                 LCv[j]=0;
//                 continue;
//             }
//             LC[j] = (LLC - PLC)/LCf[j];
//             LCv[j] = LCvar(coef[i], freq, sig * rtk->opt.eratio[0]);
//         }
//         amb = rtk->ambc + sat - 1;
//         if (ssat->slip[0] || ssat->slip[1] || ssat->slip[2] || amb->n[0] == 0.0 ||
//             fabs(timediff(amb->epoch[0], obs[0].time)) > MIN_ARC_GAP) {
//             amb->n[0] = amb->n[1] = amb->n[2] = 0.0;
//             amb->LC[0] = amb->LC[1] = amb->LC[2] = 0.0;
//             amb->LCf[0] = amb->LCf[1] = amb->LCf[2] = 0.0;
//             amb->LCv[0] = amb->LCv[1] = amb->LCv[2] = 0.0;
//             amb->fixcnt                             = 0;
//             for (j = 0; j < MAXSAT; j++)
//                 amb->flags[j] = 0;
//         }

//         for (j = 0; j < NLC; ++j) {
//             if (LC[j]) {
//                 amb->n[j] += 1.0;
//                 amb->LCf[j] = LCf[j];
//                 amb->LC[j] += (LC[j] - amb->LC[j]) / amb->n[j];
//                 amb->LCv[j] += (LCv[j] - amb->LCv[j]) / amb->n[j];
//             }
//         }
//         amb->epoch[0]=obs[0].time;
//     }

// #undef NLC
// }
// /* complementaty error function (ref [1] p.227-229) --------------------------*/
// static double q_gamma(double a, double x, double log_gamma_a);
// static double p_gamma(double a, double x, double log_gamma_a)
// {
//     double y,w;
//     int i;
    
//     if (x==0.0) return 0.0;
//     if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);
    
//     y=w=exp(a*log(x)-x-log_gamma_a)/a;
    
//     for (i=1;i<100;i++) {
//         w*=x/(a+i);
//         y+=w;
//         if (fabs(w)<1E-15) break;
//     }
//     return y;
// }
// static double q_gamma(double a, double x, double log_gamma_a)
// {
//     double y,w,la=1.0,lb=x+1.0-a,lc;
//     int i;
    
//     if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
//     w=exp(-x+a*log(x)-log_gamma_a);
//     y=w/lb;
//     for (i=2;i<100;i++) {
//         lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
//         la=lb; lb=lc;
//         w*=(i-1-a)/i;
//         y+=w/la/lb;
//         if (fabs(w/la/lb)<1E-15) break;
//     }
//     return y;
// }
// static double f_erfc(double x)
// {
//     return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
// }
// /* confidence function of integer ambiguity ----------------------------------*/
// static double conffunc(int N, double B, double sig)
// {
//     double x,p=1.0;
//     int i;
    
//     x=fabs(B-N);
//     for (i=1;i<8;i++) {
//         p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
//     }
//     return p;
// }
// /* fix wide-lane ambiguity ---------------------------------------------------*/
// static int fix_amb_WL(rtk_t *rtk, const nav_t *nav, int sat1, int sat2, int *NW)
// {
//     ambc_t *amb1,*amb2;
//     double BW,vW,confidence;
    
//     amb1=rtk->ambc+sat1-1;
//     amb2=rtk->ambc+sat2-1;
//     if (!amb1->n[0] || !amb2->n[0])
//         return 0;
//     if(amb1->LCf[0] != amb2->LCf[0])
//         return 0;

//     /* wide-lane ambiguity */
//     BW=(amb1->LC[0]-amb2->LC[0])/(CLIGHT/amb1->LCf[0]);
//     *NW=ROUND(BW);
    
//     /* variance of wide-lane ambiguity */
//     vW=(amb1->LCv[0]/amb1->n[0]+amb2->LCv[0]/amb2->n[0])/SQR((CLIGHT/amb1->LCf[0]));
//     confidence = conffunc(*NW,BW,sqrt(vW));
//     /* validation of integer wide-lane ambigyity */
//     return fabs(*NW-BW)<=rtk->opt.thresar[2]&&
//            confidence>=rtk->opt.thresar[1];
// }

// static int testsat(rtk_t *rtk, int sat1, int sat2, int *sys1, int *sys2) {
//     int s1, s2;

//     if ((s1 = satsys(sat1, NULL)) == SYS_NONE || (s2 = satsys(sat2, NULL)) == SYS_NONE) {
//         return 0;
//     }
//     if (!rtk->ssat[sat1 - 1].vsat[0] || !rtk->ssat[sat2 - 1].vsat[0]) {
//         return 0;
//     }

//     if (sys1) {
//         *sys1 = s1;
//     }
//     if (sys2) {
//         *sys2 = s2;
//     }
//     if (s1 != s2) {
//         return 0;
//     }

//     if (s1 == SYS_GPS) {
//         return 1;
//     }
//     if (s1 == SYS_GLO) {
//         return rtk->opt.glomodear;
//     }
//     if (s2 == SYS_GAL) {
//         return rtk->opt.galmodear;
//     }
//     if (s1 == SYS_CMP) {
//         s1 = bds2or3(sat1);
//         s2 = bds2or3(sat2);
//         if (s1 != s2) {
//             return 0;
//         }
//         if (s1 == 2) {
//             return rtk->opt.bdsmodear & 1;
//         }
//         if (s1 == 3) {
//             return rtk->opt.bdsmodear & 2;
//         }
//     }
//     return 0;
// }

// /* linear dependency check ---------------------------------------------------*/
// static int is_depend(int sat1, int sat2, int *flgs, int *max_flg)
// {
//     int i;
    
//     if (flgs[sat1-1]==0&&flgs[sat2-1]==0) {
//         flgs[sat1-1]=flgs[sat2-1]=++(*max_flg);
//     }
//     else if (flgs[sat1-1]==0&&flgs[sat2-1]!=0) {
//         flgs[sat1-1]=flgs[sat2-1];
//     }
//     else if (flgs[sat1-1]!=0&&flgs[sat2-1]==0) {
//         flgs[sat2-1]=flgs[sat1-1];
//     }
//     else if (flgs[sat1-1]>flgs[sat2-1]) {
//         for (i=0;i<MAXSAT;i++) if (flgs[i]==flgs[sat2-1]) flgs[i]=flgs[sat1-1];
//     }
//     else if (flgs[sat1-1]<flgs[sat2-1]) {
//         for (i=0;i<MAXSAT;i++) if (flgs[i]==flgs[sat1-1]) flgs[i]=flgs[sat2-1];
//     }
//     else return 0; /* linear depenent */
//     return 1;
// }

// /* select fixed ambiguities --------------------------------------------------*/
// static int sel_amb(int *sat1, int *sat2, double *N, double *var, int n)
// {
//     int i,j,flgs[MAXSAT]={0},max_flg=0;
    
//     /* sort by variance */
//     for (i=0;i<n;i++) for (j=1;j<n-i;j++) {
//         if (var[j]>=var[j-1]) continue;
//         SWAP_I(sat1[j],sat1[j-1]);
//         SWAP_I(sat2[j],sat2[j-1]);
//         SWAP_D(N[j],N[j-1]);
//         SWAP_D(var[j],var[j-1]);
//     }
//     /* select linearly independent satellite pair */
//     for (i=j=0;i<n;i++) {
//         if (!is_depend(sat1[i],sat2[i],flgs,&max_flg)) continue;
//         sat1[j]=sat1[i];
//         sat2[j]=sat2[i];
//         N[j]=N[i];
//         var[j++]=var[i];
//     }
//     return j;
// }

// /* fixed solution ------------------------------------------------------------*/
// static int fix_sol(rtk_t *rtk, const int *sat1, const int *sat2,
//                    const double *NC, int n)
// {
//     double *v,*H,*R;
//     int i,j,k,info;
    
//     if (n<=0) return 0;
    
//     v=zeros(n,1); H=zeros(rtk->nx,n); R=zeros(n,n);
    
//     /* constraints to fixed ambiguities */
//     for (i=0;i<n;i++) {
//         j=IB(sat1[i],0,&rtk->opt);
//         k=IB(sat2[i],0,&rtk->opt);
//         v[i]=NC[i]-(rtk->x[j]-rtk->x[k]);
//         H[j+i*rtk->nx]= 1.0;
//         H[k+i*rtk->nx]=-1.0;
//         R[i+i*n]=SQR(CONST_AMB);
//     }
//     /* update states with constraints */
//     if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,n))) {
//         trace(1,"filter error (info=%d)\n",info);
//         free(v); free(H); free(R);
//         return 0;
//     }
//     /* set solution */
//     for (i=0;i<rtk->na;i++) {
//         rtk->xa[i]=rtk->x[i];
//         for (j=0;j<rtk->na;j++) {
//             rtk->Pa[i+j*rtk->na]=rtk->Pa[j+i*rtk->na]=rtk->P[i+j*rtk->nx];
//         }
//     }
//     /* set flags */
//     for (i=0;i<n;i++) {
//         rtk->ambc[sat1[i]-1].flags[sat2[i]-1]=1;
//         rtk->ambc[sat2[i]-1].flags[sat1[i]-1]=1;
//     }
//     free(v); free(H); free(R);
//     return 1;
// }

// /* fix narrow-lane ambiguity by ILS ------------------------------------------*/
// static int fix_amb_ROUND(rtk_t *rtk, int *sat1, int *sat2, int *NW, int n)
// {
//     double C1,C2,B1,v1,BC,v,vc,*NC,*var,lam_NL;
//     int i,j,k,m=0,N1,stat;

//     int sys1,sys2,isys;
//     const double *freq; 
//     const int coef[][NFREQ]={{1,1}};
//     double confidence;

//     NC=zeros(n,1); var=zeros(n,1);

//     for (i = 0; i < n; i++) {
//         j = IB(sat1[i], 0, &rtk->opt);
//         k = IB(sat2[i], 0, &rtk->opt);

//         sys1 = satsys(sat1[i],NULL);
//         sys2 = satsys(sat2[i],NULL);
//         assert(sys1==sys2 && sys1!=SYS_NONE);
//         isys=(sys2idx(sys1));
//         assert(isys>=0);
//         freq = FreqVals[isys];
//         assert(freq[0] && freq[1]);

//         lam_NL = CLIGHT / (LCfreq(coef[0], freq));
//         C1     = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
//         C2     = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

//         /* narrow-lane ambiguity */
//         B1 = (rtk->x[j] - rtk->x[k] + C2 * (CLIGHT / freq[1]) * NW[i]) / lam_NL;
//         N1 = ROUND(B1);

//         /* variance of narrow-lane ambiguity */
//         var[m]=rtk->P[j+j*rtk->nx]+rtk->P[k+k*rtk->nx]-2.0*rtk->P[j+k*rtk->nx];
//         v1=var[m]/SQR(lam_NL);
//         confidence = conffunc(N1,B1,sqrt(v1));
//         /* validation of narrow-lane ambiguity */
//         if (fabs(N1-B1)>rtk->opt.thresar[2]||
//             confidence<rtk->opt.thresar[1]) {
//             continue;
//         }

//         /* iono-free ambiguity (m) */
//         BC=C1*(CLIGHT/freq[0])*N1+C2*(CLIGHT/freq[1])*(N1-NW[i]);

//         /* check residuals */
//         v=rtk->ssat[sat1[i]-1].resc[0]-rtk->ssat[sat2[i]-1].resc[0];
//         vc=v+(BC-(rtk->x[j]-rtk->x[k]));
//         if (fabs(vc)>THRES_RES) continue;

//         sat1[m]=sat1[i];
//         sat2[m]=sat2[i];
//         NC[m++]=BC;
//     }

//    /* select fixed ambiguities by dependancy check */
//     m=sel_amb(sat1,sat2,NC,var,m);

//     /* fixed solution */
//     stat=fix_sol(rtk,sat1,sat2,NC,m);

//     free(NC); free(var);
//     return stat&&m>=3;
// }

// /* fix narrow-lane ambiguity by ILS ------------------------------------------*/
// static int fix_amb_ILS(rtk_t *rtk, int *sat1, int *sat2, int *NW, int n)
// {
//     double C1,C2,*B1,*N1,*NC,*D,*E,*Q,s[2],lam_NL;
//     int i,j,k,m=0,info,stat,flgs[MAXSAT]={0},max_flg=0;

//     int sys1,sys2,isys;
//     const double *freq; 
//     const int coef[][NFREQ]={{1,1}};
    
//     B1=zeros(n,1); N1=zeros(n,2); D=zeros(rtk->nx,n); E=mat(n,rtk->nx);
//     Q=mat(n,n); NC=mat(n,1);
    
//     for (i=0;i<n;i++) {
        
//         /* check linear independency */
//         if (!is_depend(sat1[i],sat2[i],flgs,&max_flg)) continue;
        
//         j=IB(sat1[i],0,&rtk->opt);
//         k=IB(sat2[i],0,&rtk->opt);

//         sys1 = satsys(sat1[i], NULL);
//         sys2 = satsys(sat2[i], NULL);
//         assert(sys1 == sys2 && sys1 != SYS_NONE);
//         isys = (sys2idx(sys1));
//         assert(isys >= 0);
//         freq = FreqVals[isys];
//         assert(freq[0] && freq[1]);

//         lam_NL = CLIGHT / (LCfreq(coef[0], freq));
//         C1     = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
//         C2     = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

//         /* float narrow-lane ambiguity (cycle) */
//         B1[m]=(rtk->x[j]-rtk->x[k]+C2*(CLIGHT/freq[1])*NW[i])/lam_NL;
//         N1[m]=ROUND(B1[m]);
        
//         /* validation of narrow-lane ambiguity */
//         if (fabs(N1[m]-B1[m])>rtk->opt.thresar[2]) continue;
        
//         /* narrow-lane ambiguity transformation matrix */
//         D[j+m*rtk->nx]= 1.0/lam_NL;
//         D[k+m*rtk->nx]=-1.0/lam_NL;
        
//         sat1[m]=sat1[i];
//         sat2[m]=sat2[i];
//         NW[m++]=NW[i];
//     }
//     if (m<3) 
//         return 0;
    
//     /* covariance of narrow-lane ambiguities */
//     matmul("TN",m,rtk->nx,rtk->nx,1.0,D,rtk->P,0.0,E);
//     matmul("NN",m,m,rtk->nx,1.0,E,D,0.0,Q);
    
//     /* integer least square */
//     if ((info=lambda(m,2,B1,Q,N1,s))) {
//         trace(2,"lambda error: info=%d\n",info);
//         return 0;
//     }
//     if (s[0]<=0.0) return 0;
    
//     rtk->sol.ratio=(float)(MIN(s[1]/s[0],999.9));
    
//     /* varidation by ratio-test */
//     if (rtk->opt.thresar[0]>0.0&&rtk->sol.ratio<rtk->opt.thresar[0]) {
//         trace(2,"varidation error: n=%2d ratio=%8.3f\n",m,rtk->sol.ratio);
//         return 0;
//     }
//     trace(2,"varidation ok: %s n=%2d ratio=%8.3f\n",time_str(rtk->sol.time,0),m,
//           rtk->sol.ratio);
    
//     /* narrow-lane to iono-free ambiguity */
//     for (i=0;i<m;i++) {

//         j=IB(sat1[i],0,&rtk->opt);
//         k=IB(sat2[i],0,&rtk->opt);

//         sys1 = satsys(sat1[i], NULL);
//         sys2 = satsys(sat2[i], NULL);
//         assert(sys1 == sys2 && sys1 != SYS_NONE);
//         isys = (sys2idx(sys1));
//         assert(isys >= 0);
//         freq = FreqVals[isys];
//         assert(freq[0] && freq[1]);

//         lam_NL = CLIGHT / (LCfreq(coef[0], freq));
//         C1     = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
//         C2     = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

//         NC[i]=C1*(CLIGHT/freq[0])*N1[i]+C2*(CLIGHT/freq[1])*(N1[i]-NW[i]);
//     }
//     /* fixed solution */
//     stat=fix_sol(rtk,sat1,sat2,NC,m);
    
//     free(B1); free(N1); free(D); free(E); free(Q); free(NC);
    
//     return stat;
// }

/* ambiguity resolution in ppp -----------------------------------------------*/
// static int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc, const nav_t *nav, const double *azel, double *x,
//                   double *P) {
//     return 0;

//     double elmask;
//     int i, j, m = 0, stat = 0, *NW, *sat1, *sat2 ,sys1,sys2;

//     if(rtk->opt.modear == ARMODE_OFF){
//         return 0;
//     }

//     if (n <= 0 || rtk->opt.ionoopt != IONOOPT_IFLC || rtk->opt.nf < 2)
//         return 0;

//     trace(3, "ppp_ar: time=%s n=%d\n", time_str(obs[0].time, 0), n);

//     elmask = rtk->opt.elmaskar > 0.0 ? rtk->opt.elmaskar : rtk->opt.elmin;

//     sat1 = imat(n * n, 1); 
//     sat2 = imat(n * n, 1);
//     NW   = imat(n * n, 1);

//     memset(sat1,0,n*n*sizeof(int));
//     memset(sat2,0,n*n*sizeof(int));

//     average_LC(rtk,obs,n,nav,azel);

//     // wide-lane
//     for (i = 0; i < n - 1; i++) {
//         for (j = i + 1; j < n; j++) {
//             if (!testsat(rtk, obs[i].sat, obs[j].sat, &sys1, &sys2)) {
//                 continue;
//             }
//             if (azel[1 + i * 2] < elmask || azel[1 + j * 2] < elmask)
//                 continue;
// #if 0
//         /* test already fixed */
//         if (rtk->ambc[obs[i].sat-1].flags[obs[j].sat-1]&&
//             rtk->ambc[obs[j].sat-1].flags[obs[i].sat-1]) continue;
// #endif
//             if (fix_amb_WL(rtk, nav, obs[i].sat, obs[j].sat, NW + m)) {
//                 sat1[m] = obs[i].sat;
//                 sat2[m] = obs[j].sat;
//                 m++;
//             }
//         }
//     }
//     /* fix narrow-lane ambiguity */
//     if (rtk->opt.methodar == ARMETHOD_ROUND) {
//         stat=fix_amb_ROUND(rtk,sat1,sat2,NW,m);
//     }else if(rtk->opt.methodar == ARMETHOD_ILS){
//         stat=fix_amb_ILS(rtk,sat1,sat2,NW,m);
//     }
//     free(sat1); free(sat2); free(NW);
    
//     return stat;
// }

static int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc, const nav_t *nav, const double *azel, double *x,
                  double *P) {
    return 0;
}

/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    const prcopt_t *opt = &rtk->opt;
    double *rs, *dts, *var, *v, *H, *R, *azel, *xp, *Pp, dr[3] = {0}, std[3];
    char str[32];
    int i, j, nv, info, svh[MAXOBS], exc[MAXOBS] = {0}, stat = SOLQ_SINGLE;

    time2str(obs[0].time, str, 2);
    trace(3, "pppos   : time=%s nx=%d n=%d\n", str, rtk->nx, n);

    rs   = mat(6, n);
    dts  = mat(2, n);
    var  = mat(1, n);
    azel = zeros(2, n);

    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < opt->nf; j++)
            rtk->ssat[i].fix[j] = 0;

    /* temporal update of ekf states */
    udstate_ppp(rtk, obs, n, nav);

    /* satellite positions and clocks */
    satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        testeclipse(obs, n, nav, rs);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp, opt->odisp[0], dr);
    }
    nv = n * rtk->opt.nf * 2 + MAXSAT + 3 + MAXNVOBS;
    xp = mat(rtk->nx, 1);
    Pp = zeros(rtk->nx, rtk->nx);
    v  = mat(nv, 1);
    H  = mat(rtk->nx, nv);
    R  = mat(nv, nv);

    for (i = 0; i < rtk->opt.niter; i++) {

        matcpy(xp, rtk->x, rtk->nx, 1);
        matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

        /* prefit residuals */
        if (!(nv = ppp_res(0, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel))) {
            trace(2, "%s ppp (%d) no valid obs data\n", str, i + 1);
            break;
        }
        // tracemat(1,H,rtk->nx, nv,20,10);
        /* measurement update of ekf states */
        if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
            trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
            break;
        }
        /* postfit residuals */
        if (ppp_res(i + 1, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)) {
            matcpy(rtk->x, xp, rtk->nx, 1);
            matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
            stat = SOLQ_PPP;
            break;
        }
    }
    if (i >= rtk->opt.niter) {
        trace(2, "%s ppp (%d) iteration overflows\n", str, i);
    }
    if (stat == SOLQ_PPP) {

        /* ambiguity resolution in ppp */
        if (ppp_ar(rtk, obs, n, exc, nav, azel, xp, Pp) &&
            ppp_res(-1, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)) {

            matcpy(rtk->xa, xp, rtk->nx, 1);
            matcpy(rtk->Pa, Pp, rtk->nx, rtk->nx);

            for (i = 0; i < 3; i++)
                std[i] = sqrt(Pp[i + i * rtk->nx]);
            if (norm(std, 3) < MAX_STD_FIX)
                stat = SOLQ_FIX;
        } else {
            rtk->nfix = 0;
        }
        /* update solution status */
        update_stat(rtk, obs, n, stat);

        /* hold fixed ambiguities */
        if (stat == SOLQ_FIX && test_hold_amb(rtk)) {
            matcpy(rtk->x, xp, rtk->nx, 1);
            matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
            trace(2, "%s hold ambiguity\n", str);
            rtk->nfix = 0;
        }
    }
    free(rs);
    free(dts);
    free(var);
    free(azel);
    free(xp);
    free(Pp);
    free(v);
    free(H);
    free(R);
}
