/*------------------------------------------------------------------------------
 * preceph.c : precise ephemeris and clock functions
 *
 *          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
 *
 * references :
 *     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
 *         12 February, 2007
 *     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
 *         27 August, 1998
 *     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
 *     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
 *         Space Technology Library, 2004
 *     [5] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-d),
 *         February 21, 2016
 *
 * version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
 * history : 2009/01/18 1.0  new
 *           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
 *           2009/05/15 1.2  support glonass,galileo,qzs
 *           2009/12/11 1.3  support wild-card expansion of file path
 *           2010/07/21 1.4  added api:
 *                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
 *                               readdcb()
 *                           changed api:
 *                               readsp3()
 *                           deleted api:
 *                               eph2posp()
 *           2010/09/09 1.5  fix problem when precise clock outage
 *           2011/01/23 1.6  support qzss satellite code
 *           2011/09/12 1.7  fix problem on precise clock outage
 *                           move sunmmonpos() to rtkcmn.c
 *           2011/12/01 1.8  modify api readsp3()
 *                           precede later ephemeris if ephemeris is NULL
 *                           move eci2ecef() to rtkcmn.c
 *           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
 *           2013/11/20 1.10 modify option for api readsp3()
 *           2014/04/03 1.11 accept extenstion including sp3,eph,SP3,EPH
 *           2014/05/23 1.12 add function to read sp3 velocity records
 *                           change api: satantoff()
 *           2014/08/31 1.13 add member cov and vco in peph_t sturct
 *           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
 *           2015/05/10 1.15 add api readfcb()
 *                           modify api readdcb()
 *           2017/04/11 1.16 fix bug on antenna offset correction in peph2pos()
 *           2020/11/30 1.17 support SP3-d [5] to accept more than 85 satellites
 *                           support NavIC/IRNSS in API peph2pos()
 *                           LC defined GPS/QZS L1-L2, GLO G1-G2, GAL E1-E5b,
 *                            BDS B1I-B2I and IRN L5-S for API satantoff()
 *                           fix bug on reading SP3 file extension
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define NMAX 10         /* order of polynomial interpolation */
#define MAXDTE 900.0    /* max time difference to ephem time (s) */
#define EXTERR_CLK 1E-3 /* extrapolation error for clock (m/s) */
#define EXTERR_EPH 5E-7 /* extrapolation error for ephem (m/s^2) */

/* read SP3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats, double *bfact, char *tsys) {
    int i, j, k = 0, ns = 0, sys, prn;
    char buff[1024];

    trace(3, "readsp3h:\n");

    tsys[0] = '\0';

    for (i = 0;; i++) {
        if (!fgets(buff, sizeof(buff), fp))
            break;

        if (i == 0) {
            *type = buff[2];
            if (str2time(buff, 3, 28, time))
                return 0;
        } else if (!strncmp(buff, "+ ", 2)) { /* satellite id */
            if (ns == 0) {
                ns = (int) str2num(buff, 3, 3);
            }
            for (j = 0; j < 17 && k < ns; j++) {
                sys = code2sys(buff[9 + 3 * j]);
                prn = (int) str2num(buff, 10 + 3 * j, 2);
                if (k < MAXSAT)
                    sats[k++] = satno(sys, prn);
            }
        } else if (!strncmp(buff, "++", 2)) { /* orbit accuracy */
            continue;
        } else if (!strncmp(buff, "%c", 2) && tsys[0] == '\0') { /* time system */
            strncpy(tsys, buff + 9, 3);
            tsys[3] = '\0';
        } else if (!strncmp(buff, "%f", 2) && bfact[0] == 0.0) { /* fp base number */
            bfact[0] = str2num(buff, 3, 10);
            bfact[1] = str2num(buff, 14, 12);
        } else if (!strncmp(buff, "%i", 2)) {
            continue;
        } else if (!strncmp(buff, "/*", 2)) { /* comment */
            continue;
        } else if (!strncmp(buff, "* ", 2)) { /* first record */
            /* roll back file pointer */
            fseek(fp, -(long) strlen(buff), SEEK_CUR);
            break;
        }
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph) {
    peph_t *nav_peph;

    if (nav->ne >= nav->nemax) {
        nav->nemax += 256;
        if (!(nav_peph = (peph_t *) realloc(nav->peph, sizeof(peph_t) * nav->nemax))) {
            trace(1, "readsp3b malloc error n=%d\n", nav->nemax);
            free(nav->peph);
            nav->peph = NULL;
            nav->ne = nav->nemax = 0;
            return 0;
        }
        nav->peph = nav_peph;
    }
    nav->peph[nav->ne++] = *peph;
    return 1;
}
/* read SP3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact, char *tsys, int index, int opt,
                     nav_t *nav) {
    peph_t peph;
    gtime_t time;
    double val, std, base;
    int i, j, sat, sys, prn, n = ns * (type == 'P' ? 1 : 2), pred_o, pred_c, v;
    char buff[1024];

    trace(3, "readsp3b: type=%c ns=%d index=%d opt=%d\n", type, ns, index, opt);

    while (fgets(buff, sizeof(buff), fp)) {

        if (!strncmp(buff, "EOF", 3))
            break;

        if (buff[0] != '*' || str2time(buff, 3, 28, &time)) {
            trace(2, "sp3 invalid epoch %31.31s\n", buff);
            continue;
        }
        if (!strcmp(tsys, "UTC"))
            time = utc2gpst(time); /* utc->gpst */
        peph.time  = time;
        peph.index = index;

        for (i = 0; i < MAXSAT; i++) {
            for (j = 0; j < 4; j++) {
                peph.pos[i][j] = 0.0;
                peph.std[i][j] = 0.0f;
                peph.vel[i][j] = 0.0;
                peph.vst[i][j] = 0.0f;
            }
            for (j = 0; j < 3; j++) {
                peph.cov[i][j] = 0.0f;
                peph.vco[i][j] = 0.0f;
            }
        }
        for (i = pred_o = pred_c = v = 0; i < n && fgets(buff, sizeof(buff), fp); i++) {

            if (strlen(buff) < 4 || (buff[0] != 'P' && buff[0] != 'V'))
                continue;

            sys = buff[1] == ' ' ? SYS_GPS : code2sys(buff[1]);
            prn = (int) str2num(buff, 2, 2);
            if (sys == SYS_SBS)
                prn += 100;
            else if (sys == SYS_QZS)
                prn += 192; /* extension to sp3-c */

            if (!(sat = satno(sys, prn)))
                continue;

            if (buff[0] == 'P') {
                pred_c = strlen(buff) >= 76 && buff[75] == 'P';
                pred_o = strlen(buff) >= 80 && buff[79] == 'P';
            }
            for (j = 0; j < 4; j++) {

                /* read option for predicted value */
                if (j < 3 && (opt & 1) && pred_o)
                    continue;
                if (j < 3 && (opt & 2) && !pred_o)
                    continue;
                if (j == 3 && (opt & 1) && pred_c)
                    continue;
                if (j == 3 && (opt & 2) && !pred_c)
                    continue;

                val = str2num(buff, 4 + j * 14, 14);
                std = str2num(buff, 61 + j * 3, j < 3 ? 2 : 3);

                if (buff[0] == 'P') { /* position */
                    if (val != 0.0 && fabs(val - 999999.999999) >= 1E-6) {
                        peph.pos[sat - 1][j] = val * (j < 3 ? 1000.0 : 1E-6);
                        v                    = 1; /* valid epoch */
                    }
                    if ((base = bfact[j < 3 ? 0 : 1]) > 0.0 && std > 0.0) {
                        peph.std[sat - 1][j] = (float) (pow(base, std) * (j < 3 ? 1E-3 : 1E-12));
                    }
                } else if (v) { /* velocity */
                    if (val != 0.0 && fabs(val - 999999.999999) >= 1E-6) {
                        peph.vel[sat - 1][j] = val * (j < 3 ? 0.1 : 1E-10);
                    }
                    if ((base = bfact[j < 3 ? 0 : 1]) > 0.0 && std > 0.0) {
                        peph.vst[sat - 1][j] = (float) (pow(base, std) * (j < 3 ? 1E-7 : 1E-16));
                    }
                }
            }
        }
        if (v) {
            if (!addpeph(nav, &peph))
                return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2) {
    peph_t *q1 = (peph_t *) p1, *q2 = (peph_t *) p2;
    double tt = timediff(q1->time, q2->time);
    return tt < -1E-9 ? -1 : (tt > 1E-9 ? 1 : q1->index - q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav, int opt) {
    int i, j, k, m;

    trace(3, "combpeph: ne=%d\n", nav->ne);

    qsort(nav->peph, nav->ne, sizeof(peph_t), cmppeph);

    if (opt & 4)
        return;

    for (i = 0, j = 1; j < nav->ne; j++) {

        if (fabs(timediff(nav->peph[i].time, nav->peph[j].time)) < 1E-9) {

            for (k = 0; k < MAXSAT; k++) {
                if (norm(nav->peph[j].pos[k], 4) <= 0.0)
                    continue;
                for (m = 0; m < 4; m++)
                    nav->peph[i].pos[k][m] = nav->peph[j].pos[k][m];
                for (m = 0; m < 4; m++)
                    nav->peph[i].std[k][m] = nav->peph[j].std[k][m];
                for (m = 0; m < 4; m++)
                    nav->peph[i].vel[k][m] = nav->peph[j].vel[k][m];
                for (m = 0; m < 4; m++)
                    nav->peph[i].vst[k][m] = nav->peph[j].vst[k][m];
            }
        } else if (++i < j)
            nav->peph[i] = nav->peph[j];
    }
    nav->ne = i + 1;

    trace(4, "combpeph: ne=%d\n", nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
 * read sp3 precise ephemeris/clock files and set them to navigation data
 * args   : char   *file       I   sp3-c precise ephemeris file
 *                                 (wind-card * is expanded)
 *          nav_t  *nav        IO  navigation data
 *          int    opt         I   options (1: only observed + 2: only predicted +
 *                                 4: not combined)
 * return : none
 * notes  : see ref [1]
 *          precise ephemeris is appended and combined
 *          nav->peph and nav->ne must by properly initialized before calling the
 *          function
 *          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
 *-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt) {
    FILE *fp;
    gtime_t time    = {0};
    double bfact[2] = {0};
    int i, j, n, ns, sats[MAXSAT] = {0};
    char *efiles[MAXEXFILE], *ext, type = ' ', tsys[4] = "";

    trace(3, "readpephs: file=%s\n", file);

    for (i = 0; i < MAXEXFILE; i++) {
        if (!(efiles[i] = (char *) malloc(1024))) {
            for (i--; i >= 0; i--)
                free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n = expath(file, efiles, MAXEXFILE);

    for (i = j = 0; i < n; i++) {
        if (!(ext = strrchr(efiles[i], '.')))
            continue;

        if (!strstr(ext, ".sp3") && !strstr(ext, ".SP3") && !strstr(ext, ".eph") && !strstr(ext, ".EPH"))
            continue;

        if (!(fp = fopen(efiles[i], "r"))) {
            trace(2, "sp3 file open error %s\n", efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns = readsp3h(fp, &time, &type, sats, bfact, tsys);

        /* read sp3 body */
        readsp3b(fp, type, sats, ns, bfact, tsys, j++, opt, nav);

        fclose(fp);
    }
    for (i = 0; i < MAXEXFILE; i++)
        free(efiles[i]);

    /* combine precise ephemeris */
    if (nav->ne > 0)
        combpeph(nav, opt);
}
/* read satellite antenna parameters -------------------------------------------
 * read satellite antenna parameters
 * args   : char   *file       I   antenna parameter file
 *          gtime_t time       I   time
 *          nav_t  *nav        IO  navigation data
 * return : status (1:ok,0:error)
 * notes  : only support antex format for the antenna parameter file
 *-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav) {
    pcvs_t pcvs = {0};
    pcv_t pcv0  = {0}, *pcv;
    int i;

    trace(3, "readsap : file=%s time=%s\n", file, time_str(time, 0));

    if (!readpcv(file, &pcvs))
        return 0;

    for (i = 0; i < MAXSAT; i++) {
        pcv          = searchpcv(i + 1, "", time, &pcvs);
        nav->pcvs[i] = pcv ? *pcv : pcv0;
    }
    free(pcvs.pcv);
    return 1;
}

static int readosbf(const char *file, nav_t *nav, const sta_t *sta) {
    FILE *fp;
    char buff[1024], tobs = '\0';
    int scope = 0, year = 9, doy = 0, sod = 0, isat = 0, idx1 = CODE_NONE;
    double epoch[6] = {0, 1, 1, 0, 0, 0}, unit = 0, bias = 0, std = 0;

    if (!(fp = fopen(file, "r"))) {
        trace(2, "osb parameters file open error: %s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {
        if (0 == strncmp(buff, "+BIAS/SOLUTION", 14)) {
            scope = 1; // bias data
            continue;
        } else if (0 == strncmp(buff, "-BIAS/SOLUTION", 14)) {
            scope = 0;
            break;
        }

        if (scope == 0) {
            continue;
        }

        if (buff[0] == '*') {
            continue;
        }

        if (0 == strncmp(buff + 1, "OSB", 3)) {
            tobs     = '\0';
            unit     = 0;
            bias     = 0;
            std      = 0;
            epoch[0] = 0;

            isat = satid2no(buff + 11);
            if (isat == 0) {
                continue;
            }
            idx1 = obs2code(buff + 25 + 1);
            if (idx1 == CODE_NONE) {
                trace(4, "osb parameters code not supported yet: %2s\n", buff + 25 + 1);
                continue;
            }
            tobs = buff[25];
            if (tobs != 'C' && tobs != 'L') {
                trace(4, "osb parameters tobs not supported yet: %1s\n", buff + 25);
                continue;
            }

            if (3 > sscanf(buff + 35, "%4d:%3d:%5d", &year, &doy, &sod)) {
                trace(4, "osb parameters time read error: %14s\n", buff + 35);
                continue;
            }

            if (3 > sscanf(buff + 50, "%4d:%3d:%5d", &year, &doy, &sod)) {
                trace(4, "osb parameters time read error: %14s\n", buff + 50);
                continue;
            }

            if (0 == strncasecmp(buff + 65, "ns", 2)) {
                unit = CLIGHT / 1e9;
            } else if (0 == strncasecmp(buff + 65, "cyc", 3)) {
                trace(4, "osb parameters unit not supported yet: %4s\n", buff + 65);
                continue;
            } else {
                trace(4, "osb parameters unit not supported yet: %4s\n", buff + 65);
                continue;
            }

            bias = str2num(buff, 70, 21);
            std  = str2num(buff, 92, 11);

            if (tobs == 'C') {
                nav->satosb[isat - 1][0][idx1 - 1] = bias * unit;
            } else if (tobs == 'L') {
                nav->satosb[isat - 1][1][idx1 - 1] = bias * unit;
            } else {

            }
        } else if (0 == strncmp(buff + 1, "DCB", 3)) {

        } else {
        }
    }

    return 1;
}

extern int readosb(const char *file, nav_t *nav, const sta_t *sta) {
    int i, n;
    char *efiles[MAXEXFILE] = {0};

    trace(3, "readosb : file=%s\n", file);

    memset(nav->satosb, 0, sizeof(nav->satosb));
    memset(nav->rcvosb, 0, sizeof(nav->rcvosb));

    for (i = 0; i < MAXEXFILE; i++) {
        if (!(efiles[i] = (char *) malloc(1024))) {
            for (i--; i >= 0; i--)
                free(efiles[i]);
            return 0;
        }
    }
    n = expath(file, efiles, MAXEXFILE);

    for (i = 0; i < n; i++) {
        readosbf(efiles[i], nav, sta);
    }
    for (i = 0; i < MAXEXFILE; i++)
        free(efiles[i]);

    return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n) {
    int i, j;

    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            y[i] = (x[i + j] * y[i] - x[i] * y[i + 1]) / (x[i + j] - x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs, double *dts, double *vare, double *varc) {
    double t[NMAX + 1], p[3][NMAX + 1], c[2], *pos, std = 0.0, s[3], sinl, cosl;
    int i, j, k, index;

    trace(4, "pephpos : time=%s sat=%2d\n", time_str(time, 3), sat);

    rs[0] = rs[1] = rs[2] = dts[0] = 0.0;

    if (nav->ne < NMAX + 1 || timediff(time, nav->peph[0].time) < -MAXDTE ||
        timediff(time, nav->peph[nav->ne - 1].time) > MAXDTE) {
        trace(3, "no prec ephem %s sat=%2d\n", time_str(time, 0), sat);
        return 0;
    }
    /* binary search */
    for (i = 0, j = nav->ne - 1; i < j;) {
        k = (i + j) / 2;
        if (timediff(nav->peph[k].time, time) < 0.0)
            i = k + 1;
        else
            j = k;
    }
    index = i <= 0 ? 0 : i - 1;

    /* polynomial interpolation for orbit */
    i = index - (NMAX + 1) / 2;
    if (i < 0)
        i = 0;
    else if (i + NMAX >= nav->ne)
        i = nav->ne - NMAX - 1;

    for (j = 0; j <= NMAX; j++) {
        t[j] = timediff(nav->peph[i + j].time, time);
        if (norm(nav->peph[i + j].pos[sat - 1], 3) <= 0.0) {
            trace(3, "prec ephem outage %s sat=%2d\n", time_str(time, 0), sat);
            return 0;
        }
    }
    for (j = 0; j <= NMAX; j++) {
        pos = nav->peph[i + j].pos[sat - 1];
        /* correciton for earh rotation ver.2.4.0 */
        sinl    = sin(OMGE * t[j]);
        cosl    = cos(OMGE * t[j]);
        p[0][j] = cosl * pos[0] - sinl * pos[1];
        p[1][j] = sinl * pos[0] + cosl * pos[1];
        p[2][j] = pos[2];
    }
    for (i = 0; i < 3; i++) {
        rs[i] = interppol(t, p[i], NMAX + 1);
    }
    if (vare) {
        for (i = 0; i < 3; i++)
            s[i] = nav->peph[index].std[sat - 1][i];
        std = norm(s, 3);

        /* extrapolation error for orbit */
        if (t[0] > 0.0)
            std += EXTERR_EPH * SQR(t[0]) / 2.0;
        else if (t[NMAX] < 0.0)
            std += EXTERR_EPH * SQR(t[NMAX]) / 2.0;
        *vare = SQR(std);
    }
    /* linear interpolation for clock */
    t[0] = timediff(time, nav->peph[index].time);
    t[1] = timediff(time, nav->peph[index + 1].time);
    c[0] = nav->peph[index].pos[sat - 1][3];
    c[1] = nav->peph[index + 1].pos[sat - 1][3];

    if (t[0] <= 0.0) {
        if ((dts[0] = c[0]) != 0.0) {
            std = nav->peph[index].std[sat - 1][3] * CLIGHT - EXTERR_CLK * t[0];
        }
    } else if (t[1] >= 0.0) {
        if ((dts[0] = c[1]) != 0.0) {
            std = nav->peph[index + 1].std[sat - 1][3] * CLIGHT + EXTERR_CLK * t[1];
        }
    } else if (c[0] != 0.0 && c[1] != 0.0) {
        dts[0] = (c[1] * t[0] - c[0] * t[1]) / (t[0] - t[1]);
        i      = t[0] < -t[1] ? 0 : 1;
        std    = nav->peph[index + i].std[sat - 1][3] + EXTERR_CLK * fabs(t[i]);
    } else {
        dts[0] = 0.0;
    }
    if (varc)
        *varc = SQR(std);
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts, double *varc) {
    double t[2], c[2], std;
    int i, j, k, index;

    trace(4, "pephclk : time=%s sat=%2d\n", time_str(time, 3), sat);

    if (nav->nc < 2 || timediff(time, nav->pclk[0].time) < -MAXDTE ||
        timediff(time, nav->pclk[nav->nc - 1].time) > MAXDTE) {
        trace(3, "no prec clock %s sat=%2d\n", time_str(time, 0), sat);
        return 1;
    }
    /* binary search */
    for (i = 0, j = nav->nc - 1; i < j;) {
        k = (i + j) / 2;
        if (timediff(nav->pclk[k].time, time) < 0.0)
            i = k + 1;
        else
            j = k;
    }
    index = i <= 0 ? 0 : i - 1;

    /* linear interpolation for clock */
    t[0] = timediff(time, nav->pclk[index].time);
    t[1] = timediff(time, nav->pclk[index + 1].time);
    c[0] = nav->pclk[index].clk[sat - 1][0];
    c[1] = nav->pclk[index + 1].clk[sat - 1][0];

    if (t[0] <= 0.0) {
        if ((dts[0] = c[0]) == 0.0)
            return 0;
        std = nav->pclk[index].std[sat - 1][0] * CLIGHT - EXTERR_CLK * t[0];
    } else if (t[1] >= 0.0) {
        if ((dts[0] = c[1]) == 0.0)
            return 0;
        std = nav->pclk[index + 1].std[sat - 1][0] * CLIGHT + EXTERR_CLK * t[1];
    } else if (c[0] != 0.0 && c[1] != 0.0) {
        dts[0] = (c[1] * t[0] - c[0] * t[1]) / (t[0] - t[1]);
        i      = t[0] < -t[1] ? 0 : 1;
        std    = nav->pclk[index + i].std[sat - 1][0] * CLIGHT + EXTERR_CLK * fabs(t[i]);
    } else {
        trace(3, "prec clock outage %s sat=%2d\n", time_str(time, 0), sat);
        return 0;
    }
    if (varc)
        *varc = SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
 * compute satellite antenna phase center offset in ecef
 * args   : gtime_t time       I   time (gpst)
 *          double *rs         I   satellite position and velocity (ecef)
 *                                 {x,y,z,vx,vy,vz} (m|m/s)
 *          int    sat         I   satellite number
 *          nav_t  *nav        I   navigation data
 *          int    opt         I   output mode (0: UC, 1: IF12)
 *          double *dant       O   satellite antenna phase center offset (ecef)
 *                                 {dx,dy,dz} (m) (iono-free LC value)
 * return : none
 * notes  : iono-free LC frequencies defined as follows:
 *            GPS/QZSS : L1-L2
 *            GLONASS  : G1-G2
 *            Galileo  : E1-E5b
 *            BDS      : B1I-B2I
 *            NavIC    : L5-S
 * 
 *          fixed: add UC mode [dx1 dy1 dz1 dx2 dy2 dz2 ...]
 *-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,int opt, double *dant) {
    const pcv_t *pcv = nav->pcvs + sat - 1;
    double ex[3], ey[3], ez[3], es[3], r[3], rsun[3], gmst, erpv[5] = {0}, freq[2];
    double C1, C2, dant1, dant2;
    int i,j, sys;

    trace(4, "satantoff: time=%s sat=%2d\n", time_str(time, 3), sat);

    dant[0] = dant[1] = dant[2] = 0.0;
    if (opt == 0) { // UC
        for (i = 1; i < NFREQ; ++i) {
            dant[i * 3 + 0] = dant[i * 3 + 1] = dant[i * 3 + 2] = 0.0;
        }
    }

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, &gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i = 0; i < 3; i++)
        r[i] = -rs[i];
    if (!normv3(r, ez))
        return;
    for (i = 0; i < 3; i++)
        r[i] = rsun[i] - rs[i];
    if (!normv3(r, es))
        return;
    cross3(ez, es, r);
    if (!normv3(r, ey))
        return;
    cross3(ey, ez, ex);

    if (opt == 0) { // UC
        for (i = 0; i < NFREQ; ++i) { // frequency index
            for (j = 0; j < 3; ++j) { // x y z component
                dant[i * 3 + j] = pcv->off[i][0] * ex[j] + pcv->off[i][1] * ey[j] + pcv->off[i][2] * ez[j];
            }
        }
        return;
    }

    /* iono-free LC coefficients */
    sys = satsys(sat, NULL);
    freq[0] = idx2freq(sys,0);
    freq[1] = idx2freq(sys,1);
    if(freq[0]<=0 || freq[1]<=0){
        return ;
    }

    C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
    C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

    /* iono-free LC */
    for (i = 0; i < 3; i++) {
        dant1   = pcv->off[0][0] * ex[i] + pcv->off[0][1] * ey[i] + pcv->off[0][2] * ez[i];
        dant2   = pcv->off[1][0] * ex[i] + pcv->off[1][1] * ey[i] + pcv->off[1][2] * ez[i];
        dant[i] = C1 * dant1 + C2 * dant2;
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
 * compute satellite position/clock with precise ephemeris/clock
 * args   : gtime_t time       I   time (gpst)
 *          int    sat         I   satellite number
 *          nav_t  *nav        I   navigation data
 *          int    opt         I   sat postion option
 *                                 (0: center of mass, 1: antenna phase center)
 *          double *rs         O   sat position and velocity (ecef)
 *                                 {x,y,z,vx,vy,vz} (m|m/s)
 *          double *dts        O   sat clock {bias,drift} (s|s/s)
 *          double *var        IO  sat position and clock error variance (m)
 *                                 (NULL: no output)
 * return : status (1:ok,0:error or data outage)
 * notes  : clock includes relativistic correction but does not contain code bias
 *          before calling the function, nav->peph, nav->ne, nav->pclk and
 *          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
 *          if precise clocks are not set, clocks in sp3 are used instead
 * 
 *          fixed: opt does nothing now
 *-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt, double *rs, double *dts, double *var) {
    gtime_t time_tt;
    double rss[3], rst[3], dtss[1], dtst[1], dant[3] = {0}, vare = 0.0, varc = 0.0, tt = 1E-3;
    int i;

    trace(4, "peph2pos: time=%s sat=%2d opt=%d\n", time_str(time, 3), sat, opt);

    if (sat <= 0 || MAXSAT < sat)
        return 0;

    /* satellite position and clock bias */
    if (!pephpos(time, sat, nav, rss, dtss, &vare, &varc) || !pephclk(time, sat, nav, dtss, &varc))
        return 0;

    time_tt = timeadd(time, tt);
    if (!pephpos(time_tt, sat, nav, rst, dtst, NULL, NULL) || !pephclk(time_tt, sat, nav, dtst, NULL))
        return 0;

    /* satellite antenna offset correction */
    if (opt) {
        // fixed: moved to ppp_res
        // satantoff(time, rss, sat, nav, 1, dant); // fixme pco(IF) to pco(UC)
    }
    for (i = 0; i < 3; i++) {
        rs[i]     = rss[i] + dant[i];
        rs[i + 3] = (rst[i] - rss[i]) / tt;
    }
    /* relativistic effect correction */
    if (dtss[0] != 0.0) {
        dts[0] = dtss[0] - 2.0 * dot(rs, rs + 3, 3) / CLIGHT / CLIGHT;
        dts[1] = (dtst[0] - dtss[0]) / tt;
    } else { /* no precise clock */
        dts[0] = dts[1] = 0.0;
    }
    if (var)
        *var = vare + varc;

    return 1;
}
