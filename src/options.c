/*------------------------------------------------------------------------------
 * options.c : options functions
 *
 *          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
 *
 * version : $Revision:$ $Date:$
 * history : 2010/07/20  1.1  moved from postpos.c
 *                            added api:
 *                                searchopt(),str2opt(),opt2str(),opt2buf(),
 *                                loadopts(),saveopts(),resetsysopts(),
 *                                getsysopts(),setsysopts()
 *           2010/09/11  1.2  add options
 *                                pos2-elmaskhold,pos1->snrmaskena
 *                                pos1-snrmask1,2,3
 *           2013/03/11  1.3  add pos1-posopt1,2,3,4,5,pos2-syncsol
 *                                misc-rnxopt1,2,pos1-snrmask_r,_b,_L1,_L2,_L5
 *           2014/10/21  1.4  add pos2-bdsarmode
 *           2015/02/20  1.4  add ppp-fixed as pos1-posmode option
 *           2015/05/10  1.5  add pos2-arthres1,2,3,4
 *           2015/05/31  1.6  add pos2-armaxiter, pos1-posopt6
 *                            add selection precise for pos1-pospot3
 *           2015/11/26  1.7  modify pos1-frequency 4:l1+l2+l5+l6 -> l1+l5
 *           2015/12/05  1.8  add misc-pppopt
 *           2016/06/10  1.9  add ant2-maxaveep,ant2-initrst
 *           2016/07/31  1.10 add out-outsingle,out-maxsolstd
 *           2017/06/14  1.11 add out-outvel
 *           2020/11/30  1.12 change options pos1-frequency, pos1-ionoopt,
 *                             pos1-tropopt, pos1-sateph, pos1-navsys,
 *                             pos2-gloarmode,
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* system options buffer -----------------------------------------------------*/
static prcopt_t PrcOpt;
static solopt_t SolOpt;
static filopt_t FilOpt;
static int AntPosType[2];
static double ElMask, ElMaskAr, ElMaskHold;
static double AntPos[2][3];
static char ExSats[1024];
static char SnrMask[NFREQ][1024];

/* system options table ------------------------------------------------------*/
#define SWTOPT "0:off,1:on"
#define MODOPT "0:single,1:dgps,2:kinematic,3:static,4:movingbase,5:fixed,6:ppp-kine,7:ppp-static,8:ppp-fixed"
#define FRQOPT "1:l1,2:l1+2,3:l1+2+3,4:l1+2+3+4,5:l1+2+3+4+5"
#define TYPOPT "0:forward,1:backward,2:combined"
#define IONOPT "0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc"
#define TRPOPT "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad"
#define EPHOPT "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define NAVOPT "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:bds+64:navic"
#define SOLOPT "0:llh,1:xyz,2:enu,3:nmea"
#define TSYOPT "0:gpst,1:utc,2:jst"
#define TFTOPT "0:tow,1:hms"
#define DFTOPT "0:deg,1:dms"
#define HGTOPT "0:ellipsoidal,1:geodetic"
#define GEOOPT "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STAOPT "0:all,1:single"
#define STSOPT "0:off,1:state,2:residual"

#define ARMOPT "0:off,1:continuous,2:instantaneous,3:fix-and-hold"
#define GLOAROPT "0:off,1:on"
#define GALAROPT "0:off,1:on"
#define BDSAROPT "0:off,1:bds2,2:bds3,3:bds2+3"
#define ARMETHODOPT "0:round,1:ils"

#define POSOPT "0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm,6:raw"
#define TIDEOPT "0:off,1:on,2:otl"
#define PHWOPT "0:off,1:on,2:precise"

// label, type(0:int,1:double,2:string,3:enum), variable, description
opt_t OptionsTable[] = {
    {"pos1-posmode", 3, (void *) &PrcOpt.mode, MODOPT},
    {"pos1-frequency", 3, (void *) &PrcOpt.nf, FRQOPT},
    {"pos1-soltype", 3, (void *) &PrcOpt.soltype, TYPOPT},
    {"pos1-elmask", 1, (void *) &ElMask, "deg"},
    {"pos1-snrmask_r", 3, (void *) &PrcOpt.snrmask.ena[0], SWTOPT},
    {"pos1-snrmask_b", 3, (void *) &PrcOpt.snrmask.ena[1], SWTOPT},
    {"pos1-snrmask_L1", 2, (void *) SnrMask[0], ""},
    {"pos1-snrmask_L2", 2, (void *) SnrMask[1], ""},
    {"pos1-snrmask_L5", 2, (void *) SnrMask[2], ""},
    {"pos1-dynamics", 3, (void *) &PrcOpt.dynamics, SWTOPT},
    {"pos1-tidecorr", 3, (void *) &PrcOpt.tidecorr, TIDEOPT},
    {"pos1-ionoopt", 3, (void *) &PrcOpt.ionoopt, IONOPT},
    {"pos1-tropopt", 3, (void *) &PrcOpt.tropopt, TRPOPT},
    {"pos1-sateph", 3, (void *) &PrcOpt.sateph, EPHOPT},
    {"pos1-satpcv", 3, (void *) &PrcOpt.posopt[0], SWTOPT},
    {"pos1-rcvpcv", 3, (void *) &PrcOpt.posopt[1], SWTOPT},
    {"pos1-phw", 3, (void *) &PrcOpt.posopt[2], PHWOPT},
    {"pos1-eclipse", 3, (void *) &PrcOpt.posopt[3], SWTOPT},
    {"pos1-raim", 3, (void *) &PrcOpt.posopt[4], SWTOPT},
    {"pos1-dayclkjump", 3, (void *) &PrcOpt.posopt[5], SWTOPT},
    {"pos1-exclsats", 2, (void *) ExSats, "prn ..."},
    {"pos1-navsys", 0, (void *) &PrcOpt.navsys, NAVOPT},

    {"pos2-armode", 3, (void *) &PrcOpt.modear, ARMOPT},
    {"pos2-gloarmode", 3, (void *) &PrcOpt.glomodear, GLOAROPT},
    {"pos2-galarmode", 3, (void *) &PrcOpt.galmodear, GALAROPT},
    {"pos2-bdsarmode", 3, (void *) &PrcOpt.bdsmodear, BDSAROPT},
    {"pos2-armethod", 3, (void *) &PrcOpt.methodar, ARMETHODOPT},

    {"pos2-arthres", 1, (void *) &PrcOpt.thresar[0], ""},
    {"pos2-arthres1", 1, (void *) &PrcOpt.thresar[1], ""},
    {"pos2-arthres2", 1, (void *) &PrcOpt.thresar[2], ""},
    {"pos2-arthres3", 1, (void *) &PrcOpt.thresar[3], ""},
    {"pos2-arthres4", 1, (void *) &PrcOpt.thresar[4], ""},
    {"pos2-arlockcnt", 0, (void *) &PrcOpt.minlock, ""},
    {"pos2-arelmask", 1, (void *) &ElMaskAr, "deg"},
    {"pos2-arminfix", 0, (void *) &PrcOpt.minfix, ""},
    {"pos2-armaxiter", 0, (void *) &PrcOpt.armaxiter, ""},
    {"pos2-elmaskhold", 1, (void *) &ElMaskHold, "deg"},
    {"pos2-aroutcnt", 0, (void *) &PrcOpt.maxout, ""},
    {"pos2-maxage", 1, (void *) &PrcOpt.maxtdiff, "s"},
    {"pos2-slipthres_GF", 1, (void *) &PrcOpt.thresslip[0], "m"},
    {"pos2-slipthres_MW", 1, (void *) &PrcOpt.thresslip[1], "m"},
    {"pos2-rejionno", 1, (void *) &PrcOpt.maxinno, "m"},
    {"pos2-rejgdop", 1, (void *) &PrcOpt.maxgdop, ""},
    {"pos2-niter", 0, (void *) &PrcOpt.niter, ""},
    {"pos2-baselen", 1, (void *) &PrcOpt.baseline[0], "m"},
    {"pos2-basesig", 1, (void *) &PrcOpt.baseline[1], "m"},
    {"pos2-maxgapiono", 1, (void *) &PrcOpt.maxgapiono, " max outage count to reset iono"},

    {"out-solformat", 3, (void *) &SolOpt.posf, SOLOPT},
    {"out-outhead", 3, (void *) &SolOpt.outhead, SWTOPT},
    {"out-outopt", 3, (void *) &SolOpt.outopt, SWTOPT},
    {"out-outvel", 3, (void *) &SolOpt.outvel, SWTOPT},
    {"out-timesys", 3, (void *) &SolOpt.times, TSYOPT},
    {"out-timeform", 3, (void *) &SolOpt.timef, TFTOPT},
    {"out-timendec", 0, (void *) &SolOpt.timeu, ""},
    {"out-degform", 3, (void *) &SolOpt.degf, DFTOPT},
    {"out-fieldsep", 2, (void *) SolOpt.sep, ""},
    {"out-outsingle", 3, (void *) &PrcOpt.outsingle, SWTOPT},
    {"out-maxsolstd", 1, (void *) &SolOpt.maxsolstd, "m"},
    {"out-height", 3, (void *) &SolOpt.height, HGTOPT},
    {"out-geoid", 3, (void *) &SolOpt.geoid, GEOOPT},
    {"out-solstatic", 3, (void *) &SolOpt.solstatic, STAOPT},
    {"out-outstat", 3, (void *) &SolOpt.sstat, STSOPT},

    {"stats-eratio1", 1, (void *) &PrcOpt.eratio[0], ""},
    {"stats-eratio2", 1, (void *) &PrcOpt.eratio[1], ""},
    {"stats-errphase", 1, (void *) &PrcOpt.err[1], "m"},
    {"stats-errphaseel", 1, (void *) &PrcOpt.err[2], "m"},
    {"stats-errphasebl", 1, (void *) &PrcOpt.err[3], "m/10km"},
    {"stats-errdoppler", 1, (void *) &PrcOpt.err[4], "Hz"},

    {"stats-stdpos", 1, (void *) &PrcOpt.std[0], "m"},
    {"stats-stdvel", 1, (void *) &PrcOpt.std[1], "m"},
    {"stats-stdacc", 1, (void *) &PrcOpt.std[2], "m"},
    {"stats-stdclk", 1, (void *) &PrcOpt.std[3], "m"},
    {"stats-stdtrop", 1, (void *) &PrcOpt.std[4], "m"},
    {"stats-stdtropg", 1, (void *) &PrcOpt.std[5], "m"},
    {"stats-stdiono", 1, (void *) &PrcOpt.std[6], "m"},
    {"stats-stddcb", 1, (void *) &PrcOpt.std[7], "m"},
    {"stats-stdbias", 1, (void *) &PrcOpt.std[8], "m"},
    {"stats-stdglo_ifb", 1, (void *) &PrcOpt.std[9], "m"},

    {"stats-prnpos", 1, (void *)    &PrcOpt.prn[0], "m/s^0.5"},
    {"stats-prnaccelh", 1, (void *) &PrcOpt.prn[1], "m/s^2"},
    {"stats-prnaccelv", 1, (void *) &PrcOpt.prn[2], "m/s^2"},
    {"stats-prnclk", 1, (void *)    &PrcOpt.prn[3], "m/s"},
    {"stats-prntrop", 1, (void *)   &PrcOpt.prn[4], "m/s"},
    {"stats-prntropg", 1, (void *)  &PrcOpt.prn[5], "m/s"},
    {"stats-prniono", 1, (void *)   &PrcOpt.prn[6], "m/s"},
    {"stats-prndcb", 1, (void *)    &PrcOpt.prn[7], "m/s"},
    {"stats-prnbias", 1, (void *)   &PrcOpt.prn[8], "m/s"},

    {"stats-clkstab", 1, (void *) &PrcOpt.sclkstab, "s/s"},

    {"ant1-postype", 3, (void *) &AntPosType[0], POSOPT},
    {"ant1-pos1", 1, (void *) &AntPos[0][0], "deg|m"},
    {"ant1-pos2", 1, (void *) &AntPos[0][1], "deg|m"},
    {"ant1-pos3", 1, (void *) &AntPos[0][2], "m|m"},
    {"ant1-anttype", 2, (void *) PrcOpt.anttype[0], ""},
    {"ant1-antdele", 1, (void *) &PrcOpt.antdel[0][0], "m"},
    {"ant1-antdeln", 1, (void *) &PrcOpt.antdel[0][1], "m"},
    {"ant1-antdelu", 1, (void *) &PrcOpt.antdel[0][2], "m"},

    {"ant2-postype", 3, (void *) &AntPosType[1], POSOPT},
    {"ant2-pos1", 1, (void *) &AntPos[1][0], "deg|m"},
    {"ant2-pos2", 1, (void *) &AntPos[1][1], "deg|m"},
    {"ant2-pos3", 1, (void *) &AntPos[1][2], "m|m"},
    {"ant2-anttype", 2, (void *) PrcOpt.anttype[1], ""},
    {"ant2-antdele", 1, (void *) &PrcOpt.antdel[1][0], "m"},
    {"ant2-antdeln", 1, (void *) &PrcOpt.antdel[1][1], "m"},
    {"ant2-antdelu", 1, (void *) &PrcOpt.antdel[1][2], "m"},

    {"misc-timeinterp", 3, (void *) &PrcOpt.intpref, SWTOPT},
    {"misc-rnxopt1", 2, (void *) PrcOpt.rnxopt[0], ""},
    {"misc-rnxopt2", 2, (void *) PrcOpt.rnxopt[1], ""},

    {"file-satantfile", 2, (void *) &FilOpt.satantp, ""},
    {"file-rcvantfile", 2, (void *) &FilOpt.rcvantp, ""},
    {"file-staposfile", 2, (void *) &FilOpt.stapos, ""},
    {"file-geoidfile", 2, (void *) &FilOpt.geoid, ""},
    {"file-ionofile", 2, (void *) &FilOpt.iono, ""},
    {"file-blqfile", 2, (void *) &FilOpt.blq, ""},
    {"file-solstatfile", 2, (void *) &FilOpt.solstat, ""},
    {"file-tracefile", 2, (void *) &FilOpt.trace, ""},

    {"", 0, NULL, ""} /* terminator */
};
/* discard space characters at tail ------------------------------------------*/
static void chop(char *str) {
    char *p;
    if ((p = strchr(str, '#')))
        *p = '\0'; /* comment */
    for (p = str + strlen(str) - 1; p >= str && !isgraph((int) *p); p--)
        *p = '\0';
}
/* string to enum ------------------------------------------------------------*/
static int str2enum(const char *str, const char *comment, int *val) {
    const char *p;
    char s[32];

    for (p = comment;; p++) {
        if (!(p = strstr(p, str)))
            break;
        if (*(p - 1) != ':')
            continue;
        for (p -= 2; '0' <= *p && *p <= '9'; p--)
            ;
        return sscanf(p + 1, "%d", val) == 1;
    }
    sprintf(s, "%.30s:", str);
    if ((p = strstr(comment, s))) { /* number */
        return sscanf(p, "%d", val) == 1;
    }
    return 0;
}
/* search option ---------------------------------------------------------------
 * search option record
 * args   : char   *name     I  option name
 *          opt_t  *opts     I  options table
 *                              (terminated with table[i].name="")
 * return : option record (NULL: not found)
 *-----------------------------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts) {
    int i;

    trace(3, "searchopt: name=%s\n", name);

    for (i = 0; *opts[i].name; i++) {
        if (strstr(opts[i].name, name))
            return (opt_t *) (opts + i);
    }
    return NULL;
}
/* string to option value ------------------------------------------------------
 * convert string to option value
 * args   : opt_t  *opt      O  option
 *          char   *str      I  option value string
 * return : status (1:ok,0:error)
 *-----------------------------------------------------------------------------*/
extern int str2opt(opt_t *opt, const char *str) {
    switch (opt->format) {
        case 0:
            *(int *) opt->var = atoi(str);
            break;
        case 1:
            *(double *) opt->var = atof(str);
            break;
        case 2:
            strcpy((char *) opt->var, str);
            break;
        case 3:
            return str2enum(str, opt->comment, (int *) opt->var);
        default:
            return 0;
    }
    return 1;
}
/* load options ----------------------------------------------------------------
 * load options from file
 * args   : char   *file     I  options file
 *          opt_t  *opts     IO options table
 *                              (terminated with table[i].name="")
 * return : status (1:ok,0:error)
 *-----------------------------------------------------------------------------*/
extern int loadopts(const char *file, opt_t *opts) {
    FILE *fp;
    opt_t *opt;
    char buff[2048], *p;
    int n = 0;

    trace(3, "loadopts: file=%s\n", file);

    if (!(fp = fopen(file, "r"))) {
        trace(1, "loadopts: options file open error (%s)\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {
        n++;
        chop(buff);

        if (buff[0] == '\0')
            continue;

        if (!(p = strstr(buff, "="))) {
            fprintf(stderr, "invalid option %s (%s:%d)\n", buff, file, n);
            continue;
        }
        *p++ = '\0';
        chop(buff);
        if (!(opt = searchopt(buff, opts)))
            continue;

        if (!str2opt(opt, p)) {
            fprintf(stderr, "invalid option value %s (%s:%d)\n", buff, file, n);
            continue;
        }
    }
    fclose(fp);

    return 1;
}
/* system options buffer to options ------------------------------------------*/
static void buff2sysopts(void) {
    double pos[3], *rr;
    char buff[1024], *p, *id;
    int i, j, sat, *ps;

    PrcOpt.elmin      = ElMask * D2R;
    PrcOpt.elmaskar   = ElMaskAr * D2R;
    PrcOpt.elmaskhold = ElMaskHold * D2R;

    for (i = 0; i < 2; i++) {
        ps = i == 0 ? &PrcOpt.rovpos : &PrcOpt.refpos;
        rr = i == 0 ? PrcOpt.ru : PrcOpt.rb;

        if (AntPosType[i] == 0) { /* lat/lon/hgt */
            *ps    = 0;
            pos[0] = AntPos[i][0] * D2R;
            pos[1] = AntPos[i][1] * D2R;
            pos[2] = AntPos[i][2];
            pos2ecef(pos, rr);
        } else if (AntPosType[i] == 1) { /* xyz-ecef */
            *ps   = 0;
            rr[0] = AntPos[i][0];
            rr[1] = AntPos[i][1];
            rr[2] = AntPos[i][2];
        } else
            *ps = AntPosType[i] - 1;
    }
    /* excluded satellites */
    for (i = 0; i < MAXSAT; i++)
        PrcOpt.exsats[i] = 0;
    if (ExSats[0] != '\0') {
        strcpy(buff, ExSats);
        for (p = strtok(buff, " "); p; p = strtok(NULL, " ")) {
            if (*p == '+')
                id = p + 1;
            else
                id = p;
            if (!(sat = satid2no(id)))
                continue;
            PrcOpt.exsats[sat - 1] = *p == '+' ? 2 : 1;
        }
    }
    /* snrmask */
    for (i = 0; i < NFREQ; i++) {
        for (j = 0; j < 9; j++)
            PrcOpt.snrmask.mask[i][j] = 0.0;
        strcpy(buff, SnrMask[i]);
        for (p = strtok(buff, ","), j = 0; p && j < 9; p = strtok(NULL, ",")) {
            PrcOpt.snrmask.mask[i][j++] = atof(p);
        }
    }
    /* number of frequency (4:L1+L5) */
    // if (PrcOpt.nf == 4) {
    //     PrcOpt.nf      = 3;
    //     // PrcOpt.freqopt = 1;
    // }
}
/* reset system options to default ---------------------------------------------
 * reset system options to default
 * args   : none
 * return : none
 *-----------------------------------------------------------------------------*/
extern void resetsysopts(void) {
    int i, j;

    trace(3, "resetsysopts:\n");

    PrcOpt            = PrcoptDefault;
    SolOpt            = SoloptDefault;
    FilOpt.satantp[0] = '\0';
    FilOpt.rcvantp[0] = '\0';
    FilOpt.stapos[0]  = '\0';
    FilOpt.geoid[0]   = '\0';
    FilOpt.osb[0]     = '\0';
    FilOpt.blq[0]     = '\0';
    FilOpt.solstat[0] = '\0';
    FilOpt.trace[0]   = '\0';
    for (i = 0; i < 2; i++)
        AntPosType[i] = 0;
    ElMask     = 15.0;
    ElMaskAr   = 0.0;
    ElMaskHold = 0.0;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++) {
            AntPos[i][j] = 0.0;
        }
    ExSats[0] = '\0';
}
/* get system options ----------------------------------------------------------
 * get system options
 * args   : prcopt_t *popt   IO processing options (NULL: no output)
 *          solopt_t *sopt   IO solution options   (NULL: no output)
 *          folopt_t *fopt   IO file options       (NULL: no output)
 * return : none
 * notes  : to load system options, use loadopts() before calling the function
 *-----------------------------------------------------------------------------*/
extern void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt) {
    trace(3, "getsysopts:\n");

    buff2sysopts();
    if (popt)
        *popt = PrcOpt;
    if (sopt)
        *sopt = SolOpt;
    if (fopt)
        *fopt = FilOpt;
}