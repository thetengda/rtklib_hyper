/*------------------------------------------------------------------------------
 * rtklib.h : RTKLIB constants, types and function prototypes
 *
 *          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
 *
 * options : -DENAGLO   enable GLONASS
 *           -DENAGAL   enable Galileo
 *           -DENAQZS   enable QZSS
 *           -DENACMP   enable BeiDou
 *           -DENAIRN   enable IRNSS
 *           -DNFREQ=n  set number of obs codes/frequencies
 *           -DNEXOBS=n set number of extended obs codes
 *           -DMAXOBS=n set max number of obs data in an epoch
 *           -DWIN32    use WIN32 API
 *           -DWIN_DLL  generate library as Windows DLL
 *
 * version : $Revision:$ $Date:$
 * history : 2007/01/13 1.0  rtklib ver.1.0.0
 *           2007/03/20 1.1  rtklib ver.1.1.0
 *           2008/07/15 1.2  rtklib ver.2.1.0
 *           2008/10/19 1.3  rtklib ver.2.1.1
 *           2009/01/31 1.4  rtklib ver.2.2.0
 *           2009/04/30 1.5  rtklib ver.2.2.1
 *           2009/07/30 1.6  rtklib ver.2.2.2
 *           2009/12/25 1.7  rtklib ver.2.3.0
 *           2010/07/29 1.8  rtklib ver.2.4.0
 *           2011/05/27 1.9  rtklib ver.2.4.1
 *           2013/03/28 1.10 rtklib ver.2.4.2
 *           2020/11/30 1.11 rtklib ver.2.4.3 b34
 *-----------------------------------------------------------------------------*/
#ifndef RTKLIB_H
#define RTKLIB_H
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef WIN32
#include <windows.h>
#include <winsock2.h>
#else
#include <pthread.h>
#include <sys/select.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT extern
#endif

/* constants -----------------------------------------------------------------*/

#define VER_RTKLIB "2.4.3" /* library version */

#define PATCH_LEVEL "b34" /* patch level */

#define COPYRIGHT_RTKLIB "Copyright (C) 2007-2020 T.Takasu\nAll rights reserved."

#define SQR(x) ((x) * (x))
#define SQRT(x) ((x) <= 0.0 || (x) != (x) ? 0.0 : sqrt(x))
#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define ROUND(x) ((int)floor((x) + 0.5))
#define SGN(x) ((x) <= 0.0 ? -1.0 : 1.0)

#define PI 3.1415926535897932  /* pi */
#define D2R (PI / 180.0)       /* deg to rad */
#define R2D (180.0 / PI)       /* rad to deg */
#define CLIGHT 299792458.0     /* speed of light (m/s) */
#define SC2RAD 3.1415926535898 /* semi-circle to radian (IS-GPS) */
#define AU 149597870691.0      /* 1 AU (m) */
#define AS2R (D2R / 3600.0)    /* arc sec to radian */

#define OMGE 7.2921151467E-5 /* earth angular velocity (IS-GPS) (rad/s) */
#define RE_WGS84 6378137.0             /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84 (1.0 / 298.257223563) /* earth flattening (WGS84) */
#define HION 350000.0 /* ionosphere height (m) */

#define MAXFREQ 7 /* max NFREQ */

#define FREQ1 1.57542E9       /* L1/E1/B1C  frequency (Hz) */
#define FREQ2 1.22760E9       /* L2         frequency (Hz) */
#define FREQ5 1.17645E9       /* L5/E5a/B2a frequency (Hz) */
#define FREQ6 1.27875E9       /* E6/L6  frequency (Hz) */
#define FREQ7 1.20714E9       /* E5b    frequency (Hz) */
#define FREQ8 1.191795E9      /* E5a+b  frequency (Hz) */
#define FREQ9 2.492028E9      /* S      frequency (Hz) */
#define FREQ1_GLO 1.60200E9   /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO 0.56250E6   /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO 1.24600E9   /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO 0.43750E6   /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO 1.202025E9  /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO 1.600995E9 /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO 1.248060E9 /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP 1.561098E9  /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP 1.20714E9   /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP 1.26852E9   /* BDS B3      frequency (Hz) */

#define EFACT_GPS 1.0 /* error factor: GPS */
#define EFACT_GLO 1.5 /* error factor: GLONASS */
#define EFACT_GAL 1.0 /* error factor: Galileo */
#define EFACT_QZS 1.0 /* error factor: QZSS */
#define EFACT_SBS 3.0 /* error factor: SBAS */
#define EFACT_CMP 1.0 /* error factor: BeiDou */
#define EFACT_IRN 1.5 /* error factor: IRNSS */
#define EFACT_LEO 1.0 /* error factor: LEO */

#define SYS_NONE 0x00 /* navigation system: none */
#define SYS_GPS 0x01  /* navigation system: GPS */
#define SYS_GLO 0x04  /* navigation system: GLONASS */
#define SYS_GAL 0x08  /* navigation system: Galileo */
#define SYS_QZS 0x10  /* navigation system: QZSS */
#define SYS_SBS 0x02  /* navigation system: SBAS */
#define SYS_CMP 0x20  /* navigation system: BeiDou */
#define SYS_IRN 0x40  /* navigation system: IRNS */
#define SYS_LEO 0x80  /* navigation system: LEO */
#define SYS_ALL 0xFF  /* navigation system: all */

#ifndef NFREQ
#define NFREQ 3 /* number of carrier frequencies */
#endif
#define NFREQGLO 2 /* number of carrier frequencies of GLONASS */

#ifndef MAXNPCV
#define MAXNPCV 100 /* max number of pcv in one line */
#endif

#ifndef NEXOBS
#define NEXOBS 0 /* number of extended obs codes */
#endif

#define SNR_UNIT 0.001 /* SNR unit (dBHz) */

#define MINPRNGPS 1                         /* min satellite PRN number of GPS */
#define MAXPRNGPS 32                        /* max satellite PRN number of GPS */
#define NSATGPS (MAXPRNGPS - MINPRNGPS + 1) /* number of GPS satellites */
#define NSYSGPS 1

#ifdef ENAGLO
#define MINPRNGLO 1                         /* min satellite slot number of GLONASS */
#define MAXPRNGLO 27                        /* max satellite slot number of GLONASS */
#define NSATGLO (MAXPRNGLO - MINPRNGLO + 1) /* number of GLONASS satellites */
#define NSYSGLO 1
#else
#define MINPRNGLO 0
#define MAXPRNGLO 0
#define NSATGLO 0
#define NSYSGLO 0
#endif
#ifdef ENAGAL
#define MINPRNGAL 1                         /* min satellite PRN number of Galileo */
#define MAXPRNGAL 36                        /* max satellite PRN number of Galileo */
#define NSATGAL (MAXPRNGAL - MINPRNGAL + 1) /* number of Galileo satellites */
#define NSYSGAL 1
#else
#define MINPRNGAL 0
#define MAXPRNGAL 0
#define NSATGAL 0
#define NSYSGAL 0
#endif
#ifdef ENAQZS
#define MINPRNQZS 193                       /* min satellite PRN number of QZSS */
#define MAXPRNQZS 202                       /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                     /* min satellite PRN number of QZSS L1S */
#define MAXPRNQZS_S 191                     /* max satellite PRN number of QZSS L1S */
#define NSATQZS (MAXPRNQZS - MINPRNQZS + 1) /* number of QZSS satellites */
#define NSYSQZS 1
#else
#define MINPRNQZS 0
#define MAXPRNQZS 0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS 0
#define NSYSQZS 0
#endif
#ifdef ENACMP
#define MINPRNCMP 1                         /* min satellite sat number of BeiDou */
#define MAXPRNCMP 63                        /* max satellite sat number of BeiDou */
#define NSATCMP (MAXPRNCMP - MINPRNCMP + 1) /* number of BeiDou satellites */
#define NSYSCMP 1
#else
#define MINPRNCMP 0
#define MAXPRNCMP 0
#define NSATCMP 0
#define NSYSCMP 0
#endif
#ifdef ENAIRN
#define MINPRNIRN 1                         /* min satellite sat number of IRNSS */
#define MAXPRNIRN 14                        /* max satellite sat number of IRNSS */
#define NSATIRN (MAXPRNIRN - MINPRNIRN + 1) /* number of IRNSS satellites */
#define NSYSIRN 1
#else
#define MINPRNIRN 0
#define MAXPRNIRN 0
#define NSATIRN 0
#define NSYSIRN 0
#endif
#ifdef ENALEO
#define MINPRNLEO 1                         /* min satellite sat number of LEO */
#define MAXPRNLEO 10                        /* max satellite sat number of LEO */
#define NSATLEO (MAXPRNLEO - MINPRNLEO + 1) /* number of LEO satellites */
#define NSYSLEO 1
#else
#define MINPRNLEO 0
#define MAXPRNLEO 0
#define NSATLEO 0
#define NSYSLEO 0
#endif
// todo add SBS switch
#define MINPRNSBS 120                       /* min satellite PRN number of SBAS */
#define MAXPRNSBS 158                       /* max satellite PRN number of SBAS */
#define NSATSBS (MAXPRNSBS - MINPRNSBS + 1) /* number of SBAS satellites */
#define NSYSSBS 1

#define NSYS (NSYSGPS + NSYSGLO + NSYSGAL + NSYSQZS + NSYSSBS + NSYSCMP + NSYSIRN + NSYSLEO) /* number of systems */

#define MAXSAT (NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN + NSATSBS + NSATLEO)
#define MAXSTA 255 /* max satellite number (1 to MAXSAT) */

#ifndef MAXOBS
#define MAXOBS 96 /* max number of obs in an epoch */
#endif
#define MAXRCV 64     /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE 64 /* max number of obs type in RINEX */
#ifdef OBS_100HZ
#define DTTOL 0.005 /* tolerance of time difference (s) */
#else
#define DTTOL 0.025 /* tolerance of time difference (s) */
#endif
#define MAXDTOE 7200.0      /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0  /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0 /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0 /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0  /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0  /* max time difference to IRNSS Toe (s) */
#define MAXDTOE_SBS 360.0   /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S 86400.0   /* max time difference to ephem toe (s) for other */

#define INT_SWAP_TRAC 86400.0 /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0 /* swap interval of solution status file (s) */

#define MAXEXFILE 1024    /* max number of expanded files */
#define MAXCOMMENT 100    /* max number of RINEX comments */
#define MAXSTRPATH 1024   /* max length of stream path */
#define MAXSOLMSG 8191    /* max length of solution message */
#define MAXERRMSG 4096    /* max length of error/warning message */
#define MAXANT 64         /* max length of station name/antenna type */
#define MAXLEAPS 64       /* max number of leap seconds table */

#define RNX2VER 2.10 /* RINEX ver.2 default output version */
#define RNX3VER 3.00 /* RINEX ver.3 default output version */

enum tsys_e {
    TSYS_GPS = 0, /* time system: GPS time */
    TSYS_UTC = 1, /* time system: UTC */
    TSYS_GLO = 2, /* time system: GLONASS time */
    TSYS_GAL = 3, /* time system: Galileo time */
    TSYS_QZS = 4, /* time system: QZSS time */
    TSYS_CMP = 5, /* time system: BeiDou time */
    TSYS_IRN = 6, /* time system: IRNSS time */
};

// must keep same order with ObsCodes
enum code_e {
    CODE_NONE = 0,  /* obs code: none or unknown */
    CODE_L1C  = 1,  /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
    CODE_L1P  = 2,  /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
    CODE_L1W  = 3,  /* obs code: L1 Z-track (GPS) */
    CODE_L1Y  = 4,  /* obs code: L1Y        (GPS) */
    CODE_L1M  = 5,  /* obs code: L1M        (GPS) */
    CODE_L1N  = 6,  /* obs code: L1codeless,B1codeless (GPS,BDS) */
    CODE_L1S  = 7,  /* obs code: L1C(D)     (GPS,QZS) */
    CODE_L1L  = 8,  /* obs code: L1C(P)     (GPS,QZS) */
    CODE_L1E  = 9,  /* (not used) */
    CODE_L1A  = 10, /* obs code: E1A,B1A    (GAL,BDS) */
    CODE_L1B  = 11, /* obs code: E1B        (GAL) */
    CODE_L1X  = 12, /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
    CODE_L1Z  = 13, /* obs code: E1A+B+C,L1S (GAL,QZS) */
    CODE_L2C  = 14, /* obs code: L2C/A,G1C/A (GPS,GLO) */
    CODE_L2D  = 15, /* obs code: L2 L1C/A-(P2-P1) (GPS) */
    CODE_L2S  = 16, /* obs code: L2C(M)     (GPS,QZS) */
    CODE_L2L  = 17, /* obs code: L2C(L)     (GPS,QZS) */
    CODE_L2X  = 18, /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
    CODE_L2P  = 19, /* obs code: L2P,G2P    (GPS,GLO) */
    CODE_L2W  = 20, /* obs code: L2 Z-track (GPS) */
    CODE_L2Y  = 21, /* obs code: L2Y        (GPS) */
    CODE_L2M  = 22, /* obs code: L2M        (GPS) */
    CODE_L2N  = 23, /* obs code: L2codeless (GPS) */
    CODE_L5I  = 24, /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
    CODE_L5Q  = 25, /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
    CODE_L5X  = 26, /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
    CODE_L7I  = 27, /* obs code: E5bI,B2bI  (GAL,BDS) */
    CODE_L7Q  = 28, /* obs code: E5bQ,B2bQ  (GAL,BDS) */
    CODE_L7X  = 29, /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
    CODE_L6A  = 30, /* obs code: E6A,B3A    (GAL,BDS) */
    CODE_L6B  = 31, /* obs code: E6B        (GAL) */
    CODE_L6C  = 32, /* obs code: E6C        (GAL) */
    CODE_L6X  = 33, /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
    CODE_L6Z  = 34, /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
    CODE_L6S  = 35, /* obs code: L6S        (QZS) */
    CODE_L6L  = 36, /* obs code: L6L        (QZS) */
    CODE_L8I  = 37, /* obs code: E5abI      (GAL) */
    CODE_L8Q  = 38, /* obs code: E5abQ      (GAL) */
    CODE_L8X  = 39, /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
    CODE_L2I  = 40, /* obs code: B1_2I      (BDS) */
    CODE_L2Q  = 41, /* obs code: B1_2Q      (BDS) */
    CODE_L6I  = 42, /* obs code: B3I        (BDS) */
    CODE_L6Q  = 43, /* obs code: B3Q        (BDS) */
    CODE_L3I  = 44, /* obs code: G3I        (GLO) */
    CODE_L3Q  = 45, /* obs code: G3Q        (GLO) */
    CODE_L3X  = 46, /* obs code: G3I+Q      (GLO) */
    CODE_L1I  = 47, /* obs code: B1I        (BDS) (obsolute) */
    CODE_L1Q  = 48, /* obs code: B1Q        (BDS) (obsolute) */
    CODE_L5A  = 49, /* obs code: L5A SPS    (IRN) */
    CODE_L5B  = 50, /* obs code: L5B RS(D)  (IRN) */
    CODE_L5C  = 51, /* obs code: L5C RS(P)  (IRN) */
    CODE_L9A  = 52, /* obs code: SA SPS     (IRN) */
    CODE_L9B  = 53, /* obs code: SB RS(D)   (IRN) */
    CODE_L9C  = 54, /* obs code: SC RS(P)   (IRN) */
    CODE_L9X  = 55, /* obs code: SB+C       (IRN) */
    CODE_L1D  = 56, /* obs code: B1D        (BDS) */
    CODE_L5D  = 57, /* obs code: L5D(L5S),B2aD (QZS,BDS) */
    CODE_L5P  = 58, /* obs code: L5P(L5S),B2aP (QZS,BDS) */
    CODE_L5Z  = 59, /* obs code: L5D+P(L5S) (QZS) */
    CODE_L6E  = 60, /* obs code: L6E        (QZS) */
    CODE_L7D  = 61, /* obs code: B2bD       (BDS) */
    CODE_L7P  = 62, /* obs code: B2bP       (BDS) */
    CODE_L7Z  = 63, /* obs code: B2bD+P     (BDS) */
    CODE_L8D  = 64, /* obs code: B2abD      (BDS) */
    CODE_L8P  = 65, /* obs code: B2abP      (BDS) */
    CODE_L4A  = 66, /* obs code: G1aL1OCd   (GLO) */
    CODE_L4B  = 67, /* obs code: G1aL1OCd   (GLO) */
    CODE_L4X  = 68, /* obs code: G1al1OCd+p (GLO) */
};
#define MAXCODE 68  /* max number of obs code */

enum pmode_e {
    PMODE_SINGLE     = 0, /* positioning mode: single */
    PMODE_DGPS       = 1, /* positioning mode: DGPS/DGNSS */
    PMODE_KINEMA     = 2, /* positioning mode: kinematic */
    PMODE_STATIC     = 3, /* positioning mode: static */
    PMODE_MOVEB      = 4, /* positioning mode: moving-base */
    PMODE_FIXED      = 5, /* positioning mode: fixed */
    PMODE_PPP_KINEMA = 6, /* positioning mode: PPP-kinemaric */
    PMODE_PPP_STATIC = 7, /* positioning mode: PPP-static */
    PMODE_PPP_FIXED  = 8, /* positioning mode: PPP-fixed */
};

enum solfmt_e {
    SOLF_LLH = 0, /* solution format: lat/lon/height */
    SOLF_XYZ = 1, /* solution format: x/y/z-ecef */
    SOLF_ENU = 2, /* solution format: e/n/u-baseline */
};

enum solstat_e {
    SOLQ_NONE   = 0, /* solution status: no solution */
    SOLQ_FIX    = 1, /* solution status: fix */
    SOLQ_FLOAT  = 2, /* solution status: float */
    SOLQ_DGPS   = 4, /* solution status: DGPS/DGNSS */
    SOLQ_SINGLE = 5, /* solution status: single */
    SOLQ_PPP    = 6, /* solution status: PPP */
    SOLQ_DR     = 7, /* solution status: dead reconing */
    MAXSOLQ     = 7, /* max number of solution status */
};

enum ionoopt_e {
    IONOOPT_OFF  = 0, /* ionosphere option: correction off */
    IONOOPT_BRDC = 1, /* ionosphere option: broadcast model */
    IONOOPT_IFLC = 3, /* ionosphere option: L1/L2 iono-free LC */
    IONOOPT_EST  = 4, /* ionosphere option: estimation */
    IONOOPT_TEC  = 5, /* ionosphere option: IONEX TEC model */
};
enum tropopt_e {
    TROPOPT_OFF  = 0, /* troposphere option: correction off */
    TROPOPT_SAAS = 1, /* troposphere option: Saastamoinen model */
    TROPOPT_EST  = 3, /* troposphere option: ZTD estimation */
    TROPOPT_ESTG = 4, /* troposphere option: ZTD+grad estimation */
};

enum ephopt_e {
    EPHOPT_BRDC = 0, /* ephemeris option: broadcast ephemeris */
    EPHOPT_PREC = 1, /* ephemeris option: precise ephemeris */
};

enum armode_e {
    ARMODE_OFF     = 0, /* AR mode: off */
    ARMODE_CONT    = 1, /* AR mode: continuous */
    ARMODE_INST    = 2, /* AR mode: instantaneous */
    ARMODE_FIXHOLD = 3, /* AR mode: fix and hold */
    // ARMODE_WLNL    = 4, /* AR mode: wide lane/narrow lane */
    // ARMODE_TCAR    = 5, /* AR mode: triple carrier ar */
};
enum armethod_e {
    ARMETHOD_ROUND=0,
    ARMETHOD_ILS=1,
};

enum posopt_e {
    POSOPT_POS    = 0, /* pos option: LLH/XYZ */
    POSOPT_SINGLE = 1, /* pos option: average of single pos */
    POSOPT_FILE   = 2, /* pos option: read from pos file */
    POSOPT_RINEX  = 3, /* pos option: rinex header pos */
};
enum geoid_e {
    GEOID_EMBEDDED    = 0, /* geoid model: embedded geoid */
    GEOID_EGM96_M150  = 1, /* geoid model: EGM96 15x15" */
    GEOID_EGM2008_M25 = 2, /* geoid model: EGM2008 2.5x2.5" */
    GEOID_EGM2008_M10 = 3, /* geoid model: EGM2008 1.0x1.0" */
    GEOID_GSI2000_M15 = 4, /* geoid model: GSI geoid 2000 1.0x1.5" */
    GEOID_RAF09       = 5, /* geoid model: IGN RAF09 for France 1.5"x2" */
};

#define COMMENTH "%"                   /* comment line indicator for solution */
#define MSG_DISCONN "$_DISCONNECT\r\n" /* disconnect message */

#define LLI_SLIP 0x01   /* LLI: cycle-slip */
#define LLI_HALFC 0x02  /* LLI: half-cycle not resovled */
#define LLI_BOCTRK 0x04 /* LLI: boc tracking of mboc signal */
#define LLI_HALFA 0x40  /* LLI: half-cycle added */
#define LLI_HALFS 0x80  /* LLI: half-cycle subtracted */

#ifdef WIN32
#define thread_t HANDLE
#define lock_t CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f) EnterCriticalSection(f)
#define unlock(f) LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#else
#define thread_t pthread_t
#define lock_t pthread_mutex_t
#define initlock(f) pthread_mutex_init(f, NULL)
#define lock(f) pthread_mutex_lock(f)
#define unlock(f) pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
#endif


enum orbtype_e{
    ORB_LEO, 
    ORB_MEO,
    ORB_IGSO,
    ORB_GEO
};

/* type definitions ----------------------------------------------------------*/

typedef struct { /* time struct */
    time_t time; /* time (s) expressed by standard time_t */
    double sec;  /* fraction of second under 1 s */
} gtime_t;

typedef struct {                  /* observation data record */
    gtime_t time;                 /* receiver sampling time (GPST) */
    uint8_t sat, rcv;             /* satellite/receiver number */
    uint16_t SNR[NFREQ + NEXOBS]; /* signal strength (0.001 dBHz) */
    uint8_t LLI[NFREQ + NEXOBS];  /* loss of lock indicator */
    enum code_e code[NFREQ + NEXOBS]; /* code indicator (CODE_???) */
    double L[NFREQ + NEXOBS];     /* observation data carrier-phase (cycle) */
    double P[NFREQ + NEXOBS];     /* observation data pseudorange (m) */
    float D[NFREQ + NEXOBS];      /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {  /* observation data */
    int n, nmax;  /* number of obervation data/allocated */
    obsd_t *data; /* observation data records */
} obs_t;

typedef struct {     /* earth rotation parameter data type */
    double mjd;      /* mjd (days) */
    double xp, yp;   /* pole offset (rad) */
    double xpr, ypr; /* pole offset rate (rad/day) */
    double ut1_utc;  /* ut1-utc (s) */
    double lod;      /* length of day (s/day) */
} erpd_t;

typedef struct {  /* earth rotation parameter type */
    int n, nmax;  /* number and max number of data */
    erpd_t *data; /* earth rotation parameter data */
} erp_t;

typedef struct {                /* antenna parameter type */
    int sat;                    /* satellite number (0:receiver) */
    char type[MAXANT];          /* antenna type */
    char code[MAXANT];          /* serial number or satellite code */
    gtime_t ts, te;             /* valid time start and end */
    double off[NFREQ][3];       /* phase center offset e/n/u or x/y/z (m), index by [ifreq][ixyz] */
    double var[NFREQ][MAXNPCV]; /* phase center variation (m) */
    double zen1, zen2, dzen;    /* column grid for pcv */
    double dazi;                /* row grid for azi-dependent(0~360) pcv*/
} pcv_t;

typedef struct { /* antenna parameters type */
    int n, nmax; /* number of data/allocated */
    pcv_t *pcv;  /* antenna parameters data */
} pcvs_t;

typedef struct {           /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;               /* satellite number */
    int iode, iodc;        /* IODE,IODC */
    int sva;               /* SV accuracy (URA index) */
    int svh;               /* SV health (0:ok) */
    int week;              /* GPS/QZS: gps week, GAL: galileo week */
    int code;              /* GPS/QZS: code on L2 */
                           /* GAL: data source defined as rinex 3.03 */
                           /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
    int flag;              /* GPS/QZS: L2 P data flag */
                           /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toe, toc, ttr; /* Toe,Toc,T_trans */
                           /* SV orbit parameters */
    double A, e, i0, OMG0, omg, M0, deln, OMGd, idot;
    double crc, crs, cuc, cus, cic, cis;
    double toes;       /* Toe (s) in week */
    double fit;        /* fit interval (h) */
    double f0, f1, f2; /* SV clock parameters (af0,af1,af2) */
    double tgd[6];     /* group delay parameters */
                       /* GPS/QZS:tgd[0]=TGD */
                       /* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
                       /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
                       /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
    double Adot, ndot; /* Adot,ndot for CNAV */
} eph_t;

typedef struct {       /* GLONASS broadcast ephemeris type */
    int sat;           /* satellite number */
    int iode;          /* IODE (0-6 bit of tb field) */
    int frq;           /* satellite frequency number */
    int svh, sva, age; /* satellite health, accuracy, age of operation */
    gtime_t toe;       /* epoch of epherides (gpst) */
    gtime_t tof;       /* message frame time (gpst) */
    double pos[3];     /* satellite position (ecef) (m) */
    double vel[3];     /* satellite velocity (ecef) (m/s) */
    double acc[3];     /* satellite acceleration (ecef) (m/s^2) */
    double taun, gamn; /* SV clock bias (s)/relative freq bias */
    double dtaun;      /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {           /* precise ephemeris type */
    gtime_t time;          /* time (GPST) */
    int index;             /* ephemeris index for multiple files */
    double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float std[MAXSAT][4];  /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float vst[MAXSAT][4];  /* satellite velocity/clk-rate std (m/s|s/s) */
    float cov[MAXSAT][3];  /* satellite position covariance (m^2) */
    float vco[MAXSAT][3];  /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {           /* precise clock type */
    gtime_t time;          /* time (GPST) */
    int index;             /* clock index for multiple files */
    double clk[MAXSAT][1]; /* satellite clock (s) */
    float std[MAXSAT][1];  /* satellite clock std (s) */
} pclk_t;

typedef struct {     /* SBAS ephemeris type */
    int sat;         /* satellite number */
    gtime_t t0;      /* reference epoch time (GPST) */
    gtime_t tof;     /* time of message frame (GPST) */
    int sva;         /* SV accuracy (URA index) */
    int svh;         /* SV health (0:ok) */
    double pos[3];   /* satellite position (m) (ecef) */
    double vel[3];   /* satellite velocity (m/s) (ecef) */
    double acc[3];   /* satellite acceleration (m/s^2) (ecef) */
    double af0, af1; /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {    /* TEC grid type */
    gtime_t time;   /* epoch time (GPST) */
    int ndata[3];   /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;      /* earth radius (km) */
    double lats[3]; /* latitude start/interval (deg) */
    double lons[3]; /* longitude start/interval (deg) */
    double hgts[3]; /* heights start/interval (km) */
    double *data;   /* TEC grid data (tecu) */
    float *rms;     /* RMS values (tecu) */
} tec_t;

typedef struct {                  /* navigation data type */
    int n, nmax;                  /* number of broadcast ephemeris */
    int ng, ngmax;                /* number of glonass ephemeris */
    int ns, nsmax;                /* number of sbas ephemeris */
    int ne, nemax;                /* number of precise ephemeris */
    int nc, ncmax;                /* number of precise clock */
    int na, namax;                /* number of almanac data */
    int nt, ntmax;                /* number of tec grid data */
    eph_t *eph;                   /* GPS/QZS/GAL/BDS/IRN ephemeris */
    geph_t *geph;                 /* GLONASS ephemeris */
    seph_t *seph;                 /* SBAS ephemeris */
    peph_t *peph;                 /* precise ephemeris */
    pclk_t *pclk;                 /* precise clock */
    tec_t *tec;                   /* tec grid data */
    erp_t erp;                    /* earth rotation parameters */
    // double utc_gps[8];            /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
    // double utc_glo[8];            /* GLONASS UTC time parameters {tau_C,tau_GPS} */
    // double utc_gal[8];            /* Galileo UTC parameters */
    // double utc_qzs[8];            /* QZS UTC parameters */
    // double utc_cmp[8];            /* BeiDou UTC parameters */
    // double utc_irn[9];            /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
    // double utc_sbs[4];            /* SBAS UTC parameters */
    double ion_gps[8];            /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];            /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];            /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];            /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];            /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int glo_fcn[32];              /* GLONASS FCN + 8 */
    // double cbias[MAXSAT][3];      /* satellite DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    // double rbias[MAXRCV][2][3];   /* receiver DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    pcv_t pcvs[MAXSAT];           /* satellite antenna pcv */

    double satosb[MAXSAT][2][MAXCODE];
    double rcvosb[NSYS][2][MAXCODE]; // GREJCIL
} nav_t;

typedef struct {           /* station parameter type */
    char name[MAXANT];     /* marker name */
    char marker[MAXANT];   /* marker number */
    char antdes[MAXANT];   /* antenna descriptor */
    char antsno[MAXANT];   /* antenna serial number */
    char rectype[MAXANT];  /* receiver type descriptor */
    char recver[MAXANT];   /* receiver firmware version */
    char recsno[MAXANT];   /* receiver serial number */
    int antsetup;          /* antenna setup id */
    int itrf;              /* ITRF realization year */
    int deltype;           /* antenna delta type (0:enu,1:xyz) */
    double pos[3];         /* station position (ecef) (m) */
    double del[3];         /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;            /* antenna height (m) */
    int glo_cp_align;      /* GLONASS code-phase alignment (0:no,1:yes) */
    double glo_cp_bias[4]; /* GLONASS code-phase biases {1C,1P,2C,2P} (m) */
} sta_t;

typedef struct {   /* solution type */
    gtime_t time;  /* time (GPST) */
    double rr[6];  /* position/velocity (m|m/s) */
                   /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    float qr[6];   /* position variance/covariance (m^2) */
                   /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                   /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    float qv[6];   /* velocity variance/covariance (m^2/s^2) */
    double dtr[NSYS+1]; /* receiver clock bias to time systems (s) */
    uint8_t type;  /* type (0:xyz-ecef,1:enu-baseline) */
    uint8_t stat;  /* solution status (SOLQ_???) */
    uint8_t ns;    /* number of valid satellites */
    float age;     /* age of differential (s) */
    float ratio;   /* AR ratio factor for valiation */
    float thres;   /* AR ratio threshold for valiation */
} sol_t;

typedef struct {                 /* solution buffer type */
    int n, nmax;                 /* number of solution/max number of buffer */
    int cyclic;                  /* cyclic buffer flag */
    int start, end;              /* start/end index */
    gtime_t time;                /* current solution time */
    sol_t *data;                 /* solution data */
    double rb[3];                /* reference position {x,y,z} (ecef) (m) */
    uint8_t buff[MAXSOLMSG + 1]; /* message buffer */
    int nb;                      /* number of byte in message buffer */
} solbuf_t;

typedef struct {    /* solution status type */
    gtime_t time;   /* time (GPST) */
    uint8_t sat;    /* satellite number */
    uint8_t frq;    /* frequency (1:L1,2:L2,...) */
    float az, el;   /* azimuth/elevation angle (rad) */
    float resp;     /* pseudorange residual (m) */
    float resc;     /* carrier-phase residual (m) */
    uint8_t flag;   /* flags: (vsat<<5)+(slip<<3)+fix */
    uint16_t snr;   /* signal strength (*SNR_UNIT dBHz) */
    uint16_t lock;  /* lock counter */
    uint16_t outc;  /* outage counter */
    uint16_t slipc; /* slip counter */
    uint16_t rejc;  /* reject counter */
} solstat_t;

typedef struct {     /* solution status buffer type */
    int n, nmax;     /* number of solution/max number of buffer */
    solstat_t *data; /* solution status data */
} solstatbuf_t;

typedef struct {         /* option type */
    const char *name;    /* option name */
    int format;          /* option format (0:int,1:double,2:string,3:enum) */
    void *var;           /* pointer to option variable */
    const char *comment; /* option comment/enum labels/unit */
} opt_t;

typedef struct {           /* SNR mask type */
    int ena[2];            /* enable flag {rover,base} */
    double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {             /* processing options type */
    int mode;                /* positioning mode (PMODE_???) */
    int soltype;             /* solution type (0:forward,1:backward,2:combined) */
    int nf;                  /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
    int navsys;              /* navigation system */
    double elmin;            /* elevation mask angle (rad) */
    snrmask_t snrmask;       /* SNR mask */
    int sateph;              /* satellite ephemeris/clock (EPHOPT_???) */
    int modear;              /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
    int glomodear;           /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int galmodear;           /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int bdsmodear;           /* BeiDou AR mode (0:off,1:on) */
    int methodar;            /* AR method 0:round,1:ILS*/
    int maxout;              /* obs outage count to reset bias */
    int minlock;             /* min lock count to fix ambiguity */
    int minfix;              /* min fix count to hold ambiguity */
    int armaxiter;           /* max iteration to resolve ambiguity */
    int ionoopt;             /* ionosphere option (IONOOPT_???) */
    int tropopt;             /* troposphere option (TROPOPT_???) */
    int dynamics;            /* dynamics model (0:none,1:velociy,2:accel) */
    int tidecorr;            /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
    int niter;               /* number of filter iteration */
    int codesmooth;          /* code smoothing window size (0:none) */
    int intpref;             /* interpolate reference obs (for post mission) */
    int rovpos;              /* rover position for fixed mode */
    int refpos;              /* base position for relative mode */
                             /* (0:pos in prcopt,  1:average of single pos, */
                             /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    double eratio[NFREQ];    /* code/phase error ratio */
    double err[5];           /* measurement error factor */
                             /* [0]:reserved */
                             /* [1-3]:error factor a/b/c of phase (m) */
                             /* [4]:doppler frequency (hz) */
    double std[10];          /* initial-state std [0]pos,[1]vel [2]acc [3]clk [4]trop [5]tropg [6]iono [7]dcb [8]bias [9]glo_ifb */
    double prn[9];           /* process-noise std [0]pos,[1]acch [2]accv [3]clk [4]trop [5]tropg [6]iono [7]dcb [8]bias */
    double sclkstab;         /* satellite clock stability (sec/sec) */
    double thresar[8];       /* AR validation threshold: ratio, confidence, fraction,  */
    double elmaskar;         /* elevation mask of AR for rising satellite (deg) */
    double elmaskhold;       /* elevation mask to hold ambiguity (deg) */
    double thresslip[2];     /* slip threshold of geometry-free phase/MW (m) */
    double maxtdiff;         /* max difference of time (sec) */
    double maxinno;          /* reject threshold of innovation (m) */
    double maxgdop;          /* reject threshold of gdop */
    double baseline[2];      /* baseline length constraint {const,sigma} (m) */
    double ru[3];            /* rover position for fixed mode {x,y,z} (ecef) (m) */
    double rb[3];            /* base position for relative mode {x,y,z} (ecef) (m) */
    char anttype[2][MAXANT]; /* antenna types {rover,base} */
    double antdel[2][3];     /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    pcv_t pcvr[2];           /* receiver antenna parameters {rov,base} */
    uint8_t exsats[MAXSAT];  /* excluded satellites (1:excluded,2:included) */
    int outsingle;           /* output single by dgps/float/fix/ppp outage */
    char rnxopt[2][256];     /* rinex options {rover,base} */
    int posopt[6];           /* positioning options */
    double odisp[2][6 * 11]; /* ocean tide loading parameters {rov,base} */
    double maxgapiono;          /* obs outage count to reset iono */
} prcopt_t;

typedef struct {        /* solution options type */
    int posf;           /* solution format (SOLF_???) */
    int times;          /* time system (TIMES_???) */
    int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;          /* time digits under decimal point */
    int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;        /* output header (0:no,1:yes) */
    int outopt;         /* output processing options (0:no,1:yes) */
    int outvel;         /* output velocity options (0:no,1:yes) */
    int datum;          /* datum (0:WGS84,1:Tokyo) */
    int height;         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;      /* solution of static mode (0:all,1:single) */
    int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;          /* debug trace level (0:off,1-5:debug) */
    char sep[64];       /* field separator */
    char prog[64];      /* program name */
    double maxsolstd;   /* max std-dev for solution output (m) (0:all) */
} solopt_t;

typedef struct {              /* file options type */
    char satantp[MAXSTRPATH]; /* satellite antenna parameters file */
    char rcvantp[MAXSTRPATH]; /* receiver antenna parameters file */
    char stapos[MAXSTRPATH];  /* station positions file */
    char geoid[MAXSTRPATH];   /* external geoid data file */
    char iono[MAXSTRPATH];    /* ionosphere data file */
    char osb[MAXSTRPATH];     /* dcb data file */
    char eop[MAXSTRPATH];     /* eop data file */
    char blq[MAXSTRPATH];     /* ocean tide loading blq file */
    char solstat[MAXSTRPATH]; /* solution statistics file */
    char trace[MAXSTRPATH];   /* debug trace file */
} filopt_t;

typedef struct {                  /* RINEX options type */
    gtime_t ts, te;               /* time start/end */
    double tint;                  /* time interval (s) */
    double ttol;                  /* time tolerance (s) */
    double tunit;                 /* time unit for multiple-session (s) */
    int rnxver;                   /* RINEX version (x100) */
    int navsys;                   /* navigation system */
    int obstype;                  /* observation type */
    int freqtype;                 /* frequency type */
    char mask[7][64];             /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid[32];               /* station id for rinex file name */
    char prog[32];                /* program */
    char runby[32];               /* run-by */
    char marker[64];              /* marker name */
    char markerno[32];            /* marker number */
    char markertype[32];          /* marker type (ver.3) */
    char name[2][32];             /* observer/agency */
    char rec[3][32];              /* receiver #/type/vers */
    char ant[3][32];              /* antenna #/type */
    double apppos[3];             /* approx position x/y/z */
    double antdel[3];             /* antenna delta h/e/n */
    double glo_cp_bias[4];        /* GLONASS code-phase biases (m) */
    char comment[MAXCOMMENT][64]; /* comments */
    char rcvopt[256];             /* receiver dependent options */
    uint8_t exsats[MAXSAT];       /* excluded satellites */
    int glofcn[32];               /* glonass fcn+8 */
    int outiono;                  /* output iono correction */
    int outtime;                  /* output time system correction */
    int outleaps;                 /* output leap seconds */
    int autopos;                  /* auto approx position */
    int phshift;                  /* phase shift correction */
    int halfcyc;                  /* half cycle correction */
    int sep_nav;                  /* separated nav files */
    gtime_t tstart;               /* first obs time */
    gtime_t tend;                 /* last obs time */
    gtime_t trtcm;                /* approx log start time for rtcm */
    char tobs[7][MAXOBSTYPE][4];  /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    double shift[7][MAXOBSTYPE];  /* phase shift (cyc) {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];                  /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {           /* satellite status type */
    uint8_t sys;           /* navigation system */
    uint8_t vs;            /* valid satellite flag single */ // fixed: spp will not influence ppp
    double azel[2];        /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ];    /* residuals of pseudorange (m) */
    double resc[NFREQ];    /* residuals of carrier-phase (m) */
    uint8_t vsat[NFREQ];   /* valid satellite flag */
    uint16_t snr[NFREQ];   /* signal strength (*SNR_UNIT dBHz) */
    uint8_t fix[NFREQ];    /* ambiguity fix flag (1:fix,2:float,3:hold) */
    uint8_t slip[NFREQ];   /* cycle-slip flag */
    uint8_t half[NFREQ];   /* half-cycle valid flag */
    int lock[NFREQ];       /* lock counter of phase */
    uint32_t outc[NFREQ];  /* obs outage counter of phase */
    uint32_t slipc[NFREQ]; /* cycle-slip counter */
    uint32_t rejc[NFREQ];  /* reject counter */
    double gf[NFREQ - 1];  /* geometry-free phase (m) */
    double mw[NFREQ - 1];  /* MW-LC (m) */
    double phw;            /* phase windup (cycle) */
    gtime_t pt[2][NFREQ];  /* previous carrier-phase time */
    double ph[2][NFREQ];   /* previous carrier-phase observable (cycle) */
} ssat_t;

typedef struct {        /* ambiguity control type */
    gtime_t epoch[4];   /* last epoch */
    int n[4];           /* number of epochs */
    int LCf[4];           /* linear combination frequency */
    double LC[4];       /* linear combination average */
    double LCv[4];      /* linear combination variance */
    int fixcnt;         /* fix count */
    char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {            /* RTK control/result type */
    sol_t sol;              /* RTK solution */
    double rb[6];           /* base position/velocity (ecef) (m|m/s) */
    int nx, na;             /* number of float states/fixed states */
    double tt;              /* time difference between current and previous (s) */
    double *x, *P;          /* float states and their covariance */
    double *xa, *Pa;        /* fixed states and their covariance */
    int nfix;               /* number of continuous fixes of ambiguity */
    ambc_t ambc[MAXSAT];    /* ambibuity control */
    ssat_t ssat[MAXSAT];    /* satellite status */
    int neb;                /* bytes in error message buffer */
    char errbuf[MAXERRMSG]; /* error message buffer */
    prcopt_t opt;           /* processing options */
} rtk_t;

typedef void fatalfunc_t(const char *); /* fatal callback function type */

/* global variables ----------------------------------------------------------*/
extern const prcopt_t PrcoptDefault;     /* default positioning options */
extern const solopt_t SoloptDefault;     /* default solution output options */
extern opt_t OptionsTable[];                   /* system options table */

extern const char SysCodes[NSYS+1];
extern const int NavSyss[NSYS+1];
extern const double FreqVals[NSYS][MAXFREQ];

/* satellites, systems, codes functions --------------------------------------*/
EXPORT int satno(int sys, int prn);
EXPORT int satsys(int sat, int *prn);

//todo 
EXPORT int bds2or3(int sat);
EXPORT int sat2orbtype(int sat);

EXPORT int code2sys(char code);
EXPORT int sys2idx(int sys);

EXPORT int satid2no(const char *id);
EXPORT void satno2id(int sat, char *id);
EXPORT uint8_t obs2code(const char *obs);
EXPORT char *code2obs(uint8_t code);
EXPORT double code2freq(int sys, uint8_t code, int fcn);
EXPORT double sat2freq(int sat, uint8_t code, const nav_t *nav);
EXPORT int code2idx(int sys, uint8_t code, double *freq);
EXPORT int cfreq2idx(int sys, const char *cfreq, double *freq);
EXPORT double idx2freq(int sys,int ifreq);
EXPORT int satexclude(int sat, double var, int svh, const prcopt_t *opt);
EXPORT int testsnr(int base, int freq, double el, double snr, const snrmask_t *mask);
EXPORT void setcodepri(int sys, int idx, const char *pri);
EXPORT int getcodepri(int sys, uint8_t code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
EXPORT double *mat(int n, int m);
EXPORT int *imat(int n, int m);
EXPORT double *zeros(int n, int m);
EXPORT double *eye(int n);
EXPORT double dot(const double *a, const double *b, int n);
EXPORT double norm(const double *a, int n);
EXPORT void cross3(const double *a, const double *b, double *c);
EXPORT int normv3(const double *a, double *b);
EXPORT void matcpy(double *A, const double *B, int n, int m);
EXPORT void matmul(const char *tr, int n, int k, int m, double alpha, const double *A, const double *B, double beta,
                   double *C);
EXPORT int matinv(double *A, int n);
EXPORT int solve(const char *tr, const double *A, const double *Y, int n, int m, double *X);
EXPORT int lsq(const double *A, const double *y, int n, int m, double *x, double *Q);
EXPORT int filter(double *x, double *P, const double *H, const double *v, const double *R, int n, int m);
EXPORT int smoother(const double *xf, const double *Qf, const double *xb, const double *Qb, int n, double *xs,
                    double *Qs);
EXPORT void matprint(const double *A, int n, int m, int p, int q);
EXPORT void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);

EXPORT void add_fatal(fatalfunc_t *func);

/* time and string functions -------------------------------------------------*/
EXPORT double str2num(const char *s, int i, int n);
EXPORT int str2time(const char *s, int i, int n, gtime_t *t);
EXPORT void time2str(gtime_t t, char *str, int n);
EXPORT gtime_t epoch2time(const double *ep);
EXPORT void time2epoch(gtime_t t, double *ep);
EXPORT gtime_t gpst2time(int week, double sec);
EXPORT double time2gpst(gtime_t t, int *week);
EXPORT gtime_t gst2time(int week, double sec);
EXPORT double time2gst(gtime_t t, int *week);
EXPORT gtime_t bdt2time(int week, double sec);
EXPORT double time2bdt(gtime_t t, int *week);
EXPORT char *time_str(gtime_t t, int n);

EXPORT gtime_t timeadd(gtime_t t, double sec);
EXPORT double timediff(gtime_t t1, gtime_t t2);
EXPORT gtime_t gpst2utc(gtime_t t);
EXPORT gtime_t utc2gpst(gtime_t t);
EXPORT gtime_t gpst2bdt(gtime_t t);
EXPORT gtime_t bdt2gpst(gtime_t t);
EXPORT gtime_t timeget(void);
EXPORT void timeset(gtime_t t);
EXPORT void timereset(void);
EXPORT double time2doy(gtime_t t);
EXPORT double utc2gmst(gtime_t t, double ut1_utc);
EXPORT int read_leaps(const char *file);

EXPORT int adjgpsweek(int week);
EXPORT uint32_t tickget(void);
EXPORT void sleepms(int ms);

EXPORT int reppath(const char *path, char *rpath, gtime_t time, const char *rov, const char *base);
EXPORT int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts, gtime_t te, const char *rov,
                    const char *base);

/* coordinates transformation ------------------------------------------------*/
EXPORT void ecef2pos(const double *r, double *pos);
EXPORT void pos2ecef(const double *pos, double *r);
EXPORT void ecef2enu(const double *pos, const double *r, double *e);
EXPORT void enu2ecef(const double *pos, const double *e, double *r);
EXPORT void covenu(const double *pos, const double *P, double *Q);
EXPORT void covecef(const double *pos, const double *Q, double *P);
EXPORT void xyz2enu(const double *pos, double *E);
EXPORT void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
EXPORT void deg2dms(double deg, double *dms, int ndec);
EXPORT double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
EXPORT void readpos(const char *file, const char *rcv, double *pos);
EXPORT int sortobs(obs_t *obs);
EXPORT void uniqnav(nav_t *nav);
EXPORT int screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
EXPORT void freeobs(obs_t *obs);
EXPORT void freenav(nav_t *nav, int opt);
EXPORT int readblq(const char *file, const char *sta, double *odisp);
EXPORT int readerp(const char *file, erp_t *erp);
EXPORT int geterp(const erp_t *erp, gtime_t time, double *val);

/* debug trace functions -----------------------------------------------------*/
EXPORT void traceopen(const char *file);
EXPORT void traceclose(void);
EXPORT void tracelevel(int level);
EXPORT void trace(int level, const char *format, ...);
EXPORT void tracet(int level, const char *format, ...);
EXPORT void tracemat(int level, const double *A, int n, int m, int p, int q);
EXPORT void traceobs(int level, const obsd_t *obs, int n);
EXPORT void tracenav(int level, const nav_t *nav);
EXPORT void tracegnav(int level, const nav_t *nav);
EXPORT void tracehnav(int level, const nav_t *nav);
EXPORT void tracepeph(int level, const nav_t *nav);
EXPORT void tracepclk(int level, const nav_t *nav);
EXPORT void traceb(int level, const uint8_t *p, int n);

/* platform dependent functions ----------------------------------------------*/
EXPORT int execcmd(const char *cmd);
EXPORT int expath(const char *path, char *paths[], int nmax);
EXPORT void createdir(const char *path);

/* positioning models --------------------------------------------------------*/
EXPORT double satazel(const double *pos, const double *e, double *azel);
EXPORT double geodist(const double *rs, const double *rr, double *e);
EXPORT void dops(int ns, const double *azel, double elmin, double *dop);

/* atmosphere models ---------------------------------------------------------*/
EXPORT double ionmodel(gtime_t t, const double *ion, const double *pos, const double *azel);
EXPORT double ionmapf(const double *pos, const double *azel);
EXPORT double ionppp(const double *pos, const double *azel, double re, double hion, double *pppos);
EXPORT double tropmodel(gtime_t time, const double *pos, const double *azel, double humi);
EXPORT double tropmapf(gtime_t time, const double *pos, const double *azel, double *mapfw);
EXPORT int iontec(gtime_t time, const nav_t *nav, const double *pos, const double *azel, int opt, double *delay,
                  double *var);
EXPORT void readtec(const char *file, nav_t *nav, int opt);
EXPORT int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos, const double *azel, int ionoopt,
                    double *ion, double *var);
EXPORT int tropcorr(gtime_t time, const nav_t *nav, const double *pos, const double *azel, int tropopt, double *trp,
                    double *var);

/* antenna models ------------------------------------------------------------*/
EXPORT int readpcv(const char *file, pcvs_t *pcvs);
EXPORT pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs);
EXPORT void antmodel(const pcv_t *pcv, const double *del, const double *azel, int opt, double *dant);
EXPORT void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
EXPORT void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun, double *rmoon, double *gmst);
EXPORT void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp, const double *odisp, double *dr);

/* geiod models --------------------------------------------------------------*/
EXPORT int opengeoid(int model, const char *file);
EXPORT void closegeoid(void);
EXPORT double geoidh(const double *pos);

/* rinex functions -----------------------------------------------------------*/
EXPORT int readrnx(const char *file, int rcv, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
EXPORT int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te, double tint, const char *opt, obs_t *obs,
                    nav_t *nav, sta_t *sta);
EXPORT int readrnxc(const char *file, nav_t *nav);
EXPORT int rtk_uncompress(const char *file, char *uncfile);

/* ephemeris and clock functions ---------------------------------------------*/
EXPORT double eph2clk(gtime_t time, const eph_t *eph);
EXPORT double geph2clk(gtime_t time, const geph_t *geph);
EXPORT double seph2clk(gtime_t time, const seph_t *seph);
EXPORT void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts, double *var);
EXPORT void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts, double *var);
EXPORT void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts, double *var);
EXPORT int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt, double *rs, double *dts, double *var);
EXPORT void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,int opt, double *dant);
EXPORT int satpos(gtime_t time, gtime_t teph, int sat, int ephopt, const nav_t *nav, double *rs, double *dts,
                  double *var, int *svh);
EXPORT void satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav, int sateph, double *rs, double *dts,
                    double *var, int *svh);
EXPORT void setseleph(int sys, int sel);
EXPORT int getseleph(int sys);
EXPORT void readsp3(const char *file, nav_t *nav, int opt);
EXPORT int readsap(const char *file, gtime_t time, nav_t *nav);
EXPORT int readosb(const char *file, nav_t *nav, const sta_t *sta);

/* receiver raw data functions -----------------------------------------------*/
EXPORT uint32_t getbitu(const uint8_t *buff, int pos, int len);
EXPORT void setbitu(uint8_t *buff, int pos, int len, uint32_t data);

EXPORT int outprcopts(uint8_t *buff, const prcopt_t *opt);
EXPORT int outsolheads(uint8_t *buff, const solopt_t *opt);
EXPORT int outsols(uint8_t *buff, const sol_t *sol, const double *rb, const solopt_t *opt);
EXPORT void outprcopt(FILE *fp, const prcopt_t *opt);
EXPORT void outsolhead(FILE *fp, const solopt_t *opt);
EXPORT void outsol(FILE *fp, const sol_t *sol, const double *rb, const solopt_t *opt);

/* options functions ---------------------------------------------------------*/
EXPORT opt_t *searchopt(const char *name, const opt_t *opts);
EXPORT int str2opt(opt_t *opt, const char *str);
EXPORT int loadopts(const char *file, opt_t *opts);
EXPORT void resetsysopts(void);
EXPORT void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);

/* integer ambiguity resolution ----------------------------------------------*/
EXPORT int lambda(int n, int m, const double *a, const double *Q, double *F, double *s);
EXPORT int lambda_reduction(int n, const double *Q, double *Z);
EXPORT int lambda_search(int n, int m, const double *a, const double *Q, double *F, double *s);

/* standard positioning ------------------------------------------------------*/
EXPORT int pntpos(const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel,
                  ssat_t *ssat, char *msg);

/* precise positioning -------------------------------------------------------*/
EXPORT void rtkinit(rtk_t *rtk, const prcopt_t *opt);
EXPORT void rtkfree(rtk_t *rtk);
EXPORT int rtkpos(rtk_t *rtk, const obsd_t *obs, int nobs, const nav_t *nav);
EXPORT int rtkopenstat(const char *file, int level);
EXPORT void rtkclosestat(void);
EXPORT int rtkoutstat(rtk_t *rtk, char *buff);

/* precise point positioning -------------------------------------------------*/
EXPORT void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
EXPORT int pppnx(const prcopt_t *opt);
EXPORT int pppoutstat(rtk_t *rtk, char *buff);

/* post-processing positioning -----------------------------------------------*/
EXPORT int postpos(gtime_t ts, gtime_t te, double ti, double tu, const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile, const char *rov, const char *base);

/* application defined functions ---------------------------------------------*/
extern int showmsg(const char *format, ...);
extern void settspan(gtime_t ts, gtime_t te);
extern void settime(gtime_t time);

#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */
