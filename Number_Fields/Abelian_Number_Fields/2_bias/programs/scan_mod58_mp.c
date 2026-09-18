/* scan_mod58_mp.c — everything needed for Q(zeta_5) (mod 5) and Q(zeta_8) (mod 8).
 * Same deterministic chunk/prefix framework as drh_scan_ak_mp.c (OpenMP inside a node, range
 * splittable across nodes, merged with merge_traj.py).
 *
 * Per prime p (u = p^{-1/2}, w = 1/p) the following sums are accumulated (Kahan):
 *   P5_1..P5_4  : sum_{p = a mod 5} u                     (class sums for Q(zeta_5); p=5 excluded)
 *   P8_1,3,5,7  : sum_{p = a mod 8} u                     (class sums for Q(zeta_8); p=2 excluded)
 *   ME          : sum 1/p            M2 : sum (1/p + log(1-1/p))        (all p)
 *   H4          : sum chi_5(p)/p    (chi_5 = Legendre mod 5; enters c(chi_2) of Q(zeta_5))
 *   K3_5, K3_m4, K3_m8, K3_8       : sum_{k>=3} sum_p chi(p^k)/(k p^{k/2}) for chi_5, chi_{-4}, chi_{-8}, chi_8
 *   K3_2re, K3_2im                 : the same for the complex chi_2 mod 5 (chi_2(2)=i)
 *   S_m8, S_8, S_m4                : log prod (1-chi(p)p^{-1/2})^{-1} for chi_{-8}, chi_8, chi_{-4}
 * Characters mod 8 (a = 1,3,5,7):  chi_{-4} = +,-,+,-   chi_{-8} = +,+,-,-   chi_8 = +,-,-,+
 *
 * build: gcc -O3 -march=native -fopenmp -o scan_mod58_mp scan_mod58_mp.c -lm
 * run  : ./scan_mod58_mp N per node num_nodes threads [chunk_log2]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

typedef struct { double s, c; } KS;
static inline void kadd(KS *k, double x){ double y = x - k->c, t = k->s + y; k->c = (t - k->s) - y; k->s = t; }

#define NV 21
enum { P5_1, P5_2, P5_3, P5_4, P8_1, P8_3, P8_5, P8_7, ME, M2, H4,
       K3_5, K3_2re, K3_2im, K3_m4, K3_m8, K3_8, S_m8, S_8, S_m4, PI5 };
/* PI5 is unused padding kept for alignment of the column count; see header string */
static const char *HEADER = "x,pi_x,P5_1,P5_2,P5_3,P5_4,P8_1,P8_3,P8_5,P8_7,ME,M2,H4,K3_5,K3_2re,K3_2im,K3_m4,K3_m8,K3_8,S_m8,S_8,S_m4,unused";

typedef struct { KS v[NV]; uint64_t cnt; } ACC;
typedef struct { uint64_t cnt; double v[NV]; } SNAP;
static inline void acc_snap(const ACC *a, SNAP *o){ o->cnt = a->cnt; for (int i = 0; i < NV; i++) o->v[i] = a->v[i].s; }

static inline void acc_prime(ACC *A, uint64_t p){
    A->cnt++;
    double w = 1.0 / (double)p, u = 1.0 / sqrt((double)p);
    kadd(&A->v[ME], w);
    if (p <= 1000000ULL) kadd(&A->v[M2], w + log1p(-w));
    else { double w2 = w*w, w3 = w2*w, w4 = w2*w2; kadd(&A->v[M2], -w2/2 - w3/3 - w4/4); }

    /* closed forms shared by all real characters */
    double odd3, evn4, oddA, evnA;                 /* real chi: atanh(u)-u, -log(1-w)/2 - w/2 ; alternating: atan(u)-u, -log(1+w)/2 + w/2 */
    if (p <= 1000000ULL){
        odd3 = atanh(u) - u;  evn4 = -0.5*log1p(-w) - 0.5*w;
        oddA = atan(u)  - u;  evnA = -0.5*log1p( w) + 0.5*w;
    } else {
        double u3 = u*w, u5 = u3*w, u7 = u5*w, w2 = w*w, w3 = w2*w, w4 = w2*w2;
        odd3 = u3/3 + u5/5 + u7/7;  evn4 = w2/4 + w3/6 + w4/8;
        oddA = -u3/3 + u5/5 - u7/7; evnA = w2/4 - w3/6 + w4/8;
    }
    double halfw = 0.5*w;

    /* ---- mod 5 ---- */
    if (p != 5){
        uint64_t r = p % 5;
        int idx = (r == 1) ? P5_1 : (r == 2) ? P5_2 : (r == 3) ? P5_3 : P5_4;
        kadd(&A->v[idx], u);
        double c5 = (r == 1 || r == 4) ? 1.0 : -1.0;          /* chi_5 = Legendre */
        kadd(&A->v[H4], c5 * w);
        kadd(&A->v[K3_5], c5*odd3 + evn4);
        if (r == 1)      kadd(&A->v[K3_2re],  odd3 + evn4);    /* chi_2 = +1 */
        else if (r == 4) kadd(&A->v[K3_2re], -odd3 + evn4);    /* chi_2 = -1 */
        else {                                                  /* chi_2 = +i (r=2) or -i (r=3): chi_2^2 = -1 */
            kadd(&A->v[K3_2re], evnA);
            kadd(&A->v[K3_2im], (r == 2) ? oddA : -oddA);
        }
    }
    /* ---- mod 8 ---- */
    if (p != 2){
        uint64_t r = p & 7;
        int idx = (r == 1) ? P8_1 : (r == 3) ? P8_3 : (r == 5) ? P8_5 : P8_7;
        kadd(&A->v[idx], u);
        double cm4 = (r == 1 || r == 5) ? 1.0 : -1.0;
        double cm8 = (r == 1 || r == 3) ? 1.0 : -1.0;
        double c8  = (r == 1 || r == 7) ? 1.0 : -1.0;
        kadd(&A->v[K3_m4], cm4*odd3 + evn4);
        kadd(&A->v[K3_m8], cm8*odd3 + evn4);
        kadd(&A->v[K3_8],  c8 *odd3 + evn4);
        if (p <= 1000000ULL){
            kadd(&A->v[S_m4], -log1p(-cm4*u));
            kadd(&A->v[S_m8], -log1p(-cm8*u));
            kadd(&A->v[S_8],  -log1p(-c8 *u));
        } else {
            kadd(&A->v[S_m4], cm4*(u + odd3) + (halfw + evn4));
            kadd(&A->v[S_m8], cm8*(u + odd3) + (halfw + evn4));
            kadd(&A->v[S_8],  c8 *(u + odd3) + (halfw + evn4));
        }
    }
}

static uint64_t cps[8300]; static int ncp = 0;
static void build_cps(uint64_t N, int per){
    uint64_t prev = 0;
    for (int k = 0; k < 8192; k++){ uint64_t v = (uint64_t)(pow(10.0, 2.0 + (double)k/per) + 0.5);
        if (v >= N) break; if (v > prev){ cps[ncp++] = v; prev = v; } if (ncp >= 8200) break; }
    cps[ncp++] = N;
}
static uint64_t *bq = NULL; static uint32_t nb = 0;
static void build_base(uint64_t N){
    uint64_t lim = (uint64_t)sqrt((double)N); while ((lim+1)*(lim+1) <= N) lim++; while (lim > 1 && lim*lim > N) lim--;
    uint8_t *small = calloc(lim + 2, 1);
    for (uint64_t i = 2; i*i <= lim; i++) if (!small[i]) for (uint64_t j = i*i; j <= lim; j += i) small[j] = 1;
    bq = malloc(sizeof(uint64_t) * (lim/2 + 2));
    for (uint64_t i = 3; i <= lim; i += 2) if (!small[i]) bq[nb++] = i;
    free(small);
}
static void print_row(uint64_t x, uint64_t cnt, const double *v){
    printf("%llu,%llu", (unsigned long long)x, (unsigned long long)cnt);
    for (int i = 0; i < NV; i++) printf(",%.15g", v[i]);
    printf("\n");
}

int main(int argc, char **argv){
    if (argc < 5){ fprintf(stderr, "usage: %s N per node of [threads] [chunk_log2]\n", argv[0]); return 1; }
    uint64_t N = strtoull(argv[1], NULL, 10); int per = atoi(argv[2]); int node = atoi(argv[3]); int of = atoi(argv[4]);
    int nthreads = (argc > 5) ? atoi(argv[5]) : 0; int clog = (argc > 6) ? atoi(argv[6]) : 31;
    if (nthreads > 0) omp_set_num_threads(nthreads);
    build_cps(N, per); build_base(N);
    uint64_t span = N - 1, lo = 2 + (span*(uint64_t)node)/of, hi = 2 + (span*(uint64_t)(node+1))/of;
    const uint64_t W = 1ULL << clog; uint64_t C = (hi - lo + W - 1)/W;
    int cpA = 0; while (cpA < ncp && cps[cpA] < lo) cpA++;
    int cpB = cpA; while (cpB < ncp && cps[cpB] < hi) cpB++;
    SNAP *TOT = calloc(C, sizeof(SNAP)); SNAP *SNAPS = calloc(ncp, sizeof(SNAP));
    const uint64_t SEGODD = 1ULL << 21; double t0 = omp_get_wtime(); uint64_t chunks_done = 0;
#pragma omp parallel
    {
        uint8_t *seg = malloc(SEGODD); uint64_t *nxt = malloc(sizeof(uint64_t)*(nb+1));
#pragma omp for schedule(dynamic, 1)
        for (uint64_t c = 0; c < C; c++){
            uint64_t clo = lo + c*W, chi_ = clo + W; if (chi_ > hi) chi_ = hi;
            ACC A; memset(&A, 0, sizeof(A));
            if (clo <= 2 && 2 < chi_) acc_prime(&A, 2);
            int ci = cpA; while (ci < cpB && cps[ci] < clo) ci++;
            int ce = ci;  while (ce < cpB && cps[ce] < chi_) ce++;
            uint32_t nbq = 0; while (nbq < nb && bq[nbq]*bq[nbq] < chi_) nbq++;
            for (uint32_t i = 0; i < nbq; i++){ uint64_t q = bq[i], q2 = q*q, s = ((clo + q - 1)/q)*q;
                if (s < q2) s = q2; if (!(s & 1)) s += q; nxt[i] = s; }
            uint64_t lo_odd = (clo < 3) ? 3 : clo; if (!(lo_odd & 1)) lo_odd++;
            for (uint64_t low = lo_odd; low < chi_; low += 2*SEGODD){
                uint64_t high = low + 2*SEGODD; if (high > chi_) high = chi_;
                uint64_t nod = (high - low + 1)/2; memset(seg, 0, nod);
                for (uint32_t i = 0; i < nbq; i++){ uint64_t q = bq[i], m = nxt[i];
                    for (; m < high; m += 2*q) seg[(m - low) >> 1] = 1; nxt[i] = m; }
                for (uint64_t i = 0; i < nod; i++) if (!seg[i]){ uint64_t p = low + 2*i;
                    while (ci < ce && p > cps[ci]){ acc_snap(&A, &SNAPS[ci]); ci++; } acc_prime(&A, p); }
            }
            while (ci < ce){ acc_snap(&A, &SNAPS[ci]); ci++; }
            acc_snap(&A, &TOT[c]);
            uint64_t done_now;
#pragma omp atomic capture
            done_now = ++chunks_done;
            { uint64_t step = (C >= 100) ? C/100 : 1;
              if (done_now % step == 0 || done_now == C){ double el = omp_get_wtime() - t0, fr = (double)done_now/(double)C;
                fprintf(stderr, "node %d progress %5.1f%%  (%llu/%llu chunks)  %.0fs elapsed, ~%.0fs left\n", node, 100.0*fr,
                    (unsigned long long)done_now, (unsigned long long)C, el, el*(1.0-fr)/(fr > 0 ? fr : 1)); } }
        }
        free(seg); free(nxt);
    }
    KS R[NV]; memset(R, 0, sizeof(R)); uint64_t rcnt = 0, cidx = 0; double row[NV];
    printf("#NODE,%d,%d,%llu,%d\n#RANGE,%llu,%llu\n%s\n", node, of, (unsigned long long)N, per,
        (unsigned long long)lo, (unsigned long long)hi, HEADER);
    for (int j = cpA; j < cpB; j++){
        uint64_t owner = (cps[j] - lo)/W;
        while (cidx < owner){ for (int t = 0; t < NV; t++) kadd(&R[t], TOT[cidx].v[t]); rcnt += TOT[cidx].cnt; cidx++; }
        for (int t = 0; t < NV; t++) row[t] = R[t].s + SNAPS[j].v[t];
        print_row(cps[j], rcnt + SNAPS[j].cnt, row);
    }
    while (cidx < C){ for (int t = 0; t < NV; t++) kadd(&R[t], TOT[cidx].v[t]); rcnt += TOT[cidx].cnt; cidx++; }
    printf("#TOTAL,%llu", (unsigned long long)rcnt); for (int t = 0; t < NV; t++) printf(",%.17g", R[t].s); printf("\n");
    fprintf(stderr, "node %d/%d done: range [%llu,%llu)  count=%llu  %.1fs (%d threads)\n", node, of,
        (unsigned long long)lo, (unsigned long long)hi, (unsigned long long)rcnt, omp_get_wtime() - t0, omp_get_max_threads());
    free(TOT); free(SNAPS); free(bq); return 0;
}
