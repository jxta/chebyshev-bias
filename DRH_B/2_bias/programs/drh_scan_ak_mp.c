/* drh_scan_mp.c
 * Parallel version of drh_scan_ak.c (Aoki theorem sums) for multi-node runs (e.g. mdx: 5 nodes x 152 vCPU).
 *
 * Design: no MPI. The range [2, N] is split into `of` contiguous node ranges;
 * each node runs this binary independently (OpenMP inside the node) and writes
 * a partial CSV. merge_traj.py turns the partial files into one global
 * trajectory identical in format to the serial drh_scan.c output.
 *
 * Within a node, the range is cut into fixed-width chunks processed by OpenMP
 * threads (dynamic schedule). Each chunk accumulates its own Kahan sums and
 * records snapshots at the log-spaced checkpoints it owns; a sequential prefix
 * pass then produces node-relative cumulative values. Results are therefore
 * bitwise independent of the thread scheduling (deterministic).
 *
 * build: gcc -O3 -march=native -fopenmp -o drh_scan_mp drh_scan_mp.c -lm
 * run  : ./drh_scan_mp N per node of [threads] [chunk_log2]
 *        ./drh_scan_mp 1000000000000000 128 0 5 152 31 > part0.csv 2> log0.txt
 * merge: python3 merge_traj.py part*.csv -o drh_traj.csv
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

typedef struct { double s, c; } KS;            /* Kahan sum */
static inline void kadd(KS *k, double x){
    double y = x - k->c;
    double t = k->s + y;
    k->c = (t - k->s) - y;
    k->s = t;
}

typedef struct {
    KS SB4, ME, M2, K3, S4;
    uint64_t cnt;
} ACC;
#define NV 5

typedef struct {
    uint64_t cnt;
    double v[NV];                              /* SB4,ME,M2,K3,S4chk */
} SNAP;

static inline void acc_snap(const ACC *a, SNAP *o){
    o->cnt = a->cnt;
    o->v[0]=a->SB4.s; o->v[1]=a->ME.s; o->v[2]=a->M2.s;
    o->v[3]=a->K3.s;  o->v[4]=a->S4.s;
}

/* identical arithmetic to serial drh_scan_ak.c (do not change thresholds) */
static inline void acc_prime(ACC *A, uint64_t p){
    A->cnt++;
    if (p == 2){
        kadd(&A->ME, 0.5);
        kadd(&A->M2, 0.5 + log(0.5));
        return;                                /* chi(2)=0: not in SB4/K3/S4 */
    }
    double w = 1.0 / (double)p;
    double u = 1.0 / sqrt((double)p);
    double c = ((p & 3ULL) == 1) ? 1.0 : -1.0;
    kadd(&A->ME, w);
    if (p <= 1000000ULL){
        kadd(&A->M2, w + log1p(-w));
        kadd(&A->SB4, c * u);
        kadd(&A->K3, c * (atanh(u) - u) + (-0.5 * log1p(-w) - 0.5 * w));
        kadd(&A->S4, -log1p(-c * u));
    } else {
        double w2 = w * w, w3 = w2 * w, w4 = w2 * w2;
        double u3 = u * w, u5 = u3 * w, u7 = u5 * w;
        kadd(&A->M2, -w2 / 2 - w3 / 3 - w4 / 4);
        kadd(&A->SB4, c * u);
        double odd3 = u3 / 3 + u5 / 5 + u7 / 7;
        double evn4 = w2 / 4 + w3 / 6 + w4 / 8;
        kadd(&A->K3, c * odd3 + evn4);
        kadd(&A->S4, c * (u + odd3) + (0.5 * w + evn4));
    }
}

/* ---- checkpoints (global, log-spaced) ---- */
static uint64_t cps[8300];
static int ncp = 0;
static void build_cps(uint64_t N, int per){
    uint64_t prev = 0;
    for (int k = 0; k < 8192; k++){
        double e = 2.0 + (double)k / per;
        uint64_t v = (uint64_t)(pow(10.0, e) + 0.5);
        if (v >= N) break;
        if (v > prev){ cps[ncp++] = v; prev = v; }
        if (ncp >= 8200) break;
    }
    cps[ncp++] = N;
}

/* base primes up to sqrt(N) */
static uint64_t *bq = NULL;
static uint32_t nb = 0;
static void build_base(uint64_t N){
    uint64_t lim = (uint64_t)sqrt((double)N);
    while ((lim+1)*(lim+1) <= N) lim++;
    while (lim > 1 && lim*lim > N) lim--;
    uint8_t *small = calloc(lim + 2, 1);
    for (uint64_t i = 2; i*i <= lim; i++)
        if (!small[i]) for (uint64_t j = i*i; j <= lim; j += i) small[j] = 1;
    bq = malloc(sizeof(uint64_t) * (lim/2 + 2));
    for (uint64_t i = 3; i <= lim; i += 2) if (!small[i]) bq[nb++] = i;
    free(small);
}

int main(int argc, char **argv){
    if (argc < 5){
        fprintf(stderr, "usage: %s N per node of [threads] [chunk_log2]\n", argv[0]);
        return 1;
    }
    uint64_t N   = strtoull(argv[1], NULL, 10);
    int per      = atoi(argv[2]);
    int node     = atoi(argv[3]);
    int of       = atoi(argv[4]);
    int nthreads = (argc > 5) ? atoi(argv[5]) : 0;
    int clog     = (argc > 6) ? atoi(argv[6]) : 31;
    if (nthreads > 0) omp_set_num_threads(nthreads);

    build_cps(N, per);
    build_base(N);

    /* node range [lo, hi) tiling [2, N+1) */
    uint64_t span = N - 1;
    uint64_t lo = 2 + (span * (uint64_t)node) / of;
    uint64_t hi = 2 + (span * (uint64_t)(node + 1)) / of;

    const uint64_t W = 1ULL << clog;                    /* chunk width (numbers) */
    uint64_t C = (hi - lo + W - 1) / W;                 /* chunks in this node   */

    /* checkpoint ownership: cps[j] in [lo, hi) */
    int cpA = 0; while (cpA < ncp && cps[cpA] <  lo) cpA++;
    int cpB = cpA; while (cpB < ncp && cps[cpB] < hi) cpB++;

    SNAP *TOT   = calloc(C, sizeof(SNAP));
    SNAP *SNAPS = calloc(ncp, sizeof(SNAP));

    const uint64_t SEGODD = 1ULL << 21;                 /* odd slots per sub-segment (2 MB) */
    double t0 = omp_get_wtime();
    uint64_t chunks_done = 0;

#pragma omp parallel
    {
        uint8_t  *seg = malloc(SEGODD);
        uint64_t *nxt = malloc(sizeof(uint64_t) * (nb + 1));

#pragma omp for schedule(dynamic, 1)
        for (uint64_t c = 0; c < C; c++){
            uint64_t clo = lo + c * W;
            uint64_t chi_ = clo + W; if (chi_ > hi) chi_ = hi;

            ACC A; memset(&A, 0, sizeof(A));
            if (clo <= 2 && 2 < chi_ && N >= 2) acc_prime(&A, 2);

            /* checkpoints owned by this chunk */
            int ci = cpA; while (ci < cpB && cps[ci] < clo) ci++;
            int ce = ci;  while (ce < cpB && cps[ce] < chi_) ce++;

            /* base primes relevant to this chunk, and their first odd multiple */
            uint32_t nbq = 0;
            while (nbq < nb && bq[nbq] * bq[nbq] < chi_) nbq++;
            for (uint32_t i = 0; i < nbq; i++){
                uint64_t q = bq[i], q2 = q * q;
                uint64_t s = ((clo + q - 1) / q) * q;
                if (s < q2) s = q2;
                if (!(s & 1)) s += q;
                nxt[i] = s;
            }

            uint64_t lo_odd = (clo < 3) ? 3 : clo;
            if (!(lo_odd & 1)) lo_odd++;

            for (uint64_t low = lo_odd; low < chi_; low += 2 * SEGODD){
                uint64_t high = low + 2 * SEGODD;       /* exclusive */
                if (high > chi_) high = chi_;
                uint64_t nod = (high - low + 1) / 2;    /* odd numbers in [low, high) */
                memset(seg, 0, nod);
                for (uint32_t i = 0; i < nbq; i++){
                    uint64_t q = bq[i], m = nxt[i];
                    for (; m < high; m += 2 * q) seg[(m - low) >> 1] = 1;
                    nxt[i] = m;
                }
                for (uint64_t i = 0; i < nod; i++){
                    if (!seg[i]){
                        uint64_t p = low + 2 * i;
                        while (ci < ce && p > cps[ci]){ acc_snap(&A, &SNAPS[ci]); ci++; }
                        acc_prime(&A, p);
                    }
                }
            }
            while (ci < ce){ acc_snap(&A, &SNAPS[ci]); ci++; }
            acc_snap(&A, &TOT[c]);

uint64_t done_now;
#pragma omp atomic capture
            done_now = ++chunks_done;
            {   /* print on every 1% boundary (any thread; exact because 'capture' gives unique values) */
                uint64_t step = (C >= 100) ? C / 100 : 1;
                if (done_now % step == 0 || done_now == C) {
                    double el = omp_get_wtime() - t0;
                    double frac = (double)done_now / (double)C;
                    fprintf(stderr, "node %d progress %5.1f%%  (%llu/%llu chunks)  %.0fs elapsed, ~%.0fs left\n",
                            node, 100.0 * frac,
                            (unsigned long long)done_now, (unsigned long long)C,
                            el, el * (1.0 - frac) / (frac > 0 ? frac : 1));
                }
            }
        }
        free(seg); free(nxt);
    }

    /* sequential prefix over chunks -> node-relative cumulative values */
    KS R[NV]; memset(R, 0, sizeof(R));
    uint64_t rcnt = 0;
    uint64_t cidx = 0;
    printf("#NODE,%d,%d,%llu,%d\n", node, of, (unsigned long long)N, per);
    printf("#RANGE,%llu,%llu\n", (unsigned long long)lo, (unsigned long long)hi);
    printf("x,pi_x,SB4,ME,M2,K3,S4chk\n");
    for (int j = cpA; j < cpB; j++){
        uint64_t owner = (cps[j] - lo) / W;             /* chunk owning this cp */
        while (cidx < owner){                           /* fold finished chunks */
            for (int t = 0; t < NV; t++) kadd(&R[t], TOT[cidx].v[t]);
            rcnt += TOT[cidx].cnt;
            cidx++;
        }
        printf("%llu,%llu,%.15g,%.15g,%.15g,%.15g,%.15g\n",
            (unsigned long long)cps[j],
            (unsigned long long)(rcnt + SNAPS[j].cnt),
            R[0].s + SNAPS[j].v[0], R[1].s + SNAPS[j].v[1],
            R[2].s + SNAPS[j].v[2], R[3].s + SNAPS[j].v[3],
            R[4].s + SNAPS[j].v[4]);
    }
    while (cidx < C){
        for (int t = 0; t < NV; t++) kadd(&R[t], TOT[cidx].v[t]);
        rcnt += TOT[cidx].cnt;
        cidx++;
    }
    printf("#TOTAL,%llu,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        (unsigned long long)rcnt, R[0].s, R[1].s, R[2].s, R[3].s, R[4].s);
    fprintf(stderr, "node %d/%d done: range [%llu,%llu)  count=%llu  %.1fs (%d threads)\n",
        node, of, (unsigned long long)lo, (unsigned long long)hi,
        (unsigned long long)rcnt, omp_get_wtime() - t0, omp_get_max_threads());
    free(TOT); free(SNAPS); free(bq);
    return 0;
}
