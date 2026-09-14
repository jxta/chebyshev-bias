/* drh_scan.c
 * Partial Euler products of Dirichlet L-functions at s = 1/2:
 *     log prod_{p<=x} (1 - chi(p) p^{-1/2})^{-1}
 * accumulated simultaneously for the primitive characters
 *   chi4  : mod 4  (chi(3)=-1)                  [real, nu=1]
 *   chi3  : mod 3  (chi(2)=-1)                  [real, nu=1]
 *   chi5  : mod 5  Legendre symbol              [real, nu=1]
 *   chi5c : mod 5  order 4, chi(2)=+i           [complex, nu=0]
 *   chi7  : mod 7  Legendre symbol              [real, nu=1]
 * over primes p <= N via an odd-only segmented sieve.
 * CSV trajectory at ~32 log-spaced checkpoints per decade goes to stdout.
 *
 * build: gcc -O2 -o drh_scan drh_scan.c -lm
 * run:   ./drh_scan 10000000000 > drh_traj.csv
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

typedef struct { double s, c; } KS;            /* Kahan sum */
static inline void kadd(KS *k, double x){
    double y = x - k->c;
    double t = k->s + y;
    k->c = (t - k->s) - y;
    k->s = t;
}

static KS S4, S3, S5, S5r, S5i, S7;
static uint64_t primecount = 0;

static inline void process_prime(uint64_t p){
    primecount++;
    double u = 1.0 / sqrt((double)p);
    if (p > 100000000ULL){
        /* series fast path: u < 1e-4, truncation error < u^6/6 ~ 2e-25 per term */
        double u2 = u*u, u3 = u2*u, u4 = u2*u2, u5 = u2*u3;
        double tp =  u + 0.5*u2 + u3*(1.0/3.0) + 0.25*u4 + 0.2*u5;  /* -log(1-u) */
        double tm = -u + 0.5*u2 - u3*(1.0/3.0) + 0.25*u4 - 0.2*u5;  /* -log(1+u) */
        uint64_t r;
        r = p & 3ULL;  kadd(&S4, (r==1) ? tp : tm);
        r = p % 3ULL;  kadd(&S3, (r==1) ? tp : tm);
        r = p % 5ULL;
        kadd(&S5, (r==1 || r==4) ? tp : tm);
        if (r==1)      kadd(&S5r, tp);
        else if (r==4) kadd(&S5r, tm);
        else {         /* chi = +i (r=2) or -i (r=3): -log(1 -/+ iu) */
            kadd(&S5r, -0.5*u2 + 0.25*u4);
            double at = u - u3*(1.0/3.0) + 0.2*u5;                   /* atan(u) */
            kadd(&S5i, (r==2) ? at : -at);
        }
        r = p % 7ULL;  kadd(&S7, (r==1 || r==2 || r==4) ? tp : tm);
    } else {
        double lp = -log1p(-u);   /* chi(p) = +1 */
        double lm = -log1p( u);   /* chi(p) = -1 */
        uint64_t r;
        if (p != 2){ r = p & 3ULL; kadd(&S4, (r==1) ? lp : lm); }
        if (p != 3){ r = p % 3ULL; kadd(&S3, (r==1) ? lp : lm); }
        if (p != 5){
            r = p % 5ULL;
            kadd(&S5, (r==1 || r==4) ? lp : lm);
            if (r==1)      kadd(&S5r, lp);
            else if (r==4) kadd(&S5r, lm);
            else {
                kadd(&S5r, -0.5*log1p(u*u));
                double at = atan(u);
                kadd(&S5i, (r==2) ? at : -at);
            }
        }
        if (p != 7){ r = p % 7ULL; kadd(&S7, (r==1 || r==2 || r==4) ? lp : lm); }
    }
}

static uint64_t cps[4200];
static int ncp = 0, cpi = 0;

static int g_per = 32;
static void build_cps(uint64_t N){
    const int per = g_per;
    /* checkpoints per decade set in main via g_per */
    uint64_t prev = 0;
    for (int k = 0; k < 4096; k++){
        double e = 2.0 + (double)k / per;
        uint64_t v = (uint64_t)(pow(10.0, e) + 0.5);
        if (v >= N) break;
        if (v > prev){ cps[ncp++] = v; prev = v; }
        if (ncp >= 4100) break;
    }
    cps[ncp++] = N;
}

static void emit(uint64_t x){
    printf("%llu,%llu,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n",
        (unsigned long long)x, (unsigned long long)primecount,
        S4.s, S3.s, S5.s, S5r.s, S5i.s, S7.s);
}

static inline void maybe_emit(uint64_t p){
    while (cpi < ncp && p > cps[cpi]){ emit(cps[cpi]); cpi++; }
}

int main(int argc, char **argv){
    uint64_t N = (argc > 1) ? strtoull(argv[1], NULL, 10) : 10000000000ULL;
    if (argc > 2) g_per = atoi(argv[2]);
    build_cps(N);

    uint64_t lim = (uint64_t)sqrt((double)N);
    while ((lim+1)*(lim+1) <= N) lim++;
    while (lim > 1 && lim*lim > N) lim--;

    /* base sieve up to lim */
    uint8_t *small = calloc(lim + 2, 1);
    for (uint64_t i = 2; i*i <= lim; i++)
        if (!small[i]) for (uint64_t j = i*i; j <= lim; j += i) small[j] = 1;
    uint32_t nb = 0;
    uint64_t *bq = malloc(sizeof(uint64_t) * (lim/2 + 2));
    for (uint64_t i = 3; i <= lim; i += 2) if (!small[i]) bq[nb++] = i;
    free(small);

    printf("x,pi_x,S4,S3,S5,S5c_re,S5c_im,S7\n");
    clock_t t0 = clock();

    if (N >= 2){ maybe_emit(2); process_prime(2); }

    const uint64_t SEG = 1u << 20;            /* odd slots per segment (1 MB) */
    uint8_t *seg = malloc(SEG);
    uint64_t nseg = 0;
    for (uint64_t low = 3; low <= N; low += 2*SEG){
        uint64_t high = low + 2*SEG;          /* exclusive bound */
        if (high > N + 1) high = N + 1;
        uint64_t nod = (high - low + 1) / 2;  /* odd numbers in [low, high) */
        memset(seg, 0, nod);
        for (uint32_t i = 0; i < nb; i++){
            uint64_t q = bq[i];
            uint64_t q2 = q * q;
            if (q2 >= high) break;
            uint64_t start = ((low + q - 1) / q) * q;
            if (start < q2) start = q2;
            if (!(start & 1)) start += q;
            for (uint64_t m = start; m < high; m += 2*q)
                seg[(m - low) >> 1] = 1;
        }
        for (uint64_t i = 0; i < nod; i++){
            if (!seg[i]){
                uint64_t p = low + 2*i;
                maybe_emit(p);
                process_prime(p);
            }
        }
        if ((++nseg & 255) == 0)
            fprintf(stderr, "progress %5.1f%%  x~%.3g  %.0fs\n",
                100.0 * (double)high / (double)N, (double)high,
                (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    while (cpi < ncp) emit(cps[cpi++]);
    fprintf(stderr, "done: pi(%llu) = %llu   %.1fs\n",
        (unsigned long long)N, (unsigned long long)primecount,
        (double)(clock() - t0) / CLOCKS_PER_SEC);
    free(seg); free(bq);
    return 0;
}
