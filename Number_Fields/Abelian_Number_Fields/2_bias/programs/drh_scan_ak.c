/* drh_scan_ak.c — both sides of Aoki's theorem (Note 2026-09-04, p.19) for chi mod 4.
 *
 * Accumulates, over primes p <= x (checkpoints 128/decade):
 *   SB4 = sum_{p odd} chi(p)/sqrt(p)                      (k=1 bias sum; data side)
 *   ME  = sum_{all p} 1/p                                  (Mertens)
 *   M2  = sum_{all p} [1/p + log(1-1/p)]                   (theory side, converges)
 *   K3  = sum_{p odd} sum_{k>=3} chi(p)^k/(k p^{k/2})      (theory side, converges)
 *   S4  = sum_{p odd} -log(1-chi(p)p^{-1/2})               (independent recompute, cross-check)
 * Data side  A(x) = 1/sqrt2 - SB4 - (1/2)loglog x
 * Theory side B(x) = g/2 + 1/sqrt2 - log(2)/2 + M2/2 - log L(1/2,chi) - 1/4 + K3
 * (computed in the analysis script; this program only writes the raw sums)
 *
 * build: gcc -O2 -march=native -o drh_scan_ak drh_scan_ak.c -lm
 * run  : ./drh_scan_ak 100000000000 128 > ak_traj_1e11.csv 2> ak.log
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

typedef struct { double s, c; } KS;
static inline void kadd(KS *k, double x){
    double y = x - k->c, t = k->s + y;
    k->c = (t - k->s) - y; k->s = t;
}
static KS SB4, ME, M2, K3, S4;
static uint64_t cnt = 0;

static inline void process_prime(uint64_t p){
    cnt++;
    if (p == 2){
        kadd(&ME, 0.5);
        kadd(&M2, 0.5 + log(0.5));
        return;                                   /* chi(2)=0: not in SB4/K3/S4 */
    }
    double w = 1.0 / (double)p;
    double u = 1.0 / sqrt((double)p);
    double c = ((p & 3ULL) == 1) ? 1.0 : -1.0;    /* chi_{-4}(p) */
    kadd(&ME, w);
    if (p <= 1000000ULL){                          /* exact path */
        kadd(&M2, w + log1p(-w));
        kadd(&SB4, c * u);
        kadd(&K3, c * (atanh(u) - u) + (-0.5 * log1p(-w) - 0.5 * w));
        kadd(&S4, -log1p(-c * u));
    } else {                                       /* series path: u<=1e-3, w<=1e-6 */
        double w2 = w * w, w3 = w2 * w, w4 = w2 * w2;
        double u3 = u * w, u5 = u3 * w, u7 = u5 * w;          /* u^3=u*w etc. */
        kadd(&M2, -w2 / 2 - w3 / 3 - w4 / 4);
        kadd(&SB4, c * u);
        double odd3 = u3 / 3 + u5 / 5 + u7 / 7;               /* atanh(u)-u   */
        double evn4 = w2 / 4 + w3 / 6 + w4 / 8;               /* -log(1-w)/2 - w/2 */
        kadd(&K3, c * odd3 + evn4);
        kadd(&S4, c * (u + odd3) + (0.5 * w + evn4));
    }
}

static uint64_t cps[4200]; static int ncp = 0;
static void build_cps(uint64_t N, int per){
    uint64_t prev = 0;
    for (int k = 0; k < 4096; k++){
        uint64_t v = (uint64_t)(pow(10.0, 2.0 + (double)k / per) + 0.5);
        if (v >= N) break;
        if (v > prev){ cps[ncp++] = v; prev = v; }
        if (ncp >= 4100) break;
    }
    cps[ncp++] = N;
}
static int cpi = 0;
static inline void maybe_emit(uint64_t p){
    while (cpi < ncp && p > cps[cpi]){
        printf("%llu,%llu,%.15g,%.15g,%.15g,%.15g,%.15g\n",
            (unsigned long long)cps[cpi], (unsigned long long)cnt,
            SB4.s, ME.s, M2.s, K3.s, S4.s);
        cpi++;
    }
}

int main(int argc, char **argv){
    uint64_t N = (argc > 1) ? strtoull(argv[1], NULL, 10) : 100000000000ULL;
    int per    = (argc > 2) ? atoi(argv[2]) : 128;
    build_cps(N, per);
    printf("x,pi_x,SB4,ME,M2,K3,S4chk\n");
    clock_t t0 = clock();

    uint64_t lim = (uint64_t)sqrt((double)N);
    while ((lim + 1) * (lim + 1) <= N) lim++;
    uint8_t *small = calloc(lim + 2, 1);
    for (uint64_t i = 2; i * i <= lim; i++)
        if (!small[i]) for (uint64_t j = i * i; j <= lim; j += i) small[j] = 1;
    uint64_t nb = 0, *bq = malloc(sizeof(uint64_t) * (lim / 2 + 2));
    for (uint64_t i = 3; i <= lim; i += 2) if (!small[i]) bq[nb++] = i;

    maybe_emit(2); process_prime(2);
    const uint64_t SEGODD = 1ULL << 21;
    uint8_t *seg = malloc(SEGODD);
    uint64_t *nxt = malloc(sizeof(uint64_t) * (nb + 1));
    for (uint64_t i = 0; i < nb; i++){
        uint64_t q = bq[i], s = q * q;
        nxt[i] = s;
    }
    uint64_t done_mark = 0;
    for (uint64_t low = 3; low <= N; low += 2 * SEGODD){
        uint64_t high = low + 2 * SEGODD; if (high > N + 1) high = N + 1;
        uint64_t nod = (high - low + 1) / 2;
        memset(seg, 0, nod);
        for (uint64_t i = 0; i < nb; i++){
            uint64_t q = bq[i];
            if (q * q >= high) break;
            uint64_t m = nxt[i];
            for (; m < high; m += 2 * q) seg[(m - low) >> 1] = 1;
            nxt[i] = m;
        }
        for (uint64_t i = 0; i < nod; i++){
            if (!seg[i]){
                uint64_t p = low + 2 * i;
                maybe_emit(p);
                process_prime(p);
            }
        }
        if (++done_mark % 512 == 0)
            fprintf(stderr, "progress %5.1f%%  x~%.3g  %.0fs\n",
                100.0 * (double)(high) / (double)N, (double)high,
                (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    while (cpi < ncp){
        printf("%llu,%llu,%.15g,%.15g,%.15g,%.15g,%.15g\n",
            (unsigned long long)cps[cpi], (unsigned long long)cnt,
            SB4.s, ME.s, M2.s, K3.s, S4.s);
        cpi++;
    }
    fprintf(stderr, "done: pi(%llu) = %llu   %.1fs\n",
        (unsigned long long)N, (unsigned long long)cnt,
        (double)(clock() - t0) / CLOCKS_PER_SEC);
    free(seg); free(nxt); free(bq); free(small);
    return 0;
}
