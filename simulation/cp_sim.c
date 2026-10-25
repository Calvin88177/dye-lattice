/*
 * cp_sim.c  --  continuous-time contact process on a dye lattice, with the three
 * imperfections of the real system:
 *     field      h   : spontaneous wake-up of a dark dye            (rate h)
 *     bleaching  pb  : each ON period ends in permanent death with probability pb
 *     disorder       : missing sites (occupancy p) and quenched position jitter
 *
 * Units: time in 1/k_off (one mean ON period).  An ON dye goes dark at rate 1.
 * A dark dye i wakes at rate  h + sum_{j ON} r_ij,   r_ij = r1 * (a/d_ij)^6,
 * where a is the nearest-neighbour distance.  The pair "hand-over chance" is
 * B = r/(1+r).
 *
 * Exact event-driven (Gillespie) dynamics; a Fenwick tree holds the per-site rates.
 *
 * Usage: cp_sim key=value ...   (see parse section for keys)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ---------- RNG: xoshiro256** ---------- */
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t nextu(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}
static inline double urand(void) { return (nextu() >> 11) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t x) {
    for (int i = 0; i < 4; i++) {
        x += 0x9E3779B97F4A7C15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        s[i] = z ^ (z >> 31);
    }
}
static double gauss(void) {
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ---------- parameters ---------- */
static int    lat = 0;        /* 0 chain, 1 square, 2 triangular, 3 centred-rectangular "pixel" */
static int    Lx = 100, Ly = 1;
static int    periodic = 1;
static double cut = 1.05;     /* neighbour cutoff in units of a (1.05 = nearest neighbours only) */
static double r1 = 1.0;       /* wake rate per ON nearest neighbour, units of k_off */
static double hfield = 0.0;
static double pb = 0.0;
static double occ = 1.0;      /* site occupancy */
static double sigma = 0.0;    /* quenched position jitter per coordinate, units of a */
static double tmax = 1000.0;
static double tmin = 0.01;
static int    K = 200;        /* log-spaced time grid points */
static int    runs = 1;
static int    mode = 0;       /* 0 quench (all ON), 1 single seed, 2 avalanches (all OFF + field) */
static uint64_t seed0 = 12345;
static char   outp[512] = "out";
static double fdead_max = 0.05;   /* avalanche mode: stop when this fraction has bleached */
static long   maxbursts = 2000000;
static int    newlat = 1;     /* draw new disorder for every run */
static double px_P = 10.88, px_b = 5.30;
static int    seedrand = 0;   /* seed mode: 1 = seed at a random present site instead of the centre */
static double tstar = 30.0;   /* per-run read-out time for flash fraction / surviving fraction */   /* pixel lattice: period along helix, row spacing (nm) */

/* ---------- lattice ---------- */
static int     N;
static double *X, *Y;
static int    *nstart, *nbr;
static double *nrate;
static char   *present;
static double  Wx, Wy, dxs, dys;

static void ideal_pos(int i, int j, double *x, double *y) {
    if (lat == 0)      { *x = i;                    *y = 0; }
    else if (lat == 1) { *x = i;                    *y = j; }
    else if (lat == 2) { *x = i + 0.5 * (j & 1);    *y = j * sqrt(3.0) / 2.0; }
    else {
        double nn = sqrt(0.25 * px_P * px_P + px_b * px_b);
        *x = (i + 0.5 * (j & 1)) * px_P / nn;
        *y = j * px_b / nn;
    }
}

static void build_lattice(void) {
    static int allocated = 0;
    N = Lx * Ly;
    double x1, y1, x0, y0;
    ideal_pos(0, 0, &x0, &y0);
    ideal_pos(1, 0, &x1, &y1);
    dxs = x1 - x0;
    if (Ly > 1) { double xa, ya, xb, yb; ideal_pos(0, 0, &xa, &ya); ideal_pos(0, 2, &xb, &yb); dys = (yb - ya) / 2.0; }
    else dys = 1.0;
    Wx = Lx * dxs; Wy = Ly * dys;
    int mi = (int)ceil((cut + 6 * sigma) / dxs) + 1;
    int mj = (Ly > 1) ? (int)ceil((cut + 6 * sigma) / dys) + 1 : 0;
    if (periodic) {
        if (2 * mi + 1 > Lx) { fprintf(stderr, "Lx too small for cutoff\n"); exit(1); }
        if (Ly > 1 && 2 * mj + 1 > Ly) { fprintf(stderr, "Ly too small for cutoff\n"); exit(1); }
        if ((lat == 2 || lat == 3) && (Ly & 1)) { fprintf(stderr, "periodic staggered lattice needs even Ly\n"); exit(1); }
    }
    long maxn = (long)N * (2 * mi + 1) * (2 * mj + 1);
    if (!allocated) {
        X = malloc(sizeof(double) * N); Y = malloc(sizeof(double) * N);
        present = malloc(N);
        nstart = malloc(sizeof(int) * (N + 1));
        nbr = malloc(sizeof(int) * maxn);
        nrate = malloc(sizeof(double) * maxn);
        allocated = 1;
    }
    for (int j = 0; j < Ly; j++) for (int i = 0; i < Lx; i++) {
        int id = j * Lx + i;
        double x, y; ideal_pos(i, j, &x, &y);
        X[id] = x + sigma * gauss();
        Y[id] = y + ((Ly > 1 || sigma > 0) ? sigma * gauss() : 0.0);
        present[id] = (urand() < occ) ? 1 : 0;
    }
    long c = 0;
    for (int j = 0; j < Ly; j++) for (int i = 0; i < Lx; i++) {
        int id = j * Lx + i;
        nstart[id] = (int)c;
        if (!present[id]) continue;
        for (int dj = -mj; dj <= mj; dj++) for (int di = -mi; di <= mi; di++) {
            if (di == 0 && dj == 0) continue;
            int i2 = i + di, j2 = j + dj;
            if (periodic) { i2 = ((i2 % Lx) + Lx) % Lx; j2 = ((j2 % Ly) + Ly) % Ly; }
            else if (i2 < 0 || i2 >= Lx || j2 < 0 || j2 >= Ly) continue;
            int id2 = j2 * Lx + i2;
            if (!present[id2]) continue;
            /* ideal separation decides membership of the neighbour list; jittered one sets the rate */
            double xa, ya, xb, yb;
            ideal_pos(i, j, &xa, &ya); ideal_pos(i + di, j + dj, &xb, &yb);
            double dix = xb - xa, diy = yb - ya;
            double dideal = sqrt(dix * dix + diy * diy);
            if (dideal > cut) continue;
            /* jitter of each end = actual - ideal position (ideal taken at the wrapped index) */
            double xj_id, yj_id;
            ideal_pos(i2, j2, &xj_id, &yj_id);
            double jx = (X[id2] - xj_id) - (X[id] - xa);
            double jy = (Y[id2] - yj_id) - (Y[id] - ya);
            double ddx = dix + jx, ddy = diy + jy;
            double d = sqrt(ddx * ddx + ddy * ddy);
            if (d < 0.25) d = 0.25;
            nbr[c] = id2;
            nrate[c] = r1 / (d * d * d * d * d * d);
            c++;
        }
    }
    nstart[N] = (int)c;   /* absent sites keep an empty neighbour range */
}

/* ---------- Fenwick tree ---------- */
static double *tree, *q;
static double Q;
static int    topbit;
static inline void fen_add(int i, double d) { for (i++; i <= N; i += i & (-i)) tree[i] += d; }
static inline int fen_find(double u) {
    int pos = 0;
    for (int step = topbit; step; step >>= 1) {
        int np = pos + step;
        if (np <= N && tree[np] <= u) { pos = np; u -= tree[np]; }
    }
    return pos;   /* 0-indexed site */
}
static void fen_rebuild(void) {
    memset(tree, 0, sizeof(double) * (N + 1));
    Q = 0;
    for (int i = 0; i < N; i++) { tree[i + 1] += q[i]; Q += q[i]; int p = (i + 1) + ((i + 1) & -(i + 1)); if (p <= N) tree[p] += tree[i + 1]; }
}
static inline void setq(int i, double nq) { double d = nq - q[i]; if (d != 0.0) { q[i] = nq; fen_add(i, d); Q += d; } }

int main(int argc, char **argv) {
    for (int a = 1; a < argc; a++) {
        char *eq = strchr(argv[a], '='); if (!eq) continue; *eq = 0; char *k = argv[a], *v = eq + 1;
        if (!strcmp(k, "lat")) lat = atoi(v); else if (!strcmp(k, "Lx")) Lx = atoi(v); else if (!strcmp(k, "Ly")) Ly = atoi(v);
        else if (!strcmp(k, "periodic")) periodic = atoi(v); else if (!strcmp(k, "cut")) cut = atof(v);
        else if (!strcmp(k, "r1")) r1 = atof(v); else if (!strcmp(k, "h")) hfield = atof(v); else if (!strcmp(k, "pb")) pb = atof(v);
        else if (!strcmp(k, "occ")) occ = atof(v); else if (!strcmp(k, "sigma")) sigma = atof(v);
        else if (!strcmp(k, "tmax")) tmax = atof(v); else if (!strcmp(k, "tmin")) tmin = atof(v); else if (!strcmp(k, "K")) K = atoi(v);
        else if (!strcmp(k, "runs")) runs = atoi(v); else if (!strcmp(k, "mode")) mode = atoi(v);
        else if (!strcmp(k, "seed")) seed0 = strtoull(v, 0, 10); else if (!strcmp(k, "out")) strncpy(outp, v, 500);
        else if (!strcmp(k, "fdead")) fdead_max = atof(v); else if (!strcmp(k, "maxbursts")) maxbursts = atol(v);
        else if (!strcmp(k, "newlat")) newlat = atoi(v);
        else if (!strcmp(k, "tstar")) tstar = atof(v); else if (!strcmp(k, "seedrand")) seedrand = atoi(v);
        else if (!strcmp(k, "P")) px_P = atof(v); else if (!strcmp(k, "b")) px_b = atof(v);
        else { fprintf(stderr, "unknown key %s\n", k); return 1; }
    }
    if (lat == 0) Ly = 1;
    seed_rng(seed0);
    N = Lx * Ly;
    topbit = 1; while ((topbit << 1) <= N) topbit <<= 1;
    tree = calloc(N + 1, sizeof(double)); q = calloc(N, sizeof(double));
    char *st = malloc(N); double *w = malloc(sizeof(double) * N); int *non = malloc(sizeof(int) * N);
    double *tg = malloc(sizeof(double) * K);
    for (int k = 0; k < K; k++) tg[k] = tmin * pow(tmax / tmin, (double)k / (K - 1));
    double *s_rho = calloc(K, sizeof(double)), *s_rho2 = calloc(K, sizeof(double)), *s_surv = calloc(K, sizeof(double));
    double *s_R2 = calloc(K, sizeof(double)), *s_N = calloc(K, sizeof(double)), *s_cum = calloc(K, sizeof(double));
    double *s_alive = calloc(K, sizeof(double));
    double *cum = malloc(sizeof(double) * K), *alv = malloc(sizeof(double) * K), *onf = malloc(sizeof(double) * K);
    char fn[600];
    FILE *frun = NULL, *fb = NULL;
    snprintf(fn, 600, "%s.runs", outp); frun = fopen(fn, "w");
    fprintf(frun, "# Np Ctot_per_site periods_per_site T10 T50 T80 T90 finished tend peakfrac flashfrac_tstar alive_tstar on_tstar\n");
    if (mode == 2) { snprintf(fn, 600, "%s.bursts", outp); fb = fopen(fn, "wb"); }
    long nb_total = 0;
    long long events_total = 0;

    for (int run = 0; run < runs; run++) {
        if (run == 0 || newlat) build_lattice();
        int Np = 0; for (int i = 0; i < N; i++) Np += present[i];
        if (Np == 0) { run--; continue; }
        int Non = 0, ndead = 0; long nper = 0;
        double x0 = 0, y0 = 0; int origin = -1;
        for (int i = 0; i < N; i++) { w[i] = 0; non[i] = 0; q[i] = 0; st[i] = present[i] ? 0 : 2; }
        if (mode == 0) { for (int i = 0; i < N; i++) if (present[i]) { st[i] = 1; Non++; } }
        else if (mode == 1) {
            origin = (Ly / 2) * Lx + Lx / 2;
            if (seedrand) { do { origin = (int)(urand() * N); } while (origin >= N || !present[origin]); }
            if (!present[origin]) { /* pick nearest present site along the row */
                int o2 = origin; while (o2 < N && !present[o2]) o2++; if (o2 >= N) { run--; continue; } origin = o2; }
            st[origin] = 1; Non = 1;
        }
        if (origin >= 0) { x0 = X[origin]; y0 = Y[origin]; }
        for (int i = 0; i < N; i++) if (st[i] == 1) for (int c = nstart[i]; c < nstart[i + 1]; c++) { w[nbr[c]] += nrate[c]; non[nbr[c]]++; }
        for (int i = 0; i < N; i++) q[i] = (st[i] == 1) ? 1.0 : (st[i] == 0 ? hfield + w[i] : 0.0);
        fen_rebuild();
        double SR2 = 0;   /* sum over ON sites of squared distance from origin (seed mode) */
        double t = 0, C = 0; int kg = 0; long ev = 0; int peak = Non;
        /* avalanche bookkeeping */
        int inburst = 0; double bS = 0, bT0 = 0; long bN = 0; double bdead0 = 0; int bpeak = 0;
        int finished = 0; double tend = tmax;
        while (1) {
            if (ndead == Np || (Non == 0 && hfield <= 0.0) || Q <= 1e-12) {   /* all bleached, or absorbing */
                finished = 1; tend = t;
                while (kg < K) { cum[kg] = C; alv[kg] = (double)(Np - ndead) / Np; onf[kg] = 0; s_alive[kg] += alv[kg]; kg++; }
                break;
            }
            double dt = -log(1.0 - urand()) / Q;
            double tn = t + dt;
            while (kg < K && tg[kg] <= tn) {
                double rho = (double)Non / Np;
                s_rho[kg] += rho; s_rho2[kg] += rho * rho; s_N[kg] += Non; s_surv[kg] += (Non > 0);
                s_R2[kg] += SR2; s_alive[kg] += (double)(Np - ndead) / Np;
                cum[kg] = C + Non * (tg[kg] - t); alv[kg] = (double)(Np - ndead) / Np; onf[kg] = rho;
                kg++;
            }
            if (tn >= tmax) { C += Non * (tmax - t); t = tmax; break; }
            C += Non * dt; if (inburst) bS += Non * dt;
            t = tn;
            int i = fen_find(urand() * Q);
            if (i >= N) i = N - 1;
            if (q[i] <= 0.0) { fen_rebuild(); continue; }
            ev++;
            if (st[i] == 1) {            /* ON -> OFF or bleached */
                Non--; nper++;
                if (origin >= 0) { double ddx = X[i] - x0, ddy = Y[i] - y0; SR2 -= ddx * ddx + ddy * ddy; }
                for (int c = nstart[i]; c < nstart[i + 1]; c++) {
                    int j = nbr[c]; if (st[j] == 2) continue;
                    non[j]--; if (non[j] == 0) w[j] = 0; else w[j] -= nrate[c];
                    if (st[j] == 0) setq(j, hfield + w[j]);
                }
                if (pb > 0 && urand() < pb) { st[i] = 2; ndead++; setq(i, 0.0); }
                else { st[i] = 0; setq(i, hfield + w[i]); }
                if (inburst) { bN++; if (Non == 0) {
                    inburst = 0;
                    if (fb && bdead0 <= fdead_max && nb_total < maxbursts) { double rec[5] = {bS, t - bT0, (double)bN, bdead0, (double)bpeak}; fwrite(rec, sizeof(double), 5, fb); nb_total++; }
                } }
            } else {                     /* OFF -> ON */
                st[i] = 1; Non++; if (Non > peak) peak = Non;
                setq(i, 1.0);
                if (origin >= 0) { double ddx = X[i] - x0, ddy = Y[i] - y0; SR2 += ddx * ddx + ddy * ddy; }
                for (int c = nstart[i]; c < nstart[i + 1]; c++) {
                    int j = nbr[c]; if (st[j] == 2) continue;
                    non[j]++; w[j] += nrate[c];
                    if (st[j] == 0) setq(j, hfield + w[j]);
                }
                if (mode == 2) { if (!inburst) { inburst = 1; bS = 0; bT0 = t; bN = 0; bdead0 = (double)ndead / Np; bpeak = 1; } else if (Non > bpeak) bpeak = Non; }
            }
            if ((ev & ((1L << 22) - 1)) == 0) fen_rebuild();
            if (mode == 2 && ((double)ndead / Np > fdead_max || nb_total >= maxbursts) && !inburst) { tend = t; break; }
        }
        events_total += ev;
        /* cumulative curve and burn times */
        for (int k = kg; k < K; k++) { cum[k] = C; alv[k] = (double)(Np - ndead) / Np; onf[k] = (double)Non / Np; }
        int ks = 0; while (ks < K - 1 && tg[ks] < tstar) ks++;
        for (int k = 0; k < K; k++) s_cum[k] += cum[k] / Np;
        double Tf[4] = {-1, -1, -1, -1}, fr[4] = {0.1, 0.5, 0.8, 0.9};
        if (C > 0) for (int f = 0; f < 4; f++) {
            double target = fr[f] * C;
            for (int k = 0; k < K; k++) if (cum[k] >= target) {
                if (k == 0) Tf[f] = tg[0] * target / cum[0];
                else Tf[f] = tg[k - 1] + (tg[k] - tg[k - 1]) * (target - cum[k - 1]) / (cum[k] - cum[k - 1] + 1e-300);
                break;
            }
        }
        fprintf(frun, "%d %.6g %.6g %.6g %.6g %.6g %.6g %d %.6g %.4f %.5f %.5f %.5f\n", Np, C / Np, (double)nper / Np, Tf[0], Tf[1], Tf[2], Tf[3], finished, tend, (double)peak / Np,
                C > 0 ? cum[ks] / C : 0.0, alv[ks], onf[ks]);
    }
    fclose(frun); if (fb) fclose(fb);
    snprintf(fn, 600, "%s.ts", outp); FILE *f = fopen(fn, "w");
    fprintf(f, "# t rho rho2 surv N R2 cum_per_site alive   (runs=%d N=%d events=%lld)\n", runs, N, events_total);
    for (int k = 0; k < K; k++) {
        double Nm = s_N[k] / runs;
        fprintf(f, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", tg[k], s_rho[k] / runs, s_rho2[k] / runs, s_surv[k] / runs, Nm,
                (s_N[k] > 0 ? s_R2[k] / s_N[k] : 0.0), s_cum[k] / runs, s_alive[k] / runs);
    }
    fclose(f);
    return 0;
}
