#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double mean(double* x, int N) {
    int i;
    double m=0.0;

    for (i=0; i<N; i++) {
        m = (x[i] + (i + 1) * m) / (i + 2);
    }

    return(m);
}

int binomial(double p) {
    double r = (float) rand() / (float) RAND_MAX;
    return(r < p);
}

double sim(double s, double T, double f, double g, double tg, int N, int n_generations) {
    int generation_i, host_i;
    double p_treatment, E, E_new;
    double* H;

    p_treatment = T * tg; // p_treatment * 1 year / (tx per year) = T

    E = 0.5;
    H = (double *) malloc(N);
    for (host_i=0; host_i < N; host_i++) H[host_i] = 0.91;

    for (generation_i=0; generation_i < n_generations; generation_i++) {
        // compute new environment
        // environment gets host bugs and experiences in-environment selection
        E_new = E + f * (mean(H, N) - E) - s * E * (1 - E) / (1.0 - s * E);

        // compute & update hosts
        for (host_i=0; host_i < N; host_i++) {
            // does this host get treatment?
            if (binomial(p_treatment)) {
                // treated hosts get 100% resistance
                H[host_i] = 1.0;
            } else {
                // untreated hosts get environmental bugs and experience in-host selection
                H[host_i] = H[host_i] + g * (E - H[host_i]) - s * H[host_i] * (1.0 - H[host_i]) / (1.0 - s * H[host_i]);
            }
        }

        // update environment
        E = E_new;
    }

    free(H);

    return(mean(H, N));
}

int main(void) {
    srand(time(NULL));
    printf("%f\n", sim(0.01, 2.0, 0.05, 0.005, 1.0/219, 10000, 10000));
    return(0);
}
