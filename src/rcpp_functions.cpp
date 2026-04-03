#include <Rcpp.h>
#include <set>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector inbreeding_meta(NumericMatrix ped, NumericMatrix gamma) {

    const int n = ped.nrow();
    const int nmf = gamma.nrow();

    // Extracting parents
    NumericMatrix::Column sire = ped(_, 1);
    NumericMatrix::Column dam  = ped(_, 2);

    //-------------------------
    // Output vector of inbreeding F[i] = Aii - 1
    //-------------------------
    NumericVector F(n, NA_REAL);

    //--------------------------
    // Initialize metafounders
    // For metafounder k:
    // A_kk = gamma_kk => F_k = gamma_kk - 1
    //-------------------------
    for (int k = 0; k < nmf; ++k) {
        F[k] = gamma(k, k) - 1.0;
    }

    //---------------------------
    // Processing non-metafounders
    //---------------------------
    for (int i = nmf; i < n; ++i) {

        // Lrow stores the recursive contributions from i to its ancestors
        NumericVector Lrow(n, 0.0);
        std::set<int> active;

        // Lmf stores the total contributions that reach each metafounder
        NumericVector Lmf(nmf, 0.0);

        // Start from individual i itself
        Lrow[i] = 1.0;
        active.insert(i);

        double Aii = 0.0;

        while (!active.empty()) {

            int j = *active.rbegin();
            active.erase(j);

            double Lj = Lrow[j];
            if (Lj == 0.0) continue;
            Lrow[j] = 0.0;

            // Case 1: j is a metafounder
            if (j < nmf) {
                Lmf[j] += Lj;
                continue;
            }

            // Case 2: j is a regular individual
            int sj = static_cast<int>(sire[j]);
            int dj = static_cast<int>(dam[j]);

            // If parent = 0, use F = -1 so that unknown unrelated parents
            double F_s = (sj == 0) ? -1.0 : F[sj - 1];
            double F_d = (dj == 0) ? -1.0 : F[dj - 1];

            double Dj = 0.5 - 0.25 * (F_s + F_d);

            // Add contributions from individual j
            Aii += Lj * Lj * Dj;

            // Propagate to parents
            if (sj != 0) {
                Lrow[sj - 1] += 0.5 * Lj;
                active.insert(sj - 1);
            }
            if (dj != 0) {
                Lrow[dj - 1] += 0.5 * Lj;
                active.insert(dj - 1);
            }
        }

        //------------------------
        // Legarra correction for MF block
        // contribution = Lmf' * gamma * Lmf
        //------------------------
        double quad_form = 0.0;
        for (int r = 0; r < nmf; ++r) {
            for (int c = 0; c < nmf; ++c) {
                quad_form += Lmf[r] * gamma(r, c) * Lmf[c];
            }
        }

        Aii += quad_form;
        F[i] = Aii - 1.0;
    }

    return F;
}
