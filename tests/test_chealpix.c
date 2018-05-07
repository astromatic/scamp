/**
 * test chealpix
 *
 * Here we are only testing "neighbors_nest" because this is the one we had
 * to implement from the C++ code.
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../src/chealpix.h"

int main(int argc, char **argv)
{
    long i;
    int j, k;
    int64_t nsides = 16;
    int64_t neigh[8], testneigh[8];

    FILE *fd = fopen("extra/test_chealpix.data", "r");

    for (i=0; i<nside2npix(nsides); i++) {
        neighbours_nest64(nsides, i, neigh);
        k = fscanf(fd, "%li:%li:%li:%li:%li:%li:%li:%li:\n",
                &testneigh[0],&testneigh[1],&testneigh[2],&testneigh[3],
                &testneigh[4],&testneigh[5],&testneigh[6],&testneigh[7]);
        assert(k == 8);

        for (j=0; j<8; j++)
            assert(neigh[j] == testneigh[j]);

    }

    return 0;
}
