/* 
 * File:   main.c
 * Author: Reinhard Weyand
 *
 * Implementation of the fast evaluation algorithm for the coefficients of the cover polynomial
 * as presented in "Computing The Tutte Polynomial In Vertex-Exponential Time" by Husfeldt et. al. (Appendix p.3-4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <unistd.h>

#define MAXL 10 //max length of directed paths
#define MAXN  5 //max # of vertices
typedef int* polynom;
int n;
int adj[MAXN][MAXN];
int pow2; //2 to the power of n
int *p; //p value from bhkk-paper
int *c; //c value from bhkk-paper
int** P_U; //size at first position
int** C_U;
int P_U_I[MAXN*MAXN + 1]; //P to the power of i
int C_U_J[MAXN*MAXN + 1]; //C to the power of J
int PC_IJ[2 * MAXN*MAXN*MAXN + 1]; //(P^i)*(C*j)
int C_D[MAXN][MAXN];
int**** w;

/*
 * 
 */

void polyAdd(polynom a, polynom b, polynom r) {
    int i;
    int sizea = a[0];
    int sizeb = b[0];
    polynom x = sizea < sizeb ? b : a;
    polynom y = sizea <= sizeb ? a : b;
    sizea = x[0]; //bigger size
    sizeb = y[0]; //smaller size
    r[0] = sizea > 0 ? sizea : 0;
    for (i = 1; i <= sizea; i++) {
        r[i] = i <= sizeb ? x[i] + y[i] : x[i];
    }
}

void polyMul(int* a, int* b, int* r) {
    int i, j, sizer;
    int sizea = a[0];
    int sizeb = b[0];
    if (sizea == 0) {
        printf("polynom a empty!\n");
        r[0] = sizeb;
        for (i = 1; i <= sizeb; i++)
            r[i] = b[i];
        return;
    } else if (sizeb == 0) {
        printf("polynom b empty!\n");
        r[0] = sizea;
        for (i = 1; i <= sizea; i++)
            r[i] = a[i];
        return;
    }
    //  printf("PolyMUL  ****  sizes of polynoms: %d  and %d\n", sizea, sizeb);
    for (i = 0; i <= sizea + sizeb - 1; i++) {
        //   printf("initializing r at %d\n",i);
        r[i] = 0;
    }
    //   printf("initialized r with 0\n");
    for (i = 1; i <= sizea; i++) {
        for (j = 1; j <= sizeb; j++) {
            r[i + j - 1] += a[i] * b[j];
        }
    }
    sizer = (sizea + sizeb - 1);
    r[0] = sizer > 0 ? sizer : 0;

}

void polyPow(int* a, int exp, int* r) {
    //    printf("computing a^%d\n", exp);
    //    printf("size of polynom a is %d\n", a[0]);
    int i, j, size;
    r[0] = 0;
    /* for (i = 0; i <= a[0] * exp; i++) {
         r[i] = 0;
     }*/
    int buf[a[0] * exp + 1];
    for (i = 0; i <= a[0]; i++)
        buf[i] = a[i];
    for (i = 1; i < exp; i++) {
        polyMul(buf, a, r);
        size = r[0];
        for (j = 0; j <= size; j++) {
            buf[j] = r[j];
            r[j] = 0;
        }

    }
    for (i = 0; i <= buf[0]; i++) {
        r[i] = buf[i];
    }
}

int getNumOfNodes(int U) {
    int n = 0;
    int S = U;
    while (S > 0) {
        n += (S & 1);
        S = S >> 1;
    }
    return n;
}

int minElem(int X) {
    int ret = 0;
    while (!X & ret)
        X << 1;
    return ret;
}

int mypow(int X) {
    int i, ret;
    ret = 1;
    for (i = 0; i < X; i++)
        ret *= 2;
    return ret;
}

int fac(int x) {
    int i;
    int fac = 1;
    for (i = 1; i <= x; i++)
        fac *= i;
    return fac;
}

void computeWValue(int s) {//approach through dynamic programming. a little hard, i think

    int i, j, l, v;
    /*   int*** a = (int***) malloc(n*sizeof(int**));
     *a = (int**) malloc(n*n*sizeof(int*));
     **a = (int*) malloc(n*n*n*sizeof(int)); */
    w[s] = (int***) malloc(n * sizeof (int**));
    for (i = 0; i < n; i++) {
        w[s][i] = (int**) malloc(n * sizeof (int*));
        for (j = 0; j < n; j++) {
            w[s][i][j] = (int*) malloc(n * sizeof (int));
            if (((1 << j) & s) && ((1 << i) & s)) {
                if (i == j)
                    w[s][i][j][0] = 1;

                w[s][i][j][1] = adj[i][j] ? 1 : 0;
            }
        }

    }
    int count;
    for (l = 2; l < n; l++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (((1 << j) & s) && (1 << i) & s) { //both nodes in subgraph?
                    count = 0;
                    for (v = 0; v < n; v++) {
                        if ((1 << v) & s) {
                            if (adj[v][j])
                                count = count + w[s][i][v][l - 1]; //reccurence
                        }
                    }
                    w[s][i][j][l] = count;
                }
            }
        }

    }





}

void parseGraph(char* filename) {
    FILE* f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    int e;
    int i = 0;
    int j = 0;
    fscanf(f, "%d", &n);
    printf("got graph with %d nodes\n", n);
    while (!feof(f)) {
        fscanf(f, "%d", &e);
        if (!isspace(e)) {
            if (e < 0 && e > 1) {
                fprintf(stderr, "value different from 1 or 0 in adjacency matrix\n");
                exit(EXIT_FAILURE);
            } else {
                //            printf("parsing indice: %d, %d, val: %d, %d nodes\n", j, i, e,n);
                adj[j][i] = e;
                j = j + (i == (n - 1) ? 1 : 0);
                i = (i + 1) % n;
            }
        }
    }

}

void computeC() { //fast moebius inversion
    printf("starting to compute C\n");
    int S, j, s;
    int sum = 0;
    int**a = (int**) malloc(pow2 * sizeof (int*));
    for (S = 0; S < pow2; S++) {
        for (s = 0; s < n; s++) {
            a[S] = (int*) malloc((n + 1) * sizeof (int));
            a[S][0] = w[S][s][s][getNumOfNodes(S)];

            for (j = 1; j < n; j++) {
                if ((1 << j) & S)
                    a[S][j] = a[S][j - 1] - a[S^(1 << j)][j - 1];
                else
                    a[S][j] = a[S][j - 1];
            }
            sum += a[S][n - 1];

        }
        c[S] = sum;
        sum = 0;
        free(a[S]);
    }
    free(a);
    printf("succesfully cumputed C\n");
}

void computeP() { //fast moebius inversion
    printf("started computing P\n");
    int S, j, s, t;
    int sum = 0;
    int**a = (int**) malloc(pow2 * sizeof (int*));
    for (S = 0; S < pow2; S++) {
        a[S] = (int*) malloc((n + 1) * sizeof (int));
        for (t = 0; t < n; t++) {
            for (s = 0; s < n; s++) {
                a[S][0] = (S != 0) ? w[S][s][t][getNumOfNodes(S) - 1] : w[S][s][t][0];
                //   printf("a[%d][0] = %d\n",S,a[S][0]);
                for (j = 1; j < n; j++) {
                    if ((1 << j) & S)
                        a[S][j] = a[S][j - 1] - a[S^(1 << j)][j - 1]; //f[j](X) = - f[j-1](X\{j}) + f[j-1](X)
                    else
                        a[S][j] = a[S][j - 1]; //f[j](X) =  f[j-1](X)
                    //printf("a[%d][%d] = %d\n",S,j,a[S][j]);
                }
                sum += a[S][n - 1];

            }
        }
        p[S] = sum;
        //   printf("computed sum %d\n", sum);
        sum = 0;
        free(a[S]);
    }
    free(a);
    printf("succesfully computed P\n");
}

void computeP_U() {
    P_U = (int**) malloc(pow2 * sizeof (int));
    int U, S, i, size;
    for (U = 0; U < pow2; U++) {
        P_U[U] = (int*) malloc((MAXN + 1) * sizeof (int));
        for (S = 0; S <= U; S++) {
            size = getNumOfNodes(S);
            P_U[U][size + 1] += p[S];
            //  printf("P[S] = %d\n",p[S]);
        }
        size = MAXN;
        for (i = MAXN; i > 0; i--) {
            if (P_U[U][i] == 0)
                size--;
            else break;
        }
        P_U[U][0] = size;
        //    printf("size of P_U[%d] = %d\n",U,P_U[U][0]);
    }
}

void computeC_U() {
    C_U = (int**) malloc(pow2 * sizeof (int));
    int U, S, i, size;
    for (U = 0; U < pow2; U++) {
        C_U[U] = (int*) malloc((MAXN + 1) * sizeof (int));
        for (S = 0; S <= U; S++) {
            size = getNumOfNodes(S);
            C_U[U][size + 1] += c[S];
        }
        size = MAXN;
        for (i = MAXN; i > 0; i--) {
            if (C_U[U][i] == 0)
                size--;
            else break;
        }
    }
}

void freeAll() {
    int S, i, j, l;
    for (S = 0; S < pow2; S++) {
        free(C_U[S]);
        free(P_U[S]);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                free(w[S][i][j]);
            }
            free(w[S][i]);
        }
        free(w[S]);
    }
    free(w);
    free(C_U);
    free(P_U);
    free(c);
    free(p);
}

void computeC_D(int i, int j) {
    printf("********************************* STARTED COMPUTING C*********************************\n");
    C_D[i][j] = 0;
    int S, vorz;
    int ifac = fac(i);
    int jfac = fac(j);
    for (S = 0; S < pow2; S++) {
        printf("at set %d\n", S);
        vorz = S & 1 ? 1 : -1;
        polyPow(P_U[S], i, P_U_I);
        polyPow(C_U[S], j, C_U_J);
        polyMul(P_U_I, C_U_J, PC_IJ);
        int t;
        C_D[i][j] += vorz * PC_IJ[n + 1];
        printf("sum = %d\n", C_D[i][j]);
        printf("Polynom P_S[%d] = %d\n", t, PC_IJ[n + 1]);
    }

    printf("!!!!!!!!!      C_D[%d][%d] = %d   !!!!!!!!!!!\n", i, j, C_D[i][j]);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "no file specified\n");
        exit(EXIT_FAILURE);
    }
    int opt = -1;
    char* filename = (char*) malloc(20 * sizeof (char));
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        //   printf("opt: %c\n", opt);
        switch (opt) {
            case 'f':
                strcpy(filename, optarg);
                parseGraph(filename);
                break;
        }
    }
    pow2 = mypow(n);
    w = (int****) malloc(pow2 * sizeof (int***));
    c = (int*) malloc(pow2 * sizeof (int));
    p = (int*) malloc(pow2 * sizeof (int));
    int S;
    printf("starting to compute W\n");
    for (S = 0; S < pow2; S++) {
        computeWValue(S);
    }
    // printf("%d\n", w[3][0][1][1]);
    printf("succesfully computed W\n");
    computeC();
    computeP();
    computeP_U();
    computeC_U();
    computeC_D(1, 1);
    free(filename);
    freeAll();
    return (EXIT_SUCCESS);
}



