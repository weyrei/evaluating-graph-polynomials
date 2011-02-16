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
#define MAXN  100 //max # of vertices
int n;
int adj[MAXN][MAXN];
int pow2; //2 to the power of n
int *p; //p value from bhkk-paper
int *c; //c value from bhkk-paper
int**** w;

/*
 * 
 */

int getNumOfNodes(int U) {
    int n = 0;
    int S = U;
    while (S > 0) {
        n += (S & 1);
        S=S >> 1;
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
                else
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
    int S, j;
    int**a = (int**) malloc(pow2 * sizeof (int*));
    for (S = 0; S < pow2; S++) {
        a[S] = (int*) malloc((n + 1) * sizeof (int));
        a[S][0] = w[S][minElem(S)][minElem(S)][getNumOfNodes(S)];
        for (j = 1; j < n; j++) {
            if ((1 << j) & S)
                a[S][j] = a[S][j - 1] - a[S^(1 << j)][j - 1];
            else
                a[S][j] = a[S][j - 1];
        }
    }
    for (S = 0; S < pow2; S++) {
        c[S] = a[S][n];
        free(a[S]);
    }
    free(a);
    printf("succesfully cumputed C\n");
}

void computeP() {       //fast moebius inversion
    printf("started computing P\n");
    int S, j, s, t;
    int sum = 0;
    int**a = (int**) malloc(pow2 * sizeof (int*));
    for (S = 0; S < pow2; S++) {
        a[S] = (int*) malloc((n + 1) * sizeof (int));
        for (t = 0; t < n; t++) {
            for (s = 0; s <= t; s++) {
                a[S][0] = w[S][s][t][getNumOfNodes(S) - 1];
                for (j = 1; j < n; j++) {
                    if ((1 << j) & S)
                        a[S][j] = a[S][j - 1] - a[S^(1 << j)][j - 1]; //f[j](X) = - f[j-1](X\{j}) + f[j-1](X)
                    else
                        a[S][j] = a[S][j - 1];  //f[j](X) =  f[j-1](X)
                }
                sum += a[S][n];
            }
        }
        p[S] = sum;
        sum=0;
        free(a[S]);
    }
    free(a);
    printf("succesfully computed P\n");
}

void freeW() {
    int S, i, j, l;
    for (S = 0; S < pow2; S++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                free(w[S][i][j]);
            }
            free(w[S][i]);
        }
        free(w[S]);
    }
    free(w);
}

int main(int argc, char** argv) {
    if(argc<2){
        fprintf(stderr,"no file specified\n");
        exit(EXIT_FAILURE);
    }
    int opt = -1;
    char* filename = (char*) malloc(20 * sizeof (char));
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        printf("opt: %c\n", opt);
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
    printf("%d\n",w[pow2-1][3][10][2]);
    printf("succesfully computed W\n");
    computeC();
    computeP();

    freeW();
    free(c);
    free(p);
    free(filename);
    return (EXIT_SUCCESS);
}



