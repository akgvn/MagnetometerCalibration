// This implementation makes use of the C-language matrix source code available at
//
// http://www.mymathlib.com/matrices/
//
// It has been tested as a console application, and reproduces exactly the results given
// by the FORTRAN implementation of Magneto 1.2, using the test file mag.txt with a user
// norm of 0.569.
//
//
// Copyright (C) 2013 www.sailboatinstruments.blogspot.com
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this
// software and associated documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
// to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <float.h>

#include "math_functions.h"


// "It has been tested as a console application, and reproduces exactly
// the results given by the FORTRAN implementation of Magneto 1.2, using
// the test file mag.txt with a user norm of 0.569."


double* read_data_from_file(int* number_of_lines) {
    int nlines = 0;
    char buf[120];

    FILE* fp = fopen("mag.txt", "r");
    while (fgets(buf, 100, fp) != NULL)
        nlines++;

    rewind(fp);

    double* D = (double*) malloc(10 * nlines * sizeof(double));

    for (int i = 0; i < nlines; i++) {
        fgets(buf, 100, fp);

        double x, y, z;
        sscanf(buf, "%lf\t%lf\t%lf", &x, &y, &z);

        D[i] = x * x;
        D[nlines + i] = y * y;
        D[nlines * 2 + i] = z * z;
        D[nlines * 3 + i] = 2.0 * y * z;
        D[nlines * 4 + i] = 2.0 * x * z;
        D[nlines * 5 + i] = 2.0 * x * y;
        D[nlines * 6 + i] = 2.0 * x;
        D[nlines * 7 + i] = 2.0 * y;
        D[nlines * 8 + i] = 2.0 * z;
        D[nlines * 9 + i] = 1.0;
    }
    fclose(fp);

    *number_of_lines = nlines;
    return D;
}

int main() {
    int nlines = 0;
    double* D = read_data_from_file(&nlines);

    double S[10 * 10];

    Matrix_x_Its_Transpose(S, D, 10, nlines);
    free(D);

    // Create pre-inverted constraint matrix C

    double C[6 * 6];

    C[0] = 0.0;
    C[1] = 0.5;
    C[2] = 0.5;
    C[3] = 0.0;
    C[4] = 0.0;
    C[5] = 0.0;

    C[6] = 0.5;
    C[7] = 0.0;
    C[8] = 0.5;
    C[9] = 0.0;
    C[10] = 0.0;
    C[11] = 0.0;

    C[12] = 0.5;
    C[13] = 0.5;
    C[14] = 0.0;
    C[15] = 0.0;
    C[16] = 0.0;
    C[17] = 0.0;

    C[18] = 0.0;
    C[19] = 0.0;
    C[20] = 0.0;
    C[21] = -0.25;
    C[22] = 0.0;
    C[23] = 0.0;

    C[24] = 0.0;
    C[25] = 0.0;
    C[26] = 0.0;
    C[27] = 0.0;
    C[28] = -0.25;
    C[29] = 0.0;

    C[30] = 0.0;
    C[31] = 0.0;
    C[32] = 0.0;
    C[33] = 0.0;
    C[34] = 0.0;
    C[35] = -0.25;

    double S11[6 * 6];

    Get_Submatrix(S11, 6, 6, S, 10, 0, 0);

    double S12[6 * 4];

    Get_Submatrix(S12, 6, 4, S, 10, 0, 6);

    double S12t[4 * 6];

    Get_Submatrix(S12t, 4, 6, S, 10, 6, 0);

    double S22[4 * 4];

    Get_Submatrix(S22, 4, 4, S, 10, 6, 6);

    double S22_1[4 * 4];

    for (int i = 0; i < 16; i++)
        S22_1[i] = S22[i];

    Choleski_LU_Decomposition(S22_1, 4);

    Choleski_LU_Inverse(S22_1, 4);

    // Calculate S22a = S22_1 * S12t   4*6 = 4x4 * 4x6   C = AB

    double S22a[4 * 6];

    Multiply_Matrices(S22a, S22_1, 4, 4, S12t, 6);

    // Then calculate S22b = S12 * S22a      ( 6x6 = 6x4 * 4x6)

    double S22b[6 * 6];

    Multiply_Matrices(S22b, S12, 6, 4, S22a, 6);

    // Calculate SS = S11 - S22b

    double SS[6 * 6];

    for (int i = 0; i < 36; i++)

        SS[i] = S11[i] - S22b[i];

    double E[6 * 6];

    Multiply_Matrices(E, C, 6, 6, SS, 6);

    double SSS[6 * 6];

    Hessenberg_Form_Elementary(E, SSS, 6);

    double eigen_real[6];

    double eigen_imag[6];

    QR_Hessenberg_Matrix(E, SSS, eigen_real, eigen_imag, 6, 100);

    int index = 0;

    double maxval = eigen_real[0];

    for (int i = 1; i < 6; i++) {
        if (eigen_real[i] > maxval) {
            maxval = eigen_real[i];
            index = i;
        }
    }

    double v1[6];

    v1[0] = SSS[index];
    v1[1] = SSS[index + 6];
    v1[2] = SSS[index + 12];
    v1[3] = SSS[index + 18];
    v1[4] = SSS[index + 24];
    v1[5] = SSS[index + 30];

    // normalize v1
    double norm = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3] + v1[4] * v1[4] + v1[5] * v1[5]);
    v1[0] /= norm;
    v1[1] /= norm;
    v1[2] /= norm;
    v1[3] /= norm;
    v1[4] /= norm;
    v1[5] /= norm;

    if (v1[0] < 0.0) {
        v1[0] = -v1[0];
        v1[1] = -v1[1];
        v1[2] = -v1[2];
        v1[3] = -v1[3];
        v1[4] = -v1[4];
        v1[5] = -v1[5];
    }

    // Calculate v2 = S22a * v1      ( 4x1 = 4x6 * 6x1)
    double v2[4];

    Multiply_Matrices(v2, S22a, 4, 6, v1, 1);
    double v[10];
    v[0] = v1[0];
    v[1] = v1[1];
    v[2] = v1[2];
    v[3] = v1[3];
    v[4] = v1[4];
    v[5] = v1[5];
    v[6] = -v2[0];
    v[7] = -v2[1];
    v[8] = -v2[2];
    v[9] = -v2[3];

    double Q[3 * 3];

    Q[0] = v[0];
    Q[1] = v[5];
    Q[2] = v[4];
    Q[3] = v[5];
    Q[4] = v[1];
    Q[5] = v[3];
    Q[6] = v[4];
    Q[7] = v[3];
    Q[8] = v[2];

    double U[3];

    U[0] = v[6];
    U[1] = v[7];
    U[2] = v[8];

    double Q_1[3 * 3];

    for (int i = 0; i < 9; i++)
        Q_1[i] = Q[i];

    Choleski_LU_Decomposition(Q_1, 3);
    Choleski_LU_Inverse(Q_1, 3);

    // Calculate B = Q-1 * U   ( 3x1 = 3x3 * 3x1)
    double B[3];

    Multiply_Matrices(B, Q_1, 3, 3, U, 1);

    B[0] = -B[0]; // x-axis combined bias
    B[1] = -B[1]; // y-axis combined bias
    B[2] = -B[2]; // z-axis combined bias

    for (int i = 0; i < 3; i++)
        printf("%lf\r\n", B[i]);

    // First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1)

    double QB[3];

    Multiply_Matrices(QB, Q, 3, 3, B, 1);

    // Then calculate btqb = BT * QB    ( 1x1 = 1x3 * 3x1)
    double btqb;
    Multiply_Matrices( &btqb, B, 1, 3, QB, 1);

    // Calculate hmb = sqrt(btqb - J).
    double J = v[9];

    double hmb = sqrt(btqb - J);
    // Calculate SQ, the square root of matrix Q

    double SSSS[3 * 3];

    Hessenberg_Form_Elementary(Q, SSSS, 3);

    double eigen_real3[3];

    double eigen_imag3[3];

    QR_Hessenberg_Matrix(Q, SSSS, eigen_real3, eigen_imag3, 3, 100);

    // normalize eigenvectors

    double norm1 = sqrt(SSSS[0] * SSSS[0] + SSSS[3] * SSSS[3] + SSSS[6] * SSSS[6]);

    SSSS[0] /= norm1;
    SSSS[3] /= norm1;
    SSSS[6] /= norm1;

    double norm2 = sqrt(SSSS[1] * SSSS[1] + SSSS[4] * SSSS[4] + SSSS[7] * SSSS[7]);

    SSSS[1] /= norm2;
    SSSS[4] /= norm2;
    SSSS[7] /= norm2;

    double norm3 = sqrt(SSSS[2] * SSSS[2] + SSSS[5] * SSSS[5] + SSSS[8] * SSSS[8]);

    SSSS[2] /= norm3;
    SSSS[5] /= norm3;
    SSSS[8] /= norm3;

    double Dz[3 * 3];

    for (int i = 0; i < 9; i++)
        Dz[i] = 0.0;

    Dz[0] = sqrt(eigen_real3[0]);
    Dz[4] = sqrt(eigen_real3[1]);
    Dz[8] = sqrt(eigen_real3[2]);
    double vdz[3 * 3];

    Multiply_Matrices(vdz, SSSS, 3, 3, Dz, 3);
    Transpose_Square_Matrix(SSSS, 3);

    double SQ[3 * 3];
    Multiply_Matrices(SQ, vdz, 3, 3, SSSS, 3);

    double hm = 0.569;

    double A_1[3 * 3];

    for (int i = 0; i < 9; i++)
        A_1[i] = SQ[i] * hm / hmb;

    for (int i = 0; i < 3; i++)
        printf("%lf %lf %lf\r\n", A_1[i * 3], A_1[i * 3 + 1], A_1[i * 3 + 2]);

    return 0;
}
