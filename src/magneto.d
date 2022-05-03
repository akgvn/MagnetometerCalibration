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

import core.stdc.stdio;
import core.stdc.math;
import core.stdc.stdlib;
import core.stdc.string;


import math_functions;

Matrix read_data_matrix_from_file(string filename) {
    int nlines = 0;
    char[120] buf;

    import std.string: toStringz;
    FILE* fp = fopen(filename.toStringz(), "r");
    while (fgets(buf.ptr, 100, fp) != null)
        nlines++;

    rewind(fp);

    auto result = Matrix(10, nlines);

    foreach (col; 0..nlines) {
        fgets(buf.ptr, 100, fp);

        double x, y, z;
        sscanf(buf.ptr, "%lf\t%lf\t%lf", &x, &y, &z);

        result.get(0, col) = x * x;
        result.get(1, col) = y * y;
        result.get(2, col) = z * z;
        result.get(3, col) = 2.0 * y * z;
        result.get(4, col) = 2.0 * x * z;
        result.get(5, col) = 2.0 * x * y;
        result.get(6, col) = 2.0 * x;
        result.get(7, col) = 2.0 * y;
        result.get(8, col) = 2.0 * z;
        result.get(9, col) = 1.0;
    }
    fclose(fp);

    return result;
}

struct Result {
    double[3] combined_bias;
    double[3 * 3] correction;
}

Result calculate_the_thing(string filename, double user_norm) {
    const D = read_data_matrix_from_file(filename);
    int numberOfLines = D.cols;

    const S = Matrix_x_Its_Transpose(D);

    const S11  = S.submatrix(6, 6, 0, 0);
    const S12  = S.submatrix(6, 4, 0, 6);
    const S12t = S.submatrix(4, 6, 6, 0);
    const S22  = S.submatrix(4, 4, 6, 6);

    auto S22_1 = Matrix(4, 4); // double[4 * 4] S22_1;

    for (int i = 0; i < 16; i++)
        S22_1.get(i) = S22.get(i);

    Choleski_LU_Decomposition(S22_1.m.ptr, 4);
    Choleski_LU_Inverse(S22_1.m.ptr, 4);

    // Calculate S22a = S22_1 * S12t   4*6 = 4x4 * 4x6   C = AB

    const S22a = Multiply_Matrices(S22_1, S12t);

    // Then calculate S22b = S12 * S22a      ( 6x6 = 6x4 * 4x6)

    const S22b = Multiply_Matrices(S12, S22a);

    // Calculate SS = S11 - S22b

    auto SS = Matrix(6, 6);

    for (int i = 0; i < 36; i++)
        SS.get(i) = S11.get(i) - S22b.get(i);

    // Create pre-inverted constraint matrix C
    const C = Matrix([
        0.0, 0.5, 0.5,   0.0,   0.0,   0.0,
        0.5, 0.0, 0.5,   0.0,   0.0,   0.0,
        0.5, 0.5, 0.0,   0.0,   0.0,   0.0,
        0.0, 0.0, 0.0, -0.25,   0.0,   0.0,
        0.0, 0.0, 0.0,   0.0, -0.25,   0.0,
        0.0, 0.0, 0.0,   0.0,   0.0, -0.25,
    ], 6, 6);

    auto E = Multiply_Matrices(C, SS);

    auto SSS = Matrix(6, 6);

    Hessenberg_Form_Elementary(E.m.ptr, SSS.m.ptr, 6);

    double[6] eigen_real;
    double[6] eigen_imag;

    QR_Hessenberg_Matrix(E.m.ptr, SSS.m.ptr, eigen_real, eigen_imag, 6, 100);

    int index = 0;

    double maxval = eigen_real[0];

    for (int i = 1; i < 6; i++) {
        if (eigen_real[i] > maxval) {
            maxval = eigen_real[i];
            index = i;
        }
    }

    double[6] v1 = [
        SSS.get(index),      SSS.get(index +  6), SSS.get(index + 12),
        SSS.get(index + 18), SSS.get(index + 24), SSS.get(index + 30)
    ];

    // normalize v1
    {
        double norm = sqrt(
            v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] +
            v1[3] * v1[3] + v1[4] * v1[4] + v1[5] * v1[5]
        );

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
    }

    // Calculate v2 = S22a * v1      ( 4x1 = 4x6 * 6x1)
    double[4] v2;

    Multiply_Matrices(v2.ptr, S22a.m.ptr, 4, 6, v1.ptr, 1);

    double[10] v = [
         v1[0],  v1[1],  v1[2],  v1[3],  v1[4],
         v1[5], -v2[0], -v2[1], -v2[2], -v2[3]
    ];


    double[3 * 3] Q = [
        v[0], v[5], v[4],
        v[5], v[1], v[3],
        v[4], v[3], v[2]
    ];

    double[3] U = [ v[6], v[7], v[8] ];

    double[3 * 3] Q_1;

    for (int i = 0; i < 9; i++)
        Q_1[i] = Q[i];

    Choleski_LU_Decomposition(Q_1.ptr, 3);
    Choleski_LU_Inverse(Q_1.ptr, 3);

    // Calculate B = Q-1 * U   ( 3x1 = 3x3 * 3x1)
    double[3] B;

    Multiply_Matrices(B.ptr, Q_1.ptr, 3, 3, U.ptr, 1);

    Result result;
    result.combined_bias = [
        -B[0], // x-axis combined bias
        -B[1], // y-axis combined bias
        -B[2]  // z-axis combined bias
    ];

    // First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1)

    double[3] QB;

    Multiply_Matrices(QB.ptr, Q.ptr, 3, 3, B.ptr, 1);

    // Then calculate btqb = BT * QB    ( 1x1 = 1x3 * 3x1)
    double btqb;
    Multiply_Matrices(&btqb, B.ptr, 1, 3, QB.ptr, 1);


    // Calculate hmb = sqrt(btqb - J). (double J = v[9])
    double hmb = sqrt(btqb - v[9]);

    // Calculate SQ, the square root of matrix Q
    double[3 * 3] SSSS;

    Hessenberg_Form_Elementary(Q.ptr, SSSS.ptr, 3);

    double[3] eigen_real3, eigen_imag3;

    QR_Hessenberg_Matrix(Q.ptr, SSSS.ptr, eigen_real3, eigen_imag3, 3, 100);

    // normalize eigenvectors
    {
        double norm1 = sqrt(
            SSSS[0] * SSSS[0] +
            SSSS[3] * SSSS[3] +
            SSSS[6] * SSSS[6]
        );

        SSSS[0] /= norm1;
        SSSS[3] /= norm1;
        SSSS[6] /= norm1;

        double norm2 = sqrt(
            SSSS[1] * SSSS[1] +
            SSSS[4] * SSSS[4] +
            SSSS[7] * SSSS[7]
        );

        SSSS[1] /= norm2;
        SSSS[4] /= norm2;
        SSSS[7] /= norm2;

        double norm3 = sqrt(
            SSSS[2] * SSSS[2] +
            SSSS[5] * SSSS[5] +
            SSSS[8] * SSSS[8]
        );

        SSSS[2] /= norm3;
        SSSS[5] /= norm3;
        SSSS[8] /= norm3;
    }

    double[3 * 3] Dz;

    for (int i = 0; i < 9; i++)
        Dz[i] = 0.0;

    Dz[0] = sqrt(eigen_real3[0]);
    Dz[4] = sqrt(eigen_real3[1]);
    Dz[8] = sqrt(eigen_real3[2]);
    double[3 * 3] vdz;

    Multiply_Matrices(vdz.ptr, SSSS.ptr, 3, 3, Dz.ptr, 3);
    Transpose_Square_Matrix(SSSS.ptr, 3);

    double[3 * 3] SQ;
    Multiply_Matrices(SQ.ptr, vdz.ptr, 3, 3, SSSS.ptr, 3);

    const double hm = user_norm;

    double[3 * 3] A_1;

    for (int i = 0; i < 9; i++)
        A_1[i] = SQ[i] * hm / hmb;

    result.correction = [
        A_1[0 * 3], A_1[0 * 3 + 1], A_1[0 * 3 + 2],
        A_1[1 * 3], A_1[1 * 3 + 1], A_1[1 * 3 + 2],
        A_1[2 * 3], A_1[2 * 3 + 1], A_1[2 * 3 + 2]
    ];

    return result;
}

void main()
{
    // import std.stdio: writeln;
    // writeln(calculate_the_thing("mag.txt", 0.569));

    return;
}

unittest
{
    import std.stdio: writeln;
    import std.math.operations: isClose;

    Result expected;
    expected.combined_bias = [-0.021659, 0.013250, -0.026167];
    expected.correction = [0.956973, -0.017809, 0.006564, -0.017809, 0.964533, 0.003304, 0.006564, 0.003304, 1.036145];


    auto calc = calculate_the_thing("mag.txt", 0.569);

    writeln("calc", calc);
    writeln("expected", expected);

    for (auto i = 0; i < expected.combined_bias.length; i++)
        assert(isClose(calc.combined_bias[i], expected.combined_bias[i], 0.0001));

    for (auto i = 0; i < expected.correction.length; i++)
        assert(isClose(calc.correction[i], expected.correction[i], 0.0001));
}
