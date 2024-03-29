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
    const S = (D * D.transpose());

    const S11  = S.submatrix(6, 6, 0, 0);

    const S12  = S.submatrix(6, 4, 0, 6);
    const S12t = S.submatrix(4, 6, 6, 0);
    const S22  = S.submatrix(4, 4, 6, 6);

    auto S22_1 = Matrix(4, 4); // double[4 * 4] S22_1;

    for (int i = 0; i < 16; i++)
        S22_1.get(i) = S22.get(i);

    S22_1.choleski_lu_decomposition();
    Choleski_LU_Inverse(S22_1.m.ptr, 4);

    // Calculate S22a = S22_1 * S12t   4*6 = 4x4 * 4x6   C = AB

    const S22a = (S22_1 * S12t);

    // Then calculate S22b = S12 * S22a      ( 6x6 = 6x4 * 4x6)

    const S22b = (S12 * S22a);

    // Calculate SS = S11 - S22b

    auto SS = Matrix(6, 6);

    for (int i = 0; i < 36; i++)
        SS.get(i) = S11.get(i) - S22b.get(i);

    // Create pre-inverted constraint matrix C
    static immutable C = Matrix([
        0.0, 0.5, 0.5,   0.0,   0.0,   0.0,
        0.5, 0.0, 0.5,   0.0,   0.0,   0.0,
        0.5, 0.5, 0.0,   0.0,   0.0,   0.0,
        0.0, 0.0, 0.0, -0.25,   0.0,   0.0,
        0.0, 0.0, 0.0,   0.0, -0.25,   0.0,
        0.0, 0.0, 0.0,   0.0,   0.0, -0.25,
    ], 6, 6);

    auto E = (C * SS);

    auto SSS = Matrix(6, 6);

    Hessenberg_Form_Elementary(E.m.ptr, SSS.m.ptr, 6);

    auto eigen_real = Matrix(6, 1);
    auto eigen_imag = Matrix(6, 1);

    QR_Hessenberg_Matrix(E.m.ptr, SSS.m.ptr, eigen_real.m, eigen_imag.m, 6, 100);

    int index = 0;
    double maxval = eigen_real.get(0);
    foreach (i; 1..6) {
        if (eigen_real.get(i) > maxval) {
            maxval = eigen_real.get(i);
            index = i;
        }
    }

    auto v1 = Matrix([
        SSS.get(index),      SSS.get(index +  6), SSS.get(index + 12),
        SSS.get(index + 18), SSS.get(index + 24), SSS.get(index + 30)
    ], 6, 1).normalized();

    // Calculate v2 = S22a * v1      ( 4x1 = 4x6 * 6x1)
    const v2 = (S22a * v1);

    auto v = Matrix([
         v1.get(0),  v1.get(1),  v1.get(2),  v1.get(3),  v1.get(4),
         v1.get(5), -v2.get(0), -v2.get(1), -v2.get(2), -v2.get(3)
    ], 10, 1);


    auto Q = Matrix([
        v.get(0), v.get(5), v.get(4),
        v.get(5), v.get(1), v.get(3),
        v.get(4), v.get(3), v.get(2)
    ], 3, 3);

    auto U = Matrix([ v.get(6), v.get(7), v.get(8) ], 3, 1);

    auto Q_1 = Matrix(Q);

    Q_1.choleski_lu_decomposition();
    Choleski_LU_Inverse(Q_1.m.ptr, 3);

    // Calculate B = Q-1 * U ( 3x1 = 3x3 * 3x1)
    const B = (Q_1 * U);

    Result result;
    result.combined_bias = [
        -B.get(0), // x-axis combined bias
        -B.get(1), // y-axis combined bias
        -B.get(2)  // z-axis combined bias
    ];

    // First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1)

    const QB = (Q * B);
    const BT = B.transpose();

    // Then calculate btqb = BT * QB    ( 1x1 = 1x3 * 3x1)
    const double btqb = (BT * QB).flattened_value();

    // Calculate hmb = sqrt(btqb - J). (double J = v[9])
    const double hmb = sqrt(btqb - v.get(9));

    // Calculate SQ, the square root of matrix Q
    auto SSSS = Matrix(3, 3);

    Hessenberg_Form_Elementary(Q.m.ptr, SSSS.m.ptr, 3);

    auto eigen_real3 = Matrix(3, 1), eigen_imag3 = Matrix(3, 1);

    QR_Hessenberg_Matrix(Q.m.ptr, SSSS.m.ptr, eigen_real3.m, eigen_imag3.m, 3, 100);

    // normalize eigenvectors
    {
        double norm1 = sqrt(
            SSSS.get(0) * SSSS.get(0) +
            SSSS.get(3) * SSSS.get(3) +
            SSSS.get(6) * SSSS.get(6)
        );

        SSSS.get(0) /= norm1;
        SSSS.get(3) /= norm1;
        SSSS.get(6) /= norm1;

        double norm2 = sqrt(
            SSSS.get(1) * SSSS.get(1) +
            SSSS.get(4) * SSSS.get(4) +
            SSSS.get(7) * SSSS.get(7)
        );

        SSSS.get(1) /= norm2;
        SSSS.get(4) /= norm2;
        SSSS.get(7) /= norm2;

        double norm3 = sqrt(
            SSSS.get(2) * SSSS.get(2) +
            SSSS.get(5) * SSSS.get(5) +
            SSSS.get(8) * SSSS.get(8)
        );

        SSSS.get(2) /= norm3;
        SSSS.get(5) /= norm3;
        SSSS.get(8) /= norm3;
    }

    const _00 = sqrt(eigen_real3.get(0));
    const _11 = sqrt(eigen_real3.get(1));
    const _22 = sqrt(eigen_real3.get(2));
    const Dz = Matrix([
        _00, 0.0, 0.0,
        0.0, _11, 0.0,
        0.0, 0.0, _22
    ], 3, 3);

    const vdz = (SSSS * Dz);

    const SQ = (vdz * SSSS.transpose());

    const double hm = user_norm;

    auto A_1 = Matrix(3, 3);

    for (int i = 0; i < 9; i++)
        A_1.get(i) = SQ.get(i) * hm / hmb;

    result.correction = A_1.m;
    return result;
}

void main()
{
    import std.stdio: writeln;
    writeln(calculate_the_thing("mag.txt", 0.569));

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
