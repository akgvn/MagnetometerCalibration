const std = @import("std");
const assert = std.debug.assert;
const types = @import("types.zig");
const Matrix = types.Matrix;

// Doesn't work for some reason!
// const c = @cImport({
//     @cInclude("math_functions.h");
// });

extern fn Matrix_x_Its_Transpose(result: [*]f64, matrix: [*]const f64, nrows: c_int, ncols: c_int) void;
pub fn multiplyMatrixWithTranspose(matrix: []f64, comptime nrows: i32, ncols: i32) Matrix(nrows, nrows) {
    var result: [nrows * nrows]f64 = undefined;
    Matrix_x_Its_Transpose(&result, matrix.ptr, nrows, ncols);
    return Matrix(nrows, nrows){ .data = result };
}


const isMatrix = types.isMatrix;
const areMultipliableMatrices = types.areMultipliableMatrices;
const multipliedMatrixType = types.multipliedMatrixType;

extern fn Multiply_Matrices(C: [*]f64, A: [*]const f64, nrows: c_int, ncols: c_int, B: [*]const f64, mcols: c_int) void;
pub fn multiplyMatrices(A: anytype, B: anytype) multipliedMatrixType(@TypeOf(A), @TypeOf(B)) {
    var result = multipliedMatrixType(@TypeOf(A), @TypeOf(B)).zeroed();
    const typeA = comptime if (std.meta.trait.isPtrTo(.Struct)(@TypeOf(A))) std.meta.Child(@TypeOf(A)) else @TypeOf(A);
    const typeB = comptime if (std.meta.trait.isPtrTo(.Struct)(@TypeOf(B))) std.meta.Child(@TypeOf(B)) else @TypeOf(B);
    Multiply_Matrices(&result.data, &A.data, typeA.rows, typeA.cols, &B.data, typeB.cols);
    return result;
}

extern fn Hessenberg_Form_Elementary(A: [*]f64, S: [*]f64, n: c_int) c_int;
// pub fn hessenbergFormElementary(A: Matrix, S: Matrix, n: i32) i32;
extern fn QR_Hessenberg_Matrix(H: [*]f64, S: [*]f64, eigen_real: [*]f64, eigen_imag: [*]f64, n: c_int, max_iteration_count: c_int) c_int;
// pub fn QR_Hessenberg_Matrix(H: Matrix, S: Matrix, eigen_real: Matrix, eigen_imag: Matrix, n: i32, max_iteration_count: i32) i32;
extern fn Transpose_Square_Matrix(A: [*]f64, n: c_int) void;
// pub fn Transpose_Square_Matrix(A: Matrix, n: i32) void;
