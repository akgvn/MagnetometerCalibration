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
const isMatrixRef = types.isMatrixRef;
const areMultipliableMatrices = types.areMultipliableMatrices;
const multipliedMatrixType = types.multipliedMatrixType;

extern fn Multiply_Matrices(C: [*]f64, A: [*]const f64, nrows: c_int, ncols: c_int, B: [*]const f64, mcols: c_int) void;
pub fn multiplyMatrices(A: anytype, B: anytype) multipliedMatrixType(@TypeOf(A), @TypeOf(B)) {
    var result = multipliedMatrixType(@TypeOf(A), @TypeOf(B)).zeroed();
    const typeA = types.getDereferencedOrSame(@TypeOf(A));
    const typeB = types.getDereferencedOrSame(@TypeOf(B));
    Multiply_Matrices(&result.data, &A.data, typeA.rows, typeA.cols, &B.data, typeB.cols);
    return result;
}

fn hessenbergReturnType(comptime T: type) type {
    assert(isMatrixRef(T));
    const n = @typeInfo(T).Pointer.child.rows;
    return Matrix(n, n);
}

extern fn Hessenberg_Form_Elementary(A: [*]const f64, S: [*]f64, n: c_int) c_int; // a is square, n is it's col count, s is output
pub fn hessenbergFormElementary(matrix: anytype) hessenbergReturnType(@TypeOf(matrix)) {
    const returnType = hessenbergReturnType(@TypeOf(matrix));
    var result: returnType = undefined;
    _ = Hessenberg_Form_Elementary(&matrix.data, &result.data, returnType.rows);
    return result;
}

extern fn QR_Hessenberg_Matrix(H: [*]const f64, S: [*]const f64, eigen_real: [*]f64, eigen_imag: [*]f64, n: c_int, max_iteration_count: c_int) c_int;
pub fn hessenbergQRMatrix(H: anytype, S: anytype, comptime n: i32, eigen_real: *Matrix(n, 1), eigen_imag: *Matrix(n, 1), comptime max_iteration_count: i32) i32 {
    comptime assert(isMatrixRef(@TypeOf(H)) and isMatrixRef(@TypeOf(S)));
    return QR_Hessenberg_Matrix(&H.data, &S.data, &eigen_real.data, &eigen_imag.data, n, max_iteration_count);
}