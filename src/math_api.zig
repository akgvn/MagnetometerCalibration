const Matrix = @import("types.zig").Matrix;

// Doesn't work for some reason!
// const c = @cImport({
//     @cInclude("math_functions.h");
// });

extern fn Matrix_x_Its_Transpose(result: [*]f64, matrix: [*]const f64, nrows: c_int, ncols: c_int) void;
extern fn Choleski_LU_Decomposition(A: [*]f64, n: c_int) c_int;
extern fn Choleski_LU_Inverse(LU: [*]f64, n: c_int) c_int;
pub extern fn Get_Submatrix(S: [*]f64, mrows: c_int, mcols: c_int, A: [*]const f64, ncols: c_int, row: c_int, col: c_int) void;
extern fn Hessenberg_Form_Elementary(A: [*]f64, S: [*]f64, n: c_int) c_int;
extern fn Multiply_Matrices(C: [*]f64, A: [*]f64, nrows: c_int, ncols: c_int, B: [*]f64, mcols: c_int) void;
extern fn QR_Hessenberg_Matrix(H: [*]f64, S: [*]f64, eigen_real: [*]f64, eigen_imag: [*]f64, n: c_int, max_iteration_count: c_int) c_int;
extern fn Transpose_Square_Matrix(A: [*]f64, n: c_int) void;

pub fn multiplyMatrixWithTranspose(matrix: []f64, comptime nrows: i32, ncols: i32) Matrix(nrows, nrows) {
    var result: [nrows * nrows]f64 = undefined;
    Matrix_x_Its_Transpose(&result, matrix.ptr, nrows, ncols);
    return Matrix(nrows, nrows){ .data = result };
}

pub fn choleskiDecomposition(comptime n: i32, A: Matrix(n, n)) Matrix(n, n) {
    var result = A;
    Choleski_LU_Decomposition(&result.data, n);
    return result;
}

pub fn choleskiInverse(comptime n: i32, LU: Matrix(n, n)) Matrix(n, n) {
    var result = LU;
    Choleski_LU_Inverse(&result.data, n);
    return result;
}


// pub fn hessenbergFormElementary(A: Matrix, S: Matrix, n: i32) i32;
// pub fn multiplyMatrices(C: Matrix, A: Matrix, nrows: i32, ncols: i32, B: Matrix, mcols: i32) void;
// pub fn QR_Hessenberg_Matrix(H: Matrix, S: Matrix, eigen_real: Matrix, eigen_imag: Matrix, n: i32, max_iteration_count: i32) i32;
// pub fn Transpose_Square_Matrix(A: Matrix, n: i32) void;
