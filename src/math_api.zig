// Doesn't work for some reason!
// const c = @cImport({
//     @cInclude("math_functions.h");
// });

pub extern fn Matrix_x_Its_Transpose(result: [*c]f64, matrix: [*c]f64, nrows: c_int, ncols: c_int) void;
pub extern fn Choleski_LU_Decomposition(A: *f64, n: c_int) c_int;
pub extern fn Choleski_LU_Inverse(LU: *f64, n: c_int) c_int;
pub extern fn Get_Submatrix(S: *f64, mrows: c_int, mcols: c_int, A: *f64, ncols: c_int, row: c_int, col: c_int) void;
pub extern fn Hessenberg_Form_Elementary(A: *f64, S: *f64, n: c_int) c_int;
pub extern fn Multiply_Matrices(C: *f64, A: *f64, nrows: c_int, ncols: c_int, B: *f64, mcols: c_int) void;
pub extern fn QR_Hessenberg_Matrix(H: *f64, S: *f64, eigen_real: *f64, eigen_imag: *f64, n: c_int, max_iteration_count: c_int) c_int;
pub extern fn Transpose_Square_Matrix(A: *f64, n: c_int) void;
