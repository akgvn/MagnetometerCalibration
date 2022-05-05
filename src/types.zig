const std = @import("std");
const assert = std.debug.assert;

pub fn Matrix(comptime rowCount: u32, comptime colCount: u32) type {
    return struct {
        data: [rowCount * colCount]f64,

        pub const rows = rowCount;
        pub const cols = colCount;

        const Self = @This();

        pub fn zeroed() Self {
            return .{ .data = .{0} ** (rows * cols) };
        }

        pub fn init(data: [rowCount * colCount]f64) Self {
            return .{ .data = data };
        }

        pub fn getMut(self: *Self, row: u64, col: u64) *f64 {
            return &self.data[row * Self.rows + col];
        }

        pub fn get(self: *const Self, row: u64, col: u64) f64 {
            return self.data[row * Self.rows + col];
        }

        extern fn Get_Submatrix(S: [*]f64, mrows: c_int, mcols: c_int, A: [*]const f64, ncols: c_int, row: c_int, col: c_int) void;
        pub fn getSubmatrix(self: *const Self, comptime mrows: i32, comptime mcols: i32, row: i32, col: i32) Matrix(mrows, mcols) {
            var result = Matrix(mrows, mcols).zeroed();
            Get_Submatrix(&result.data, mrows, mcols, &self.data, cols, row, col);
            return result;
        }

        pub fn isSquare() bool {
            return rows == cols;
        }

        pub fn substract(self: *const Self, rhs: *const Self) Self {
            var result = self.*;

            var row: u32 = 0;
            var col: u32 = 0;
            while (row < rows) : (row += 1) {
                while (col < cols) : (col += 1) {
                    result.getMut(row, col).* = self.get(row, col) - rhs.get(row, col);
                }
            }

            return result;
        }

        extern fn Choleski_LU_Decomposition(A: [*]f64, n: c_int) c_int;
        pub fn choleskiDecomposed(self: *const Self) Self {
            comptime assert(Self.isSquare());
            var result = self.*;
            _ = Choleski_LU_Decomposition(&result.data, Self.cols);
            return result;
        }

        extern fn Choleski_LU_Inverse(LU: [*]f64, n: c_int) c_int;
        pub fn choleskiInversed(self: *const Self) Self {
            comptime assert(Self.isSquare());
            var result = self.*;
            _ = Choleski_LU_Inverse(&result.data, Self.cols);
            return result;
        }

        extern fn Transpose_Square_Matrix(A: [*]f64, n: c_int) void;
        pub fn transposed(self: *Self) Self {
            comptime assert(Self.isSquare()); // Only square matrices are supported for now!
            var result = self.*;
            _ = Transpose_Square_Matrix(&result.data, Self.cols);
            return result;
        }
    };
}

pub const CalibrationResult = struct { bias: [3]f64, corr: [9]f64 };

pub const FileReadResult = struct { list: []f64, line_count: i32 };

/// Only returns true for types produced by calls to `Matrix`.
pub fn isMatrix(comptime T: type) bool {
    if (@typeInfo(T) == .Pointer) return isMatrix(@typeInfo(T).Pointer.child);
    if (@typeInfo(T) != .Struct) return false;

    if (!@hasDecl(T, "rows")) return false;
    if (!@hasDecl(T, "cols")) return false;

    if (!std.meta.declarationInfo(T, "rows").is_pub) return false;
    if (!std.meta.declarationInfo(T, "cols").is_pub) return false;

    if (@TypeOf(T.rows) != u32) return false;
    if (@TypeOf(T.cols) != u32) return false;

    return Matrix(T.rows, T.cols) == T;
}

pub fn areMultipliableMatrices(comptime T: type, comptime K: type) bool {
    return (isMatrix(T) and isMatrix(K) and (@field(T, "cols") == @field(K, "rows")));
}

pub fn multipliedMatrixType(comptime A: type, comptime B: type) type {
    const dereferencedAType = if (std.meta.trait.isPtrTo(.Struct)(A)) std.meta.Child(A) else A;
    const dereferencedBType = if (std.meta.trait.isPtrTo(.Struct)(B)) std.meta.Child(B) else B;
    comptime assert(areMultipliableMatrices(dereferencedAType, dereferencedBType));

    return Matrix(dereferencedAType.rows, dereferencedBType.cols);
}
