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

        pub fn isSquare() bool {
            return rows == cols;
        }

        pub fn substract(self: *const Self, rhs: *const Self) Self {
            var result = self.*;

            comptime var row: u32 = 0;
            comptime var col: u32 = 0;
            inline while (row < rows) : (row += 1) {
                inline while (col < cols) : (col += 1) {
                    result.getMut(row, col).* = self.get(row, col) - rhs.get(row, col);
                }
            }

            return result;
        }

        pub fn norm(self: *const Self) f64 {
            comptime assert(Self.cols == 1); // Must be Vector

            var sum: f64 = 0.0;

            comptime var col = 0;
            inline while (col < Self.cols) : (col += 1) {
                const val = self.get(0, col);
                sum += val * val;
            }

            return std.math.sqrt(sum); // TODO ??? is this right
        }

        pub fn normalized(self: *const Self) Self {
            comptime assert(Self.cols == 1); // Must be Vector
            var result = self.*;
            const n = self.norm();

            var col: u32 = 0;
            while (col < @TypeOf(result).cols) : (col += 1) {
                result.getMut(0, col).* = self.get(0, col) / n;
            }

            if (self.get(0, 0) < 0) {
                col = 0;
                while (col < @TypeOf(result).cols) : (col += 1) {
                    result.getMut(0, col).* = -1;
                }
            }

            return result;
        }

        // ---------- External function bindings ----------

        extern fn Get_Submatrix(S: [*]f64, mrows: c_int, mcols: c_int, A: [*]const f64, ncols: c_int, row: c_int, col: c_int) void;
        pub fn getSubmatrix(self: *const Self, comptime mrows: i32, comptime mcols: i32, row: i32, col: i32) Matrix(mrows, mcols) {
            var result = Matrix(mrows, mcols).zeroed();
            Get_Submatrix(&result.data, mrows, mcols, &self.data, cols, row, col);
            return result;
        }

        extern fn Choleski_LU_Decomposition(A: [*]f64, n: c_int) c_int;
        pub fn choleskiDecomposed(self: *const Self) !Self {
            comptime assert(Self.isSquare());
            var result = self.*;
            const cResult = Choleski_LU_Decomposition(&result.data, Self.cols);
            if (cResult != 0) {
                return error.CholeskiLUDecompositionFailed;
            }
            return result;
        }

        extern fn Choleski_LU_Inverse(LU: [*]f64, n: c_int) c_int;
        pub fn choleskiInversed(self: *const Self) !Self {
            comptime assert(Self.isSquare());
            var result = self.*;
            const cResult = Choleski_LU_Inverse(&result.data, Self.cols);
            if (cResult != 0) {
                return error.CholeskiLUInverseFailed;
            }
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
pub fn isMatrixValue(comptime T: type) bool {
    if (@typeInfo(T) != .Struct) return false;

    if (!@hasDecl(T, "rows")) return false;
    if (!@hasDecl(T, "cols")) return false;

    if (@TypeOf(T.rows) != u32) return false;
    if (@TypeOf(T.cols) != u32) return false;

    return Matrix(T.rows, T.cols) == T;
}

pub fn isMatrixRef(comptime T: type) bool {
    if (@typeInfo(T) == .Pointer)
        return isMatrixValue(@typeInfo(T).Pointer.child)
    else if (isMatrixValue(T))
        @compileError("Expected Matrix pointer, got Matrix!");
    return false;
}

pub fn isMatrix(comptime T: type) bool {
    return isMatrixValue(T) or isMatrixRef(T);
}

pub fn getDereferencedOrSame(comptime T: type) type {
    return switch (@typeInfo(T)) {
        .Pointer => |p| p.child,
        else => T,
    };
}

pub fn multipliedMatrixType(comptime A: type, comptime B: type) type {
    const dereferencedA = getDereferencedOrSame(A);
    const a1 = switch (@typeInfo(dereferencedA)) {
        .ErrorUnion => |eu| eu.payload,
        else => dereferencedA,
    };

    const dereferencedB = getDereferencedOrSame(B);
    const b1 = switch (@typeInfo(dereferencedB)) {
        .ErrorUnion => |eu| eu.payload,
        else => dereferencedB,
    };

    comptime assert(@typeInfo(a1) != .ErrorUnion);
    comptime assert(@typeInfo(b1) != .ErrorUnion);

    comptime assert(isMatrix(a1));
    comptime assert(isMatrix(b1));
    comptime assert(@field(a1, "cols") == @field(b1, "rows"));

    return Matrix(a1.rows, b1.cols);
}
