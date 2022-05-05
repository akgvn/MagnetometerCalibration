const std = @import("std");

pub fn Matrix(comptime rowCount: u32, comptime colCount: u32) type {
    return struct {
        data: [rowCount * colCount]f64,

        pub const rows = rowCount;
        pub const cols = colCount;

        const Self = @This();

        pub fn zeroed() Self {
            return .{ .data = .{0} ** (rows * cols) };
        }

        pub fn init(data: []f64) Self {
            return .{ .data = data };
        }

        pub fn get(self: *Self, row: u64, col: u64) *f64 {
            return &self.data[row * self.rows + col];
        }

        extern fn Get_Submatrix(S: [*]f64, mrows: c_int, mcols: c_int, A: [*]const f64, ncols: c_int, row: c_int, col: c_int) void;
        pub fn getSubmatrix(self: *const Self, comptime mrows: i32, comptime mcols: i32, row: i32, col: i32) Matrix(mrows, mcols) {
            var result = Matrix(mrows, mcols).zeroed();

            // const Get_Submatrix = @import("math_api.zig").Get_Submatrix;
            Get_Submatrix(&result.data, mrows, mcols, &self.data, cols, row, col);

            return result;
        }

        pub fn isSquare() bool {
            return rows == cols;
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
    const assert = std.debug.assert;
    const dereferencedAType = if (std.meta.trait.isPtrTo(.Struct)(A)) std.meta.Child(A) else A;
    const dereferencedBType = if (std.meta.trait.isPtrTo(.Struct)(B)) std.meta.Child(B) else B;
    comptime assert(areMultipliableMatrices(dereferencedAType, dereferencedBType));

    return Matrix(dereferencedAType.cols, dereferencedBType.cols);
}
