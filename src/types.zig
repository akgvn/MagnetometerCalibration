pub fn Matrix(comptime rows: u32, comptime cols: u32) type {
    return struct {
        data: [rows * cols]f64,

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

        pub fn getSubmatrix(self: *const Self, comptime mrows: i32, comptime mcols: i32, row: i32, col: i32) Matrix(mrows, mcols) {
            var result = Matrix(mrows, mcols).zeroed();

            const Get_Submatrix = @import("math_api.zig").Get_Submatrix;
            Get_Submatrix(&result.data, mrows, mcols, &self.data, cols, row, col);

            return result;
        }
    };
}

pub const CalibrationResult = struct { bias: [3]f64, corr: [9]f64 };

pub const FileReadResult = struct { list: []f64, line_count: i32 };
