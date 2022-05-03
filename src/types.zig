pub fn Matrix(comptime rows: u32, comptime cols: u32) type {
    return struct {
        data: [rows * cols]f64,

        const Self = @This();

        pub fn init(data: []f64) Self {
            return .{ .data = data };
        }

        pub fn get(self: *Self, row: u64, col: u64) *f64 {
            return &self.data[row * self.rows + col];
        }
    };
}

pub const CalibrationResult = struct { bias: [3]f64, corr: [9]f64 };

pub const FileReadResult = struct { list: []f64, line_count: i32 };