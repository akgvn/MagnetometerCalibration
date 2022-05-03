pub const Matrix = struct {
    data: []f64,
    cols: u64,
    rows: u64,

    pub fn init(data: []f64, cols: u64, rows: u64) Matrix {
        return .{ .data = data, .cols = cols, .rows = rows };
    }
};

pub const Result = struct { bias: [3]f64, corr: [9]f64 };
