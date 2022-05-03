const std = @import("std");

const Matrix = struct {
    data: []f64,
    cols: u64,
    rows: u64,

    fn init(cols: u64, rows: u64) Matrix {
        return .{ .data = null, .cols = cols, .rows = rows };
    }
};

const Result = struct { bias: [3]f64, corr: [9]f64 };

fn calculate_the_thing() Result {
    return .{ .bias = .{ 0.0, 0.0, 0.0 }, .corr = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
}

pub fn main() anyerror!void {
    std.log.info("All your codebase are belong to us.", .{});
}

test "basic test" {
    const result = calculate_the_thing();

    const combined_bias = [_]f64{ -0.021659, 0.013250, -0.026167 };

    const correction = [_]f64{
        0.956973,  -0.017809, 0.006564,
        -0.017809, 0.964533,  0.003304,
        0.006564,  0.003304,  1.036145,
    };

    try std.testing.expectEqual(result, .{ .bias = combined_bias, .corr = correction });
}
