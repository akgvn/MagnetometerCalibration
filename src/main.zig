const std = @import("std");
const math = @import("math");
const Matrix = @import("types.zig").Matrix;
const Result = @import("types.zig").Result;
const mapi = @import("math_api.zig");

fn calculate_the_thing() Result {
    return .{ .bias = .{ 0.0, 0.0, 0.0 }, .corr = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
}

pub fn main() anyerror!void {
    var mat = [_]f64{
        0.956973,  -0.017809, 0.006564,
        -0.017809, 0.964533,  0.003304,
        0.006564,  0.003304,  1.036145,
    };

    var res: [9]f64 = .{0} ** 9;

    mapi.Matrix_x_Its_Transpose(&res, &mat, 3, 3);

    std.log.info("All your codebase are belong to us. {any}", .{res});
}

test "Mag.txt results test" {
    const result = calculate_the_thing();

    const expected = Result{ .bias = [_]f64{ -0.021659, 0.013250, -0.026167 }, .corr = [_]f64{
        0.956973,  -0.017809, 0.006564,
        -0.017809, 0.964533,  0.003304,
        0.006564,  0.003304,  1.036145,
    } };

    try std.testing.expectEqual(expected, result);
}
