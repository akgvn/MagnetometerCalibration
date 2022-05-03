const std = @import("std");
const math = @import("math");
const Matrix = @import("types.zig").Matrix;
const Result = @import("types.zig").Result;
const mapi = @import("math_api.zig");

fn calculate_the_thing() Result {
    return .{ .bias = .{ 0.0, 0.0, 0.0 }, .corr = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
}

fn parse_line(line: []const u8) ![3]f64 {
    const parseFloat = std.fmt.parseFloat;
    const tokenize = std.mem.tokenize;

    var line_iterator = tokenize(u8, line, "\t\n\x0D"); // \x0D is ASCII 13, carriage return. parseFloat fails when the slice includes that.

    const first = try parseFloat(f64, line_iterator.next().?);
    const second = try parseFloat(f64, line_iterator.next().?);
    const third = try parseFloat(f64, line_iterator.next().?);

    return [3]f64{ first, second, third };
}

fn read_file_data_to_matrix() !Matrix {
    const cwd = std.fs.cwd();
    var file = try cwd.openFile("mag.txt", .{});
    defer file.close();

    const bufferedReader = std.io.bufferedReader;
    const reader = bufferedReader(file.reader()).reader();

    // var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    // defer arena.deinit();
    // const ally = &arena.allocator;

    var buf: [4096]u8 = undefined;
    while (try reader.readUntilDelimiterOrEof(&buf, '\n')) |line| {
        const line_data = try parse_line(line);
    }

    return Matrix.init(&.{}, 1, 1);
}

pub fn main() anyerror!void {
    _ = try read_file_data_to_matrix();

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
