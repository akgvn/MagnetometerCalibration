const std = @import("std");
const math = @import("math");
const types = @import("types.zig");
const mapi = @import("math_api.zig");
const Allocator = std.mem.Allocator;

const Matrix = types.Matrix;
const CalibrationResult = types.CalibrationResult;
const FileReadResult = types.FileReadResult;

fn calculateTheThing() CalibrationResult {
    return .{ .bias = .{ 0.0, 0.0, 0.0 }, .corr = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
}

fn parseLine(line: []const u8) ![3]f64 {
    const parseFloat = std.fmt.parseFloat;
    const tokenize = std.mem.tokenize;

    var line_iterator = tokenize(u8, line, "\t\n\r");

    const first = try parseFloat(f64, line_iterator.next().?);
    const second = try parseFloat(f64, line_iterator.next().?);
    const third = try parseFloat(f64, line_iterator.next().?);

    return [3]f64{ first, second, third };
}

fn readFileData(allocator: Allocator) !FileReadResult {
    var file = try std.fs.cwd().openFile("mag.txt", .{});
    defer file.close();

    const bufferedReader = std.io.bufferedReader;
    const reader = bufferedReader(file.reader()).reader();

    var list = std.ArrayList(f64).init(allocator);

    var buf: [128]u8 = undefined;
    var line_count: i32 = 0;
    while (try reader.readUntilDelimiterOrEof(&buf, '\n')) |line| {
        const line_data = try parseLine(line);
        const x = line_data[0];
        const y = line_data[1];
        const z = line_data[2];

        try list.append(x * x);
        try list.append(y * y);
        try list.append(z * z);
        try list.append(2.0 * y * z);
        try list.append(2.0 * x * z);
        try list.append(2.0 * x * y);
        try list.append(2.0 * x);
        try list.append(2.0 * y);
        try list.append(2.0 * z);
        try list.append(1.0);

        line_count += 1;
    }

    return FileReadResult{ .list = list.toOwnedSlice(), .line_count = line_count };
}

pub fn main() anyerror!void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    var fileData = try readFileData(arena.allocator());

    const S = mapi.multiplyMatrixWithTranspose(fileData.list, 10, fileData.line_count);
    arena.deinit();

    std.log.info("In main, result of mapi.Matrix_x_Its_Transpose is {any}", .{S});
}

test "Mag.txt test" {
    const result = calculateTheThing();

    const expected = CalibrationResult{ .bias = [_]f64{ -0.021659, 0.013250, -0.026167 }, .corr = [_]f64{
        0.956973,  -0.017809, 0.006564,
        -0.017809, 0.964533,  0.003304,
        0.006564,  0.003304,  1.036145,
    } };

    try std.testing.expectEqual(expected, result);
}
