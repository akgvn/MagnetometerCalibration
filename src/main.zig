const std = @import("std");
const math = @import("math");
const types = @import("types.zig");
const mapi = @import("math_api.zig");
const Allocator = std.mem.Allocator;

const Matrix = types.Matrix;
const CalibrationResult = types.CalibrationResult;
const FileReadResult = types.FileReadResult;

fn calculateTheThing() !CalibrationResult {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    var fileData = try readFileData(arena.allocator());

    var S = mapi.multiplyMatrixWithTranspose(fileData.list, 10, fileData.line_count);
    arena.deinit();

    const S11 = S.getSubmatrix(6, 6, 0, 0);
    const S12 = S.getSubmatrix(6, 4, 0, 6);
    const S12t = S.getSubmatrix(4, 6, 6, 0);
    const S22 = S.getSubmatrix(4, 4, 6, 6);

    const S22_1 = S22.choleskiDecomposed().choleskiInversed();

    const S22a = mapi.multiplyMatrices(&S22_1, &S12t);
    const S22b = mapi.multiplyMatrices(&S12, &S22a);

    const SS = S11.substract(&S22b);

    // Create pre-inverted constraint matrix C
    const C = Matrix(6, 6).init(.{
        0.0, 0.5, 0.5, 0.0,   0.0,   0.0,
        0.5, 0.0, 0.5, 0.0,   0.0,   0.0,
        0.5, 0.5, 0.0, 0.0,   0.0,   0.0,
        0.0, 0.0, 0.0, -0.25, 0.0,   0.0,
        0.0, 0.0, 0.0, 0.0,   -0.25, 0.0,
        0.0, 0.0, 0.0, 0.0,   0.0,   -0.25,
    });

    const E = mapi.multiplyMatrices(&C, &SS);

    const SSS = mapi.hessenbergFormElementary(&E);

    var eigen_real: Matrix(1, 6) = undefined;
    var eigen_imag: Matrix(1, 6) = undefined;
    _ = mapi.hessenbergQRMatrix(&E, &SSS, 6, &eigen_real, &eigen_imag, 100);

    std.log.info("{}", .{SSS});

    return CalibrationResult{ .bias = .{ 0.0, 0.0, 0.0 }, .corr = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
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

fn get(item: []f64, row: i32, col: i32, rows: i32) *f64 {
    const pos = @intCast(usize, row * rows + col);
    return &item[pos];
}

fn transpose(mat: []f64, rows: i32, cols: i32) []f64 {
    var row: i32 = 0;
    var col: i32 = 0;
    while (row < rows) : (row += 1) {
        while (col < cols) : (col += 1) {
            get(mat, row, col, rows).* = get(mat, col, row, rows).*;
        }
    }
    return mat;
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

    return FileReadResult{ .list = transpose(list.toOwnedSlice(), 10, line_count), .line_count = line_count };
}

pub fn main() anyerror!void {
    _ = try calculateTheThing();
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
