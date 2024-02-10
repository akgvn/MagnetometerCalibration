const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    // const mathLib = b.addStaticLibrarySource("math", null);
    // mathLib.linkLibC();
    // mathLib.install();

    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    var exe = b.addExecutable(.{
        .name = "MagnetometerCalibration",
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });

    const CFlags = &.{};
    exe.addCSourceFile(.{
        .file = .{
            .path = "src/math_functions.c",
        },
        .flags = CFlags,
    });
    exe.linkLibC();

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    const unit_tests = b.addTest(.{
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    unit_tests.addCSourceFile(.{
        .file = .{
            .path = "src/math_functions.c",
        },
        .flags = CFlags,
    });
    unit_tests.linkLibC();

    const run_unit_tests = b.addRunArtifact(unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&exe.step);
    test_step.dependOn(&run_unit_tests.step);
}
