"""Build isolated scanner with receiver shifts [-3,1], physical delays [-1,3]."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess


def build(output):
    receiver = Path("/home/zoxoai/chh/GPU/High/gnss-sdr-cuda-scan")
    source = receiver / "src/algorithms/tracking/libs/frozen_reference_open_loop_scan.cc"
    library = receiver / "build-cuda-scan/src/algorithms/tracking/libs"
    target = library / "CMakeFiles/frozen_reference_open_loop_scan.dir"
    original = source.read_text()
    changed = original
    replacements = {
        "-2.0 + 0.01 * static_cast<double>(tap)": "-3.0 + 0.01 * static_cast<double>(tap)",
        "-2.0F + 0.01F * static_cast<float>(tap)": "-3.0F + 0.01F * static_cast<float>(tap)",
    }
    counts = {}
    for old, new in replacements.items():
        counts[old] = changed.count(old)
        changed = changed.replace(old, new)
    if list(counts.values()) != [2, 1]:
        raise RuntimeError(f"unexpected scanner source: {counts}")
    output.mkdir(parents=True, exist_ok=True)
    copied = output / "gss_wide_scan.cc"
    copied.write_text(changed)
    obj = output / "gss_wide_scan.o"
    binary = output / "gss_wide_scan"
    flags = {}
    for line in (target / "flags.make").read_text().splitlines():
        if " = " in line:
            key, value = line.split(" = ", 1)
            flags[key] = shlex.split(value)
    compile_cmd = ["/usr/bin/c++", *flags["CXX_DEFINES"], *flags["CXX_INCLUDES"],
                   *flags["CXX_FLAGS"], "-c", str(copied), "-o", str(obj)]
    subprocess.run(compile_cmd, cwd=library, check=True)
    link = shlex.split((target / "link.txt").read_text())
    link[link.index("-o") + 1] = str(binary)
    link = [str(obj) if t.endswith("frozen_reference_open_loop_scan.cc.o") else t
            for t in link]
    subprocess.run(link, cwd=library, check=True)
    (output / "build_manifest.json").write_text(json.dumps({
        "source": str(source), "source_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "modified_sha256": hashlib.sha256(changed.encode()).hexdigest(),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "changes": counts, "compile": compile_cmd, "link": link,
        "receiver_shift_grid": [-3, 1, .01], "physical_delay_grid": [-1, 3, .01],
        "original_binary_modified": False,
    }, indent=2))


if __name__ == "__main__":
    build(Path("results/gss_replay_20261020/build").resolve())
