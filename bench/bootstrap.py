"""Locate the project, manage .venv, and incrementally build the Fortran C ABI."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ENV_READY = "GJKEPA_LAUNCHER_READY"
LIB_NAMES = (
    "gjkepa_c.dll",
    "libgjkepa_c.dll",
    "libgjkepa_c.so",
    "libgjkepa_c.dylib",
    "gjkepa_c.so",
)


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def venv_python(root: Path) -> Path:
    if os.name == "nt":
        return root / ".venv" / "Scripts" / "python.exe"
    return root / ".venv" / "bin" / "python"


def _run(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> int:
    print(">", " ".join(cmd), flush=True)
    completed = subprocess.run(cmd, cwd=cwd, env=env, check=False)
    return int(completed.returncode)


def ensure_venv(root: Path) -> Path:
    py = venv_python(root)
    if py.is_file():
        return py
    print(f"Creating virtualenv at {root / '.venv'} (requires network only if pip later installs).", flush=True)
    creator = Path(sys.executable)
    code = _run([str(creator), "-m", "venv", str(root / ".venv")])
    if code != 0 or not py.is_file():
        raise SystemExit(
            "Failed to create .venv. Install Python 3.10+ and ensure `python -m venv` works."
        )
    return py


def relaunch_in_venv(root: Path) -> None:
    py = ensure_venv(root)
    current = Path(sys.executable).resolve()
    target = py.resolve()
    if current == target:
        os.environ[ENV_READY] = "1"
        return
    if os.environ.get(ENV_READY) == "1":
        return
    os.environ[ENV_READY] = "1"
    argv = [str(target), str((root / "launcher.py").resolve()), *sys.argv[1:]]
    print(f"Re-launching inside {target}", flush=True)
    os.execv(str(target), argv)


def ensure_python_deps(root: Path, *, headless: bool) -> None:
    req = root / "bench" / "requirements.txt"
    try:
        import numpy  # noqa: F401
        if not headless:
            import pyvista  # noqa: F401
        return
    except ImportError:
        pass
    print("Python packages missing; installing from bench/requirements.txt (needs network once).", flush=True)
    code = _run([sys.executable, "-m", "pip", "install", "-r", str(req)])
    if code != 0:
        raise SystemExit("pip install failed. Check network or install numpy and pyvista into .venv.")
    try:
        import numpy  # noqa: F401
        if not headless:
            import pyvista  # noqa: F401
    except ImportError as exc:
        raise SystemExit(f"Dependencies still missing after pip install: {exc}") from exc


def _which(name: str) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    extras = []
    local = os.environ.get("LOCALAPPDATA")
    if local:
        extras.append(Path(local) / "Microsoft" / "WinGet" / "Packages")
    extras.extend(
        [
            Path(r"C:\Program Files\mingw64\bin"),
            Path(r"C:\msys64\mingw64\bin"),
            Path("/usr/bin"),
        ]
    )
    for folder in extras:
        if folder.is_file():
            continue
        direct = folder / name
        if folder.is_dir() and (direct.is_file() or (direct.with_suffix(".exe")).is_file()):
            cand = direct if direct.is_file() else direct.with_suffix(".exe")
            return str(cand)
        if folder.is_dir() and folder.name == "Packages":
            matches = list(folder.glob(f"**/bin/{name}.exe")) + list(folder.glob(f"**/bin/{name}"))
            if matches:
                return str(matches[0])
    return None


def find_tools() -> dict[str, str]:
    cmake = _which("cmake")
    if cmake is None:
        raise SystemExit(
            "CMake not found on PATH.\n"
            "Install CMake >= 3.20 (https://cmake.org/download/) and reopen the terminal."
        )
    fortran = _which("gfortran") or _which("ifx") or _which("ifort")
    if fortran is None:
        raise SystemExit(
            "No Fortran compiler on PATH.\n"
            "Install gfortran (WinLibs / MinGW / MSYS2) or Intel ifx/ifort, then retry."
        )
    ninja = _which("ninja")
    generator = "Ninja" if ninja else None
    return {
        "cmake": cmake,
        "fortran": fortran,
        "ninja": ninja or "",
        "generator": generator or "",
    }


def build_dir(root: Path) -> Path:
    return root / "build"


def configure_and_build(root: Path) -> Path:
    tools = find_tools()
    bdir = build_dir(root)
    cache = bdir / "CMakeCache.txt"
    env = os.environ.copy()
    fortran_dir = str(Path(tools["fortran"]).resolve().parent)
    env["PATH"] = fortran_dir + os.pathsep + env.get("PATH", "")

    if not cache.is_file():
        cmd = [tools["cmake"], "-S", str(root), "-B", str(bdir)]
        if tools["generator"]:
            cmd += ["-G", tools["generator"]]
        cmd += [
            f"-DCMAKE_Fortran_COMPILER={tools['fortran']}",
            "-DCMAKE_BUILD_TYPE=Release",
        ]
        if _run(cmd, env=env) != 0:
            raise SystemExit("CMake configure failed. Fix the compiler/CMake error above and retry.")

    build_cmd = [tools["cmake"], "--build", str(bdir), "--target", "gjkepa_c"]
    if (bdir / "build.ninja").is_file() or tools["generator"] == "Ninja":
        pass
    else:
        build_cmd += ["--config", "Release"]
    if _run(build_cmd, env=env) != 0:
        raise SystemExit(
            "CMake build failed. Refusing to load a previous shared library. "
            "See the compiler output above."
        )
    lib = find_shared_library(bdir)
    if lib is None:
        raise SystemExit(f"Build reported success but no {LIB_NAMES} was found under {bdir}.")
    return lib


def find_shared_library(bdir: Path) -> Path | None:
    hits: list[Path] = []
    for name in LIB_NAMES:
        hits.extend(p for p in bdir.rglob(name) if p.is_file())
    if not hits:
        return None
    hits.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return hits[0]


def prepare(root: Path, *, headless: bool = False) -> Path:
    relaunch_in_venv(root)
    ensure_python_deps(root, headless=headless)
    return configure_and_build(root)
