# Modkit Installation for macOS

**Prerequisites:** Apple Silicon Mac, macOS 12.3+ (Monterey or later)

> **Why 12.3?** PyTorch MPS (Metal GPU backend) requires macOS 12.3+. The script exits early on older systems.

## Quick Start

```bash
bash mac_compile_modkit.sh ~/tools
```

Installs the latest modkit to `~/tools` using system Python. Takes 10–15 minutes.

---

## What Gets Installed

- Xcode Command Line Tools, Homebrew, Rust & Cargo
- Python virtual environment with PyTorch (GPU-enabled via Metal Performance Shaders)
- Compiled modkit binary at `~/tools/modkit/target/release/modkit`

---

## After Installation

The installer automatically adds the environment setup to `~/.zprofile`, so modkit is available in every new terminal session without any manual steps.

To activate in the **current** session immediately after install:

```bash
source ~/tools/setup_modkit_env.sh ~/tools
modkit --version
```

`setup_modkit_env.sh` configures:
- `LIBTORCH`, `DYLD_LIBRARY_PATH`, `LD_LIBRARY_PATH` — required for the modkit binary to find libtorch
- `PATH` — so you can type `modkit` directly
- `RAYON_NUM_THREADS` — automatically set to the number of **Performance cores** on your Mac

### RAYON_NUM_THREADS

Rayon is modkit's parallel processing library. On Apple Silicon, the script detects P-cores via:

```bash
sysctl -n hw.perflevel0.logicalcpu   # P-cores (used by default)
sysctl -n hw.perflevel1.logicalcpu   # E-cores (for reference)
```

To override for a single run:

```bash
RAYON_NUM_THREADS=8 modkit pileup input.bam output.bed
```

---

## Python Version Control

Three environment variables control Python selection:

| Variable | Values | Default | Purpose |
|----------|--------|---------|---------|
| `MODKIT_PYTHON_PROVIDER` | `auto`, `system`, `pyenv`, `uv` | `auto` | Which Python manager to use |
| `MODKIT_PYTHON_VERSION` | e.g. `3.11.9`, `3.12` | _(empty)_ | Request a specific Python version |
| `MODKIT_USE_UV` | `auto`, `0`, `1` | `auto` | Use `uv` for venv and pip operations |

### Provider Notes

- **`auto`**: Uses `uv` or `pyenv` if `MODKIT_PYTHON_VERSION` is set and they are available; otherwise falls back to system `python3`. If a version is requested but neither tool is available, the script warns and prompts before falling back.
- **`system`**: Uses `python3` from `$PATH`. Version cannot be controlled.
- **`pyenv`**: Installs the requested version if needed. Does not change your global pyenv version.
- **`uv`**: Fastest option. Installs Python and manages the venv/pip. Auto-installed via Homebrew if needed.

**Venv reuse:** If a virtual environment already exists from a different Python configuration, the script detects the mismatch and prompts you to recreate it.

### Examples

```bash
# Default: system Python
bash mac_compile_modkit.sh ~/tools

# Specific version via pyenv
MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash mac_compile_modkit.sh ~/tools

# Specific version via uv (fastest)
MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash mac_compile_modkit.sh ~/tools

# Force standard pip (no uv)
MODKIT_USE_UV=0 bash mac_compile_modkit.sh ~/tools

# Specific modkit version
bash mac_compile_modkit.sh ~/tools v0.5.0
```

---

## Common Issues & Fixes

### "Could not resolve a usable Python executable"

```bash
python3 --version           # check if Python is available
brew install python@3.11    # install if missing
```

### "PyTorch verification failed"

```bash
source ~/tools/setup_modkit_env.sh ~/tools
~/tools/venv_modkit/bin/python -m pip install --upgrade pip torch numpy
```

### "Compilation failed"

> Always source the environment before rebuilding. Without it, `LIBTORCH` is unset and the build will fail again.

```bash
source ~/tools/setup_modkit_env.sh ~/tools
cd ~/tools/modkit
cargo clean
cargo build --release --features accelerate,tch
```

### "MPS available: false"

- Confirm Apple Silicon: `uname -m` should print `arm64`
- Confirm macOS 12.3+: `sw_vers -productVersion`
- CPU-only fallback is used automatically if MPS is unavailable

---

## Help & Details

```bash
bash mac_compile_modkit.sh --help   # all options
cat ~/tools/installation_info.txt   # paths, versions, reproducible commands
modkit --help                       # modkit usage
```

