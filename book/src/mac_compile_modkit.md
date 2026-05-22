# Modkit Installation for macOS

**Prerequisites:** Apple Silicon Mac, macOS 11.0+

> **MPS GPU builds only:** macOS 12.3+ (Monterey or later) is required when `MODKIT_MPS_SUPPORT=1`. The script exits early on older systems for MPS builds only. Builds without MPS (`MODKIT_MPS_SUPPORT=0`, the default) run on macOS 11.0+.

## Quick Start

```bash
# Without MPS GPU support (default)
bash mac_compile_modkit.sh ~/tools

# With MPS GPU support (requires Python + PyTorch)
MODKIT_MPS_SUPPORT=1 bash mac_compile_modkit.sh ~/tools
```

The default (no-MPS) build installs modkit without Python or PyTorch. Takes 5–10 minutes. The MPS build adds Metal GPU acceleration and takes 10–15 minutes.

---

## What Gets Installed

**Always installed:**
- Xcode Command Line Tools, Homebrew, Rust & Cargo
- Compiled modkit binary at `~/tools/modkit/target/release/modkit`

**Additionally installed when `MODKIT_MPS_SUPPORT=1`:**
- Python virtual environment with PyTorch (GPU-enabled via Metal Performance Shaders)

---

## After Installation

The installer automatically adds the environment setup to `~/.zprofile`, so modkit is available in every new terminal session **silently** — no output is printed on shell startup.

To activate in the **current** session immediately after install:

```bash
source ~/tools/setup_modkit_env.sh ~/tools
modkit --version
```

To see full environment details:

```bash
source ~/tools/setup_modkit_env.sh ~/tools --verbose
```

`setup_modkit_env.sh` configures:
- `PATH` — so you can type `modkit` directly
- `LIBTORCH`, `DYLD_LIBRARY_PATH`, `LD_LIBRARY_PATH` — only set for MPS builds; required for the binary to find libtorch

---

## MPS GPU Support

`MODKIT_MPS_SUPPORT=1` enables Metal Performance Shaders GPU acceleration. This requires Python and PyTorch, which the script installs automatically.

| Variable | Values | Default | Purpose |
|----------|--------|---------|--------|
| `MODKIT_MPS_SUPPORT` | `0`, `1` | `0` | Enable MPS GPU support |

## Python Version Control

> These variables are only used when `MODKIT_MPS_SUPPORT=1`.

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
# No MPS (default): no Python or PyTorch needed
bash mac_compile_modkit.sh ~/tools

# MPS build with system Python
MODKIT_MPS_SUPPORT=1 bash mac_compile_modkit.sh ~/tools

# MPS build, specific Python version via pyenv
MODKIT_MPS_SUPPORT=1 MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash mac_compile_modkit.sh ~/tools

# MPS build, specific Python version via uv (fastest)
MODKIT_MPS_SUPPORT=1 MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash mac_compile_modkit.sh ~/tools

# MPS build, force standard pip (no uv)
MODKIT_MPS_SUPPORT=1 MODKIT_USE_UV=0 bash mac_compile_modkit.sh ~/tools

# Specific modkit version (no MPS)
bash mac_compile_modkit.sh ~/tools v0.5.0

# Specific modkit version (with MPS)
MODKIT_MPS_SUPPORT=1 bash mac_compile_modkit.sh ~/tools v0.5.0
```

---

## Common Issues & Fixes

### "Could not resolve a usable Python executable"

> Only relevant for `MODKIT_MPS_SUPPORT=1` builds.

```bash
python3 --version           # check if Python is available
brew install python@3.11    # install if missing
```

### "PyTorch verification failed"

> Only relevant for `MODKIT_MPS_SUPPORT=1` builds.

```bash
source ~/tools/setup_modkit_env.sh ~/tools
~/tools/venv_modkit/bin/python -m pip install --upgrade pip torch numpy
```

### "Compilation failed"

**No-MPS build:**

```bash
cd ~/tools/modkit
cargo clean
cargo build --release --features accelerate
```

**MPS build:**

> Always source the environment before rebuilding. Without it, `LIBTORCH` is unset and the build will fail again.

```bash
source ~/tools/setup_modkit_env.sh ~/tools
cd ~/tools/modkit
cargo clean
cargo build --release --features accelerate,tch
```

### "MPS available: false"

MPS support is opt-in (`MODKIT_MPS_SUPPORT=1`). If you built without MPS, this is expected — rebuild with `MODKIT_MPS_SUPPORT=1` to enable GPU acceleration.

For an MPS build where GPU is still unavailable:
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

