# Modkit Installation for macOS

**Prerequisites:** Apple Silicon Mac, macOS 11+

## Quick Start

```bash
bash mac_compile_modkit.sh ~/tools
```

Installs modkit to `~/tools` with default settings (system Python, latest version). Takes 10–15 minutes.

---

## What Gets Installed

- Xcode Command Line Tools
- Homebrew
- Rust & Cargo
- Python virtual environment
- PyTorch with GPU (Metal Performance Shaders)
- Compiled modkit binary

---

## After Installation

Every new terminal session:

```bash
source ~/tools/setup_modkit_env.sh ~/tools
modkit --version
```

Add to `~/.zprofile` for auto-setup.

---

## Python Version Control

By default, uses system `python3`. To control which Python:

| Use Case | Command |
|----------|---------|
| System Python (default) | `bash mac_compile_modkit.sh ~/tools` |
| Specific pyenv version | `MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash mac_compile_modkit.sh ~/tools` |
| Specific uv version | `MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash mac_compile_modkit.sh ~/tools` |
| Without uv (slower) | `MODKIT_USE_UV=0 bash mac_compile_modkit.sh ~/tools` |

### Environment Variables

- `MODKIT_PYTHON_PROVIDER`: `auto` (default), `system`, `pyenv`, `uv`
- `MODKIT_PYTHON_VERSION`: e.g., `3.11.9` or `3.12` (optional)
- `MODKIT_USE_UV`: `auto` (default), `0`, `1`

---

## Common Issues & Fixes

### "Could not resolve a usable Python executable"

Check Python is installed:
```bash
python3 --version
```

If missing, install via Homebrew:
```bash
brew install python@3.11
```

### "PyTorch verification failed"

Reinstall PyTorch:
```bash
source ~/tools/setup_modkit_env.sh ~/tools
~/tools/venv_modkit/bin/python -m pip install --upgrade pip
~/tools/venv_modkit/bin/python -m pip install torch numpy
```

### "Compilation failed"

Clean and rebuild:
```bash
cd ~/tools/modkit
cargo clean
LIBTORCH_USE_PYTORCH=1 cargo build --release --features accelerate,tch
```

###  "MPS available: false"

GPU not detected. Check:
- Apple Silicon: `uname -m` → should print `arm64`
- macOS 11+: `sw_vers -productVersion`
- CPU-only fallback is used automatically

---

## All Options

```bash
bash mac_compile_modkit.sh --help
```

---

## Installation Info

After install, view details:
```bash
cat ~/tools/installation_info.txt
```

Contains paths, versions, and reproducible run commands.

---

## Next Steps

- Read modkit help: `modkit --help`
- See [modkit repository](https://github.com/nanoporetech/modkit)

