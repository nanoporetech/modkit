#!/bin/bash

# ================================================================================
# Modkit Compilation Script for macOS (Apple Silicon)
# ================================================================================
# This script automates the compilation of Oxford
# Nanopore's modkit bioinformatics tool on macOS with GPU acceleration support.
#
# Prerequisites: Apple Silicon Mac running macOS 12.3 or later
#
# What this script does:
#   0. Installs Xcode Command Line Tools
#   1. Installs Homebrew package manager
#   2. Installs rustup (Rust toolchain installer)
#   3. Installs Rust compiler (rustc) and Cargo build tool
#   4. Clones the modkit GitHub repository
#   5. Checks out the latest release version
#   6. Creates a Python virtual environment
#   7. Installs PyTorch in the virtual environment
#   8. Sets and verifies environment variables for libtorch
#   9. Builds modkit with macOS GPU (MPS) support
#
# Usage:
#   bash mac_compile_modkit.sh <installation_directory> [modkit_version]
#   bash mac_compile_modkit.sh --help
#
# Examples:
#   bash mac_compile_modkit.sh /path/to/install          # Install latest version to location
#   bash mac_compile_modkit.sh ~/tools v0.5.0            # Install specific version to location
#
# Optional Python controls:
#   MODKIT_PYTHON_PROVIDER=auto|system|pyenv|uv
#   MODKIT_PYTHON_VERSION=3.11.9
#   MODKIT_USE_UV=auto|0|1
#
# Examples:
#   MODKIT_PYTHON_PROVIDER=system bash mac_compile_modkit.sh ~/tools
#   MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash mac_compile_modkit.sh ~/tools
#   MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash mac_compile_modkit.sh ~/tools
#   MODKIT_USE_UV=0 bash mac_compile_modkit.sh ~/tools
# ================================================================================

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# ================================================================================
# Configuration
# ================================================================================

print_usage() {
    cat << EOF
Usage:
  bash $0 <installation_directory> [modkit_version]
  bash $0 --help

Arguments:
  installation_directory    Required. Path where modkit will be installed.
  modkit_version            Optional. Git tag to checkout. Default: latest.

Optional Python environment controls:
  MODKIT_PYTHON_PROVIDER    auto | system | pyenv | uv
                            Default: auto

  MODKIT_PYTHON_VERSION     Optional Python version request.
                            For uv:    3.11, 3.12, 3.12.2 are acceptable.
                            For pyenv: prefer an exact version, e.g. 3.11.9.

  MODKIT_USE_UV             auto | 0 | 1
                            auto: use uv for venv/pip if uv is available.
                            0:    use python -m venv and python -m pip.
                            1:    require/use uv venv and uv pip.

Examples:
  bash $0 ~/tools
  bash $0 ~/tools v0.5.0
  MODKIT_PYTHON_PROVIDER=system bash $0 ~/tools
  MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash $0 ~/tools
  MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash $0 ~/tools
  MODKIT_USE_UV=0 bash $0 ~/tools
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    print_usage
    exit 0
fi

# Installation directory mandatory
if [[ $# -lt 1 ]]; then
    print_usage
    exit 1
fi
INSTALL_DIR="$1"

# Modkit version to install default latest release
MODKIT_VERSION="${2:-latest}"

# Python toolchain controls
MODKIT_PYTHON_PROVIDER="${MODKIT_PYTHON_PROVIDER:-auto}"  # auto | system | pyenv | uv
MODKIT_PYTHON_VERSION="${MODKIT_PYTHON_VERSION:-}"        # optional
MODKIT_USE_UV="${MODKIT_USE_UV:-auto}"                    # auto | 0 | 1

# Resolved later by resolve_python_toolchain
MODKIT_PYTHON_BIN=""
MODKIT_PYTHON_PROVIDER_EFFECTIVE=""
MODKIT_UV_ENABLED=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ================================================================================
# Helper Functions
# ================================================================================

print_step() {
    echo ""
    echo -e "${BLUE}===================================================================${NC}"
    echo -e "${BLUE}STEP $1: $2${NC}"
    echo -e "${BLUE}===================================================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ Error: $1${NC}"
}

command_exists() {
    command -v "$1" &> /dev/null
}

press_enter_to_continue() {
    echo ""
    read -p "Press Enter to continue..."
}

ensure_command_or_brew_install() {
    local CMD="$1"
    local FORMULA="${2:-$1}"

    if command_exists "${CMD}"; then
        return 0
    fi

    if command_exists brew; then
        echo "Installing ${FORMULA} via Homebrew..."
        brew install "${FORMULA}"
    else
        print_error "${CMD} is required but was not found, and Homebrew is unavailable"
        exit 1
    fi
}

resolve_python_toolchain() {
    echo ""
    echo "Resolving Python toolchain..."

    case "${MODKIT_PYTHON_PROVIDER}" in
        auto|system|pyenv|uv) ;;
        *)
            print_error "Invalid MODKIT_PYTHON_PROVIDER='${MODKIT_PYTHON_PROVIDER}'. Use auto, system, pyenv, or uv."
            exit 1
            ;;
    esac

    case "${MODKIT_PYTHON_PROVIDER}" in
        system)
            MODKIT_PYTHON_PROVIDER_EFFECTIVE="system"
            MODKIT_PYTHON_BIN="$(command -v python3 || true)"
            ;;

        pyenv)
            ensure_command_or_brew_install pyenv
            MODKIT_PYTHON_PROVIDER_EFFECTIVE="pyenv"
            if [[ -n "${MODKIT_PYTHON_VERSION}" ]]; then
                pyenv install -s "${MODKIT_PYTHON_VERSION}"
                local PYENV_PREFIX
                PYENV_PREFIX="$(pyenv prefix "${MODKIT_PYTHON_VERSION}")"
                MODKIT_PYTHON_BIN="${PYENV_PREFIX}/bin/python3"
            else
                MODKIT_PYTHON_BIN="$(pyenv which python3)"
            fi
            ;;

        uv)
            ensure_command_or_brew_install uv
            MODKIT_PYTHON_PROVIDER_EFFECTIVE="uv"
            if [[ -n "${MODKIT_PYTHON_VERSION}" ]]; then
                uv python install "${MODKIT_PYTHON_VERSION}"
                MODKIT_PYTHON_BIN="$(uv python find "${MODKIT_PYTHON_VERSION}")"
            else
                MODKIT_PYTHON_BIN="$(uv python find python3 2>/dev/null || uv python find)"
            fi
            ;;

        auto)
            if [[ -n "${MODKIT_PYTHON_VERSION}" && $(command_exists uv && echo yes || echo no) == "yes" ]]; then
                MODKIT_PYTHON_PROVIDER_EFFECTIVE="uv"
                uv python install "${MODKIT_PYTHON_VERSION}"
                MODKIT_PYTHON_BIN="$(uv python find "${MODKIT_PYTHON_VERSION}")"
            elif [[ -n "${MODKIT_PYTHON_VERSION}" && $(command_exists pyenv && echo yes || echo no) == "yes" ]]; then
                MODKIT_PYTHON_PROVIDER_EFFECTIVE="pyenv"
                pyenv install -s "${MODKIT_PYTHON_VERSION}"
                local PYENV_PREFIX
                PYENV_PREFIX="$(pyenv prefix "${MODKIT_PYTHON_VERSION}")"
                MODKIT_PYTHON_BIN="${PYENV_PREFIX}/bin/python3"
            else
                # Suggestion 1 fix: warn if version was requested but cannot be honoured
                if [[ -n "${MODKIT_PYTHON_VERSION}" ]]; then
                    print_warning "MODKIT_PYTHON_VERSION='${MODKIT_PYTHON_VERSION}' was requested but neither uv nor pyenv is available."
                    print_warning "Falling back to system python3. The installed Python may not match the requested version."
                    read -p "Continue anyway? (y/n): " -n 1 -r
                    echo ""
                    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                        echo "Aborted. Install uv (brew install uv) or pyenv (brew install pyenv) and retry."
                        exit 1
                    fi
                fi
                MODKIT_PYTHON_PROVIDER_EFFECTIVE="system"
                MODKIT_PYTHON_BIN="$(command -v python3 || true)"
            fi
            ;;
    esac

    if [[ -z "${MODKIT_PYTHON_BIN}" || ! -x "${MODKIT_PYTHON_BIN}" ]]; then
        print_error "Could not resolve a usable Python executable"
        exit 1
    fi

    case "${MODKIT_USE_UV}" in
        1|true|yes)
            ensure_command_or_brew_install uv
            MODKIT_UV_ENABLED=1
            ;;
        0|false|no)
            MODKIT_UV_ENABLED=0
            ;;
        auto)
            if command_exists uv; then
                MODKIT_UV_ENABLED=1
            else
                MODKIT_UV_ENABLED=0
            fi
            ;;
        *)
            print_error "Invalid MODKIT_USE_UV='${MODKIT_USE_UV}'. Use auto, 0, or 1."
            exit 1
            ;;
    esac

    echo "Python configuration:"
    echo "  Provider requested: ${MODKIT_PYTHON_PROVIDER}"
    echo "  Provider used:      ${MODKIT_PYTHON_PROVIDER_EFFECTIVE}"
    echo "  Version requested:  ${MODKIT_PYTHON_VERSION:-default}"
    echo "  Python executable:  ${MODKIT_PYTHON_BIN}"
    echo "  Python version:     $("${MODKIT_PYTHON_BIN}" --version 2>&1)"
    echo "  uv for venv/pip:    $([[ "${MODKIT_UV_ENABLED}" == "1" ]] && echo yes || echo no)"
}

# ================================================================================
# STEP 0: Install Xcode Command Line Tools
# ================================================================================

install_xcode_tools() {
    print_step 0 "Installing Xcode Command Line Tools"
    
    if xcode-select -p &> /dev/null; then
        print_success "Xcode Command Line Tools already installed at: $(xcode-select -p)"
    else
        echo "Installing Xcode Command Line Tools..."
        echo "A dialog will appear - please click 'Install' and wait for completion."
        xcode-select --install
        
        echo ""
        echo "Waiting for Xcode Command Line Tools installation to complete..."
        echo "This may take several minutes (timeout: 10 minutes)."
        
        # Wait for installation to complete with a timeout
        local MAX_WAIT=600  # 10 minutes in seconds
        local WAITED=0
        while ! xcode-select -p &> /dev/null; do
            sleep 5
            WAITED=$((WAITED + 5))
            if [[ ${WAITED} -ge ${MAX_WAIT} ]]; then
                print_error "Timed out waiting for Xcode Command Line Tools installation"
                echo "Please try one of the following:"
                echo "  1. Re-run: xcode-select --install"
                echo "  2. Install via: System Settings > Software Update"
                echo "  3. Download from: https://developer.apple.com/download/more/"
                echo "Then re-run this script."
                exit 1
            fi
        done
        
        print_success "Xcode Command Line Tools installed successfully"
    fi
    
    # Verify installation
    if xcode-select -p &> /dev/null; then
        print_success "Verification: Xcode tools are available"
    else
        print_error "Xcode Command Line Tools installation failed"
        exit 1
    fi
}

# ================================================================================
# STEP 1: Install Homebrew
# ================================================================================

install_homebrew() {
    print_step 1 "Installing Homebrew Package Manager"
    
    if command_exists brew; then
        print_success "Homebrew already installed at: $(which brew)"
        echo "Current version: $(brew --version | head -n1)"
    else
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        
        # Add Homebrew to PATH based on architecture (only if not already in .zprofile)
        if [[ $(uname -m) == "arm64" ]]; then
            # Apple Silicon
            local BREW_LINE='eval "$(/opt/homebrew/bin/brew shellenv)"'
            if ! grep -qF "${BREW_LINE}" "${HOME}/.zprofile" 2>/dev/null; then
                echo "${BREW_LINE}" >> "${HOME}/.zprofile"
            fi
            eval "$(/opt/homebrew/bin/brew shellenv)"
        else
            # Intel
            local BREW_LINE='eval "$(/usr/local/bin/brew shellenv)"'
            if ! grep -qF "${BREW_LINE}" "${HOME}/.zprofile" 2>/dev/null; then
                echo "${BREW_LINE}" >> "${HOME}/.zprofile"
            fi
            eval "$(/usr/local/bin/brew shellenv)"
        fi
        
        print_success "Homebrew installed successfully"
    fi
    
    # Verify installation
    if command_exists brew; then
        print_success "Verification: Homebrew is available"
        brew --version
    else
        print_error "Homebrew installation failed"
        exit 1
    fi
}

# ================================================================================
# STEP 2: Install rustup using Homebrew
# ================================================================================

install_rustup() {
    print_step 2 "Installing rustup (Rust Toolchain Installer)"
    
    if command_exists rustup; then
        print_success "rustup already installed at: $(which rustup)"
        echo "Current version: $(rustup --version)"
    elif command_exists rustup-init; then
        print_success "rustup-init already available at: $(which rustup-init)"
        echo "Will be used in the next step to install the Rust toolchain."
    else
        echo "Installing rustup-init via Homebrew..."
        brew install rustup-init
        print_success "rustup-init installed via Homebrew"
    fi
}

# ================================================================================
# STEP 3: Install Rust and Cargo using rustup
# ================================================================================

install_rust_cargo() {
    print_step 3 "Installing Rust Compiler and Cargo Build Tool"
    
    if command_exists rustc && command_exists cargo; then
        print_success "Rust and Cargo already installed"
        echo "  rustc version: $(rustc --version)"
        echo "  cargo version: $(cargo --version)"
    else
        echo "Running rustup-init to install Rust toolchain..."
        rustup-init -y --default-toolchain stable
        
        # Source cargo environment
        source "${HOME}/.cargo/env"
        
        print_success "Rust toolchain installed successfully"
    fi
    
    # Ensure cargo is in PATH for current session
    if [[ -f "${HOME}/.cargo/env" ]]; then
        source "${HOME}/.cargo/env"
    fi
    
    # Verify installation
    if command_exists rustc && command_exists cargo; then
        print_success "Verification: Rust and Cargo are available"
        echo "  rustc: $(which rustc) - $(rustc --version)"
        echo "  cargo: $(which cargo) - $(cargo --version)"
    else
        print_error "Rust/Cargo installation failed"
        echo "Please run: source ${HOME}/.cargo/env"
        exit 1
    fi
}

# ================================================================================
# STEP 4: Clone Modkit GitHub Repository
# ================================================================================

clone_modkit_repo() {
    print_step 4 "Cloning Modkit GitHub Repository"

    mkdir -p "${INSTALL_DIR}"
    cd "${INSTALL_DIR}"

    MODKIT_REPO_DIR="${INSTALL_DIR}/modkit"

    if [[ -d "${MODKIT_REPO_DIR}/.git" ]]; then
        print_warning "Modkit repository already exists at: ${MODKIT_REPO_DIR}"
        echo "Updating existing repository..."
        cd "${MODKIT_REPO_DIR}"
        git fetch --all --tags
        # Reset any stale working tree changes so the subsequent checkout is clean
        git reset --hard
        git clean -fd
        print_success "Repository updated and working tree cleaned"
    else
        echo "Cloning modkit repository..."
        git clone https://github.com/nanoporetech/modkit.git
        cd modkit
        print_success "Repository cloned successfully"
    fi

    if [[ -d "${MODKIT_REPO_DIR}/.git" ]]; then
        print_success "Verification: Repository available at ${MODKIT_REPO_DIR}"
    else
        print_error "Failed to clone modkit repository"
        exit 1
    fi
}

# ================================================================================
# STEP 5: Checkout Latest Release (or specified version)
# ================================================================================

checkout_version() {
    print_step 5 "Checking Out Modkit Version"

    cd "${MODKIT_REPO_DIR}"

    if [[ "${MODKIT_VERSION}" == "latest" ]]; then
        echo "Fetching latest release version..."
        LATEST_TAG=$(git tag -l "v*" --sort=-v:refname 2>/dev/null | head -n1 || echo "")

        if [[ -z "${LATEST_TAG}" ]]; then
            LATEST_TAG=$(git tag --sort=-v:refname 2>/dev/null | head -n1 || echo "")
        fi

        if [[ -z "${LATEST_TAG}" ]]; then
            print_warning "No release tags found, using main branch"
            git checkout main
            git pull origin main
            MODKIT_VERSION="main"
        else
            echo "Latest release: ${LATEST_TAG}"
            git checkout "${LATEST_TAG}"
            MODKIT_VERSION="${LATEST_TAG}"
        fi
    else
        echo "Checking out version: ${MODKIT_VERSION}"
        git checkout "${MODKIT_VERSION}"
    fi

    print_success "Using modkit version: ${MODKIT_VERSION}"
    echo "Current commit: $(git rev-parse --short HEAD)"

    # Verify Cargo.toml version matches the git tag as a sanity check
    local CARGO_VERSION
    CARGO_VERSION=$(grep -m1 '^version' "${MODKIT_REPO_DIR}/modkit/Cargo.toml" 2>/dev/null \
        | sed 's/.*"\(.*\)".*/\1/' || echo "unknown")
    echo "Cargo crate version: ${CARGO_VERSION}"

    if [[ "${MODKIT_VERSION}" != "main" && \
          "${MODKIT_VERSION}" != "v${CARGO_VERSION}" && \
          "${MODKIT_VERSION}" != "${CARGO_VERSION}" ]]; then
        print_warning "Git tag '${MODKIT_VERSION}' does not match Cargo.toml version '${CARGO_VERSION}'."
        print_warning "The compiled binary will report version ${CARGO_VERSION}."
        print_warning "This is an upstream discrepancy in the modkit repository."
    else
        print_success "Git tag matches Cargo.toml version: ${CARGO_VERSION}"
    fi
}

# ================================================================================
# STEP 6: Create Python Virtual Environment
# ================================================================================

create_venv() {
    print_step 6 "Creating Python Virtual Environment"

    VENV_DIR="${INSTALL_DIR}/venv_modkit"
    local FINGERPRINT_FILE="${VENV_DIR}/.python_info"

    echo "Using Python: ${MODKIT_PYTHON_BIN} - $("${MODKIT_PYTHON_BIN}" --version 2>&1)"

    # Suggestion 2 fix: validate existing venv against current provider/version
    if [[ -d "${VENV_DIR}" ]]; then
        if [[ -f "${FINGERPRINT_FILE}" ]]; then
            local SAVED_PROVIDER SAVED_VERSION
            SAVED_PROVIDER=$(grep '^provider=' "${FINGERPRINT_FILE}" | cut -d= -f2)
            SAVED_VERSION=$(grep '^python_bin=' "${FINGERPRINT_FILE}" | cut -d= -f2)
            if [[ "${SAVED_PROVIDER}" != "${MODKIT_PYTHON_PROVIDER_EFFECTIVE}" || \
                  "${SAVED_VERSION}" != "${MODKIT_PYTHON_BIN}" ]]; then
                print_warning "Existing venv was created with a different Python configuration:"
                print_warning "  Saved provider:    ${SAVED_PROVIDER}  (current: ${MODKIT_PYTHON_PROVIDER_EFFECTIVE})"
                print_warning "  Saved python_bin:  ${SAVED_VERSION}"
                print_warning "  Current python_bin: ${MODKIT_PYTHON_BIN}"
                read -p "Delete and recreate the venv with the new configuration? (y/n): " -n 1 -r
                echo ""
                if [[ $REPLY =~ ^[Yy]$ ]]; then
                    rm -rf "${VENV_DIR}"
                    echo "Deleted old virtual environment."
                else
                    print_warning "Keeping existing venv. Build may use wrong Python."
                fi
            else
                print_success "Existing venv matches current configuration. Reusing."
            fi
        else
            print_warning "Virtual environment exists but has no fingerprint file."
            print_warning "Cannot verify it matches the current Python configuration. Reusing."
        fi
    fi

    if [[ ! -d "${VENV_DIR}" ]]; then
        echo "Creating virtual environment at: ${VENV_DIR}"
        if [[ "${MODKIT_UV_ENABLED}" == "1" ]]; then
            uv venv --python "${MODKIT_PYTHON_BIN}" "${VENV_DIR}"
        else
            "${MODKIT_PYTHON_BIN}" -m venv "${VENV_DIR}"
        fi
        # Write fingerprint
        printf 'provider=%s\npython_bin=%s\n' \
            "${MODKIT_PYTHON_PROVIDER_EFFECTIVE}" \
            "${MODKIT_PYTHON_BIN}" > "${FINGERPRINT_FILE}"
        print_success "Virtual environment created"
    fi

    # Verify virtual environment
    if [[ -f "${VENV_DIR}/bin/activate" && -x "${VENV_DIR}/bin/python" ]]; then
        print_success "Verification: Virtual environment ready"
        echo "  Venv Python: ${VENV_DIR}/bin/python - $("${VENV_DIR}/bin/python" --version 2>&1)"
    else
        print_error "Failed to create virtual environment"
        exit 1
    fi
}

# ================================================================================
# STEP 7: Activate Environment and Install PyTorch
# ================================================================================

install_pytorch() {
    print_step 7 "Installing PyTorch in Virtual Environment"
    
    local VENV_PYTHON="${VENV_DIR}/bin/python"

    # Deactivate conda only if an environment is actually active.
    if [[ -n "${CONDA_DEFAULT_ENV:-}" ]] && command_exists conda; then
        print_warning "Conda environment '${CONDA_DEFAULT_ENV}' is active. Deactivating to avoid conflicts..."
        eval "$(conda shell.bash hook)"
        conda deactivate 2>/dev/null || true
        print_success "Conda environment deactivated"
    fi
    
    # Activate virtual environment
    echo "Activating virtual environment..."
    source "${VENV_DIR}/bin/activate"
    
    print_success "Virtual environment activated: ${VIRTUAL_ENV}"
    
    # Install PyTorch
    echo "Installing PyTorch and NumPy and dependencies (this may take a few minutes)..."
    if [[ "${MODKIT_UV_ENABLED}" == "1" ]]; then
        uv pip install --python "${VENV_PYTHON}" torch numpy
    else
        "${VENV_PYTHON}" -m pip install --upgrade pip
        "${VENV_PYTHON}" -m pip install torch numpy
    fi
    
    print_success "PyTorch and NumPy installed successfully"
    
    # Verify PyTorch installation
    echo "Verifying PyTorch installation..."
    "${VENV_PYTHON}" -c "import torch; print(f'PyTorch version: {torch.__version__}')" || {
        print_error "PyTorch verification failed"
        exit 1
    }
    
    # Check MPS (Metal Performance Shaders) availability for Apple Silicon
    if [[ $(uname -m) == "arm64" ]]; then
        echo "Checking Metal Performance Shaders (GPU) support..."
        "${VENV_PYTHON}" -c "import torch; print(f'MPS available: {torch.backends.mps.is_available()}')"
    fi
    
    print_success "Verification: PyTorch is working correctly"
}

# ================================================================================
# STEP 8: Set and Verify Environment Variables
# ================================================================================

setup_environment_variables() {
    print_step 8 "Setting Up Environment Variables for libtorch"

    SETUP_SCRIPT="${INSTALL_DIR}/setup_modkit_env.sh"

    echo "Generating environment setup script at: ${SETUP_SCRIPT}"
    cat > "${SETUP_SCRIPT}" << 'SETUP_EOF'
#!/bin/bash
# ================================================================================
# Modkit Environment Setup Script
# ================================================================================
# IMPORTANT: This script must be SOURCED, not executed.
#
# Usage:
#   source setup_modkit_env.sh <installation_directory> [--verbose]
#
# Options:
#   --verbose    Print environment configuration details. Omit for silent setup.
#
# After sourcing:
#   modkit --version
#   modkit pileup input.bam output.bed
#   modkit open-chromatin predict --device mps -i input.bam -o output.bed
# ================================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced, not executed."
    echo "Usage: source ${0} <installation_directory> [--verbose]"
    exit 1
fi

# Parse arguments
_MODKIT_VERBOSE=0
_MODKIT_INSTALL_DIR=""
for _arg in "$@"; do
    case "${_arg}" in
        --verbose|-v) _MODKIT_VERBOSE=1 ;;
        *)            _MODKIT_INSTALL_DIR="${_arg}" ;;
    esac
done
unset _arg

if [[ -z "${_MODKIT_INSTALL_DIR}" ]]; then
    echo "Usage: source setup_modkit_env.sh <installation_directory> [--verbose]"
    unset _MODKIT_VERBOSE _MODKIT_INSTALL_DIR
    return 1
fi

MODKIT_INSTALL_DIR="${_MODKIT_INSTALL_DIR}"
unset _MODKIT_INSTALL_DIR

if [[ ! -d "${MODKIT_INSTALL_DIR}" ]]; then
    echo "Error: Installation directory not found: ${MODKIT_INSTALL_DIR}"
    unset MODKIT_INSTALL_DIR _MODKIT_VERBOSE
    return 1
fi

MODKIT_VENV_DIR="${MODKIT_INSTALL_DIR}/venv_modkit"
MODKIT_REPO_DIR="${MODKIT_INSTALL_DIR}/modkit"
MODKIT_BINARY="${MODKIT_REPO_DIR}/target/release/modkit"

if [[ ! -d "${MODKIT_VENV_DIR}" ]]; then
    echo "Error: Virtual environment not found: ${MODKIT_VENV_DIR}"
    unset MODKIT_INSTALL_DIR MODKIT_VENV_DIR MODKIT_REPO_DIR MODKIT_BINARY _MODKIT_VERBOSE
    return 1
fi

# Deactivate conda only if actually active
if command -v conda &> /dev/null && [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    [[ "${_MODKIT_VERBOSE}" == "1" ]] && \
        echo "Deactivating conda environment '${CONDA_DEFAULT_ENV}'..."
    eval "$(conda shell.bash hook)"
    conda deactivate 2>/dev/null || true
fi

# Activate the Python virtual environment (suppress output)
source "${MODKIT_VENV_DIR}/bin/activate"

# Detect Python version from venv
MODKIT_PYTHON_VER=$("${MODKIT_VENV_DIR}/bin/python" -c \
    "import sys; print(f'python{sys.version_info.major}.{sys.version_info.minor}')")

export LIBTORCH_USE_PYTORCH=1
export LIBTORCH_BYPASS_VERSION_CHECK=1
export LIBTORCH="${MODKIT_VENV_DIR}/lib/${MODKIT_PYTHON_VER}/site-packages/torch"
export DYLD_LIBRARY_PATH="${LIBTORCH}/lib${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
export LD_LIBRARY_PATH="${LIBTORCH}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PATH="${MODKIT_REPO_DIR}/target/release:${PATH}"

_PERF_CORES=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
export RAYON_NUM_THREADS="${_PERF_CORES}"
unset _PERF_CORES

if [[ "${_MODKIT_VERBOSE}" == "1" ]]; then
    echo ""
    echo "Modkit environment configured:"
    echo "  LIBTORCH_USE_PYTORCH          = ${LIBTORCH_USE_PYTORCH}"
    echo "  LIBTORCH_BYPASS_VERSION_CHECK = ${LIBTORCH_BYPASS_VERSION_CHECK}"
    echo "  LIBTORCH                      = ${LIBTORCH}"
    echo "  DYLD_LIBRARY_PATH             = ${DYLD_LIBRARY_PATH}"
    echo "  LD_LIBRARY_PATH               = ${LD_LIBRARY_PATH}"
    echo "  PATH (prepended)              = ${MODKIT_REPO_DIR}/target/release"
    echo "  RAYON_NUM_THREADS             = ${RAYON_NUM_THREADS}  (P-cores)"
    echo "  Python                        = ${MODKIT_VENV_DIR}/bin/python" \
         "($("${MODKIT_VENV_DIR}/bin/python" --version 2>&1), ${MODKIT_PYTHON_VER})"

    if [[ -d "${LIBTORCH}" ]]; then
        echo "  libtorch path                 = OK"
    else
        echo "  WARNING: LIBTORCH path does not exist: ${LIBTORCH}"
    fi

    if [[ -f "${MODKIT_BINARY}" ]]; then
        echo "  modkit binary                 = ${MODKIT_BINARY}"
        echo ""
        echo "Ready. Run 'modkit --version' to verify."
    else
        echo "  modkit binary                 = NOT FOUND"
    fi
fi

unset MODKIT_INSTALL_DIR MODKIT_VENV_DIR MODKIT_REPO_DIR MODKIT_BINARY MODKIT_PYTHON_VER _MODKIT_VERBOSE
SETUP_EOF

    chmod +x "${SETUP_SCRIPT}"
    print_success "Environment setup script generated at: ${SETUP_SCRIPT}"

    echo "Sourcing environment setup script..."
    local SAVED_MODKIT_REPO_DIR="${MODKIT_REPO_DIR}"
    local SAVED_VENV_DIR="${VENV_DIR}"

    # When sourcing during install, use --verbose so the installer still shows full output
    source "${SETUP_SCRIPT}" "${INSTALL_DIR}" --verbose

    MODKIT_REPO_DIR="${SAVED_MODKIT_REPO_DIR}"
    VENV_DIR="${SAVED_VENV_DIR}"

    echo ""
    echo "Verifying paths..."

    if [[ -d "${LIBTORCH}" ]]; then
        print_success "LIBTORCH path exists: ${LIBTORCH}"
    else
        print_error "LIBTORCH path not found: ${LIBTORCH}"
        exit 1
    fi

    if [[ -d "${LIBTORCH}/lib" ]]; then
        print_success "libtorch libraries found"
        echo "Sample libraries:"
        ls "${LIBTORCH}/lib" | grep -E "\.dylib$" | head -5
    else
        print_error "libtorch lib directory not found"
        exit 1
    fi

    echo ""
    echo "To set up this environment in a new terminal session, run:"
    echo "  source \"${SETUP_SCRIPT}\" \"${INSTALL_DIR}\""
}

# ================================================================================
# STEP 9: Build Modkit with macOS GPU Support (continued)
# ================================================================================

build_modkit() {
    print_step 9 "Building Modkit with macOS GPU (MPS) Support"
    
    cd "${MODKIT_REPO_DIR}"
    
    echo "Build configuration:"
    echo "  Features: accelerate,tch"
    echo "  Mode: release (optimized)"
    echo "  GPU Support: Metal Performance Shaders (MPS)"
    echo "  Architecture: $(uname -m)"
    echo ""
    echo "This will take several minutes (5-15 minutes depending on your Mac)..."
    echo "Please be patient while Cargo compiles modkit and its dependencies."
    echo ""
    
    # Run cargo build
    if cargo build --release --features accelerate,tch; then
        print_success "Modkit compiled successfully!"
    else
        print_error "Compilation failed"
        echo ""
        echo "Troubleshooting tips:"
        echo "  1. Set up the build environment first:"
        echo "       source \"${INSTALL_DIR}/setup_modkit_env.sh\" \"${INSTALL_DIR}\""
        echo "  2. Check that PyTorch is installed:"
        echo "       \"${VENV_DIR}/bin/python\" -m pip list | grep torch"
        echo "  3. Clean and rebuild:"
        echo "       cd \"${MODKIT_REPO_DIR}\" && cargo clean && cargo build --release --features accelerate,tch"
        exit 1
    fi
    
    # Verify binary was created
    MODKIT_BINARY="${MODKIT_REPO_DIR}/target/release/modkit"
    
    if [[ -f "${MODKIT_BINARY}" ]]; then
        print_success "Binary created at: ${MODKIT_BINARY}"
        
        # Get file size
        BINARY_SIZE=$(du -h "${MODKIT_BINARY}" | cut -f1)
        echo "  Binary size: ${BINARY_SIZE}"
        
        # Test binary
        echo ""
        echo "Testing modkit binary..."
        if "${MODKIT_BINARY}" --version; then
            print_success "Modkit is working correctly!"
        else
            print_warning "Binary exists but version check failed"
            echo "You may need to set DYLD_LIBRARY_PATH when running modkit"
        fi
        
        # Check available subcommands
        echo ""
        echo "Available modkit commands:"
        "${MODKIT_BINARY}" --help | grep -A 20 "SUBCOMMANDS:" || true
        
    else
        print_error "Binary not found at expected location: ${MODKIT_BINARY}"
        exit 1
    fi
}

# ================================================================================
# Save Installation Info
# ================================================================================

save_installation_info() {
    INFO_FILE="${INSTALL_DIR}/installation_info.txt"
    SETUP_SCRIPT="${INSTALL_DIR}/setup_modkit_env.sh"
    local UV_VERSION_INFO
    local PYENV_VERSION_INFO
    local VENV_PYTHON="${VENV_DIR}/bin/python"

    UV_VERSION_INFO="$(command_exists uv && uv --version || echo "Not available")"
    PYENV_VERSION_INFO="$(command_exists pyenv && pyenv --version || echo "Not available")"
    
    cat > "${INFO_FILE}" << EOF
 Modkit Installation Information
 ================================

 Installation Date: $(date)
 macOS Version: $(sw_vers -productVersion)
 Architecture: $(uname -m)
 Hostname: $(hostname)

 Installation Paths:
 -------------------
 Installation Directory: ${INSTALL_DIR}
 Modkit Repository: ${MODKIT_REPO_DIR}
 Modkit Binary: ${MODKIT_REPO_DIR}/target/release/modkit
 Virtual Environment: ${VENV_DIR}
 Environment Setup Script: ${SETUP_SCRIPT}

 Versions:
 ---------
 Modkit Version: ${MODKIT_VERSION}
 Rust Version: $(rustc --version)
 Cargo Version: $(cargo --version)
 Python Provider Requested: ${MODKIT_PYTHON_PROVIDER}
 Python Provider Used: ${MODKIT_PYTHON_PROVIDER_EFFECTIVE}
 Python Executable: ${MODKIT_PYTHON_BIN}
 Python Version: $("${MODKIT_PYTHON_BIN}" --version 2>&1)
 Virtual Environment Python: ${VENV_PYTHON}
 Virtual Environment Python Version: $("${VENV_PYTHON}" --version 2>&1)
 uv Enabled For Venv/Pip: $([[ "${MODKIT_UV_ENABLED}" == "1" ]] && echo "yes" || echo "no")
 uv Version: ${UV_VERSION_INFO}
 pyenv Version: ${PYENV_VERSION_INFO}
 PyTorch Version: $("${VENV_PYTHON}" -c "import torch; print(torch.__version__)" 2>/dev/null || echo "Unknown")

 Quick Start:
 -----------
 source "${SETUP_SCRIPT}" "${INSTALL_DIR}"
 modkit --version

 Python Run Examples:
 --------------------
 System Python:
 MODKIT_PYTHON_PROVIDER=system bash mac_compile_modkit.sh "${INSTALL_DIR}" "${MODKIT_VERSION}"

 pyenv Python:
 MODKIT_PYTHON_PROVIDER=pyenv MODKIT_PYTHON_VERSION=3.11.9 bash mac_compile_modkit.sh "${INSTALL_DIR}" "${MODKIT_VERSION}"

 uv Python:
 MODKIT_PYTHON_PROVIDER=uv MODKIT_PYTHON_VERSION=3.12 bash mac_compile_modkit.sh "${INSTALL_DIR}" "${MODKIT_VERSION}"

 Disable uv venv/pip usage:
 MODKIT_USE_UV=0 bash mac_compile_modkit.sh "${INSTALL_DIR}" "${MODKIT_VERSION}"

 For detailed usage instructions, run:
 modkit --help

 Runtime Threading:
 ------------------
 RAYON_NUM_THREADS is set automatically to the number of Performance cores
 detected on your Mac (hw.perflevel0.logicalcpu). This is exported by
 setup_modkit_env.sh each time it is sourced.

 To override for a single run:
   RAYON_NUM_THREADS=8 modkit pileup input.bam output.bed
EOF
    
    echo "Installation information saved to: ${INFO_FILE}"
}

# ================================================================================
# Main Installation Flow
# ================================================================================

main() {
    echo ""
    echo -e "${BLUE}###################################################################${NC}"
    echo -e "${BLUE}#                                                                 #${NC}"
    echo -e "${BLUE}#        Modkit Installation Script for macOS                     #${NC}"
    echo -e "${BLUE}#                                                                 #${NC}"
    echo -e "${BLUE}###################################################################${NC}"
    echo ""
    echo "This script will install and compile modkit with GPU support."
    echo ""
    echo "Installation directory: ${INSTALL_DIR}"
    echo "Modkit version: ${MODKIT_VERSION}"
    echo "Python provider: ${MODKIT_PYTHON_PROVIDER}"
    echo "Python version request: ${MODKIT_PYTHON_VERSION:-default}"
    echo "Use uv for venv/pip: ${MODKIT_USE_UV}"
    echo ""

    # Confirm before proceeding
    read -p "Continue with installation? (y/n): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled."
        exit 0
    fi
    
    # Record start time
    START_TIME=$(date +%s)
    
    # Execute installation steps
    install_xcode_tools
    install_homebrew
    install_rustup
    install_rust_cargo
    clone_modkit_repo
    checkout_version
    resolve_python_toolchain
    create_venv
    install_pytorch
    setup_environment_variables
    build_modkit
    
    # Save installation info
    save_installation_info

    # Add setup script to ~/.zprofile if not already present
    local ZPROFILE="${HOME}/.zprofile"
    local SOURCE_LINE="source \"${SETUP_SCRIPT}\" \"${INSTALL_DIR}\""
    if grep -qF "${SOURCE_LINE}" "${ZPROFILE}" 2>/dev/null; then
        print_success "~/.zprofile already contains the modkit environment setup"
    else
        echo "${SOURCE_LINE}" >> "${ZPROFILE}"
        print_success "Added modkit environment setup to ~/.zprofile"
        echo "  It will activate automatically in every new terminal session."
    fi

    # Calculate installation time
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    MINUTES=$((DURATION / 60))
    REMAINING_SECONDS=$((DURATION % 60))
    
    echo ""
    echo -e "${GREEN}Total installation time: ${MINUTES} minutes ${REMAINING_SECONDS} seconds${NC}"
    echo ""
    
}

# ================================================================================
# Script Entry Point
# ================================================================================

# Check if running on macOS
if [[ "$(uname -s)" != "Darwin" ]]; then
    print_error "This script is designed for macOS only"
    echo "Detected OS: $(uname -s)"
    exit 1
fi

# Check macOS version (require 11.0+)
MACOS_VERSION=$(sw_vers -productVersion)
echo "Detected macOS version: ${MACOS_VERSION}"

MACOS_MAJOR=$(echo "${MACOS_VERSION}" | cut -d. -f1)
MACOS_MINOR=$(echo "${MACOS_VERSION}" | cut -d. -f2)

if [[ "${MACOS_MAJOR}" -lt 12 ]] || \
   [[ "${MACOS_MAJOR}" -eq 12 && "${MACOS_MINOR}" -lt 3 ]]; then
    print_error "macOS 12.3 (Monterey) or later is required for PyTorch Metal (MPS) GPU support."
    echo "Detected macOS version: ${MACOS_VERSION}"
    echo "On macOS < 12.3 torch.backends.mps.is_available() returns false and GPU acceleration"
    echo "cannot be used. Please upgrade macOS before running this script."
    exit 1
fi

# Run main installation
main

exit 0
