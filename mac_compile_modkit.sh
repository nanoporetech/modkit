#!/bin/bash

################################################################################
# Modkit Compilation Script for macOS (Apple Silicon)
################################################################################
# This script automates the compilation of Oxford
# Nanopore's modkit bioinformatics tool on macOS with GPU acceleration support.
#
# Prerequisites: Apple Silicon Mac running macOS 11 or later
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
#
# Examples:
#   bash mac_compile_modkit.sh /path/to/install          # Install latest version to location
#   bash mac_compile_modkit.sh ~/tools v0.5.0            # Install specific version to location
#
################################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

################################################################################
# Configuration
################################################################################

# Installation directory (mandatory)
if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <installation_directory> [modkit_version]"
    echo "  installation_directory: Required. Path where modkit will be installed."
    echo "  modkit_version:         Optional. Git tag to checkout (default: latest)."
    exit 1
fi
INSTALL_DIR="$1"

# Modkit version to install (default: latest release)
MODKIT_VERSION="${2:-latest}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Helper Functions
################################################################################

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

################################################################################
# STEP 0: Install Xcode Command Line Tools
################################################################################

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

################################################################################
# STEP 1: Install Homebrew
################################################################################

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

################################################################################
# STEP 2: Install rustup using Homebrew
################################################################################

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

################################################################################
# STEP 3: Install Rust and Cargo using rustup
################################################################################

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

################################################################################
# STEP 4: Clone Modkit GitHub Repository
################################################################################

clone_modkit_repo() {
    print_step 4 "Cloning Modkit GitHub Repository"
    
    # Create installation directory
    mkdir -p "${INSTALL_DIR}"
    cd "${INSTALL_DIR}"
    
    MODKIT_REPO_DIR="${INSTALL_DIR}/modkit"
    
    if [[ -d "${MODKIT_REPO_DIR}/.git" ]]; then
        print_warning "Modkit repository already exists at: ${MODKIT_REPO_DIR}"
        echo "Updating existing repository..."
        cd "${MODKIT_REPO_DIR}"
        git fetch --all --tags
        print_success "Repository updated"
    else
        echo "Cloning modkit repository..."
        git clone https://github.com/nanoporetech/modkit.git
        cd modkit
        print_success "Repository cloned successfully"
    fi
    
    # Verify clone
    if [[ -d "${MODKIT_REPO_DIR}/.git" ]]; then
        print_success "Verification: Repository available at ${MODKIT_REPO_DIR}"
    else
        print_error "Failed to clone modkit repository"
        exit 1
    fi
}

################################################################################
# STEP 5: Checkout Latest Release (or specified version)
################################################################################

checkout_version() {
    print_step 5 "Checking Out Modkit Version"
    
    cd "${MODKIT_REPO_DIR}"
    
    if [[ "${MODKIT_VERSION}" == "latest" ]]; then
        # Get the latest release tag by semantic version (sorted by version, descending)
        echo "Fetching latest release version..."
        LATEST_TAG=$(git tag -l "v*" --sort=-v:refname 2>/dev/null | head -n1 || echo "")
        
        if [[ -z "${LATEST_TAG}" ]]; then
            # Fallback: try all tags if no v* tags exist
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
    
    # Show current commit
    echo "Current commit: $(git rev-parse --short HEAD)"
}

################################################################################
# STEP 6: Create Python Virtual Environment
################################################################################

create_venv() {
    print_step 6 "Creating Python Virtual Environment"
    
    VENV_DIR="${INSTALL_DIR}/venv_modkit"
    
    echo "Using Python: $(which python3) - $(python3 --version)"
    
    if [[ -d "${VENV_DIR}" ]]; then
        print_warning "Virtual environment already exists at: ${VENV_DIR}"
        echo "Using existing virtual environment"
    else
        echo "Creating virtual environment at: ${VENV_DIR}"
        python3 -m venv "${VENV_DIR}"
        print_success "Virtual environment created"
    fi
    
    # Verify virtual environment
    if [[ -f "${VENV_DIR}/bin/activate" ]]; then
        print_success "Verification: Virtual environment ready"
    else
        print_error "Failed to create virtual environment"
        exit 1
    fi
}

################################################################################
# STEP 7: Activate Environment and Install PyTorch
################################################################################

install_pytorch() {
    print_step 7 "Installing PyTorch in Virtual Environment"
    
    # Deactivate conda base environment if active, to avoid conflicts with pip
    if command_exists conda; then
        print_warning "Conda detected. Deactivating conda base environment to avoid conflicts..."
        eval "$(conda shell.bash hook)"
        conda deactivate 2>/dev/null || true
        print_success "Conda base environment deactivated"
    fi
    
    # Activate virtual environment
    echo "Activating virtual environment..."
    source "${VENV_DIR}/bin/activate"
    
    print_success "Virtual environment activated: ${VIRTUAL_ENV}"
    
    # Install PyTorch
    echo "Installing PyTorch and NumPy and dependencies (this may take a few minutes)..."
    pip3 install torch numpy
    
    print_success "PyTorch and NumPy installed successfully"
    
    # Verify PyTorch installation
    echo "Verifying PyTorch installation..."
    python3 -c "import torch; print(f'PyTorch version: {torch.__version__}')" || {
        print_error "PyTorch verification failed"
        exit 1
    }
    
    # Check MPS (Metal Performance Shaders) availability for Apple Silicon
    if [[ $(uname -m) == "arm64" ]]; then
        echo "Checking Metal Performance Shaders (GPU) support..."
        python3 -c "import torch; print(f'MPS available: {torch.backends.mps.is_available()}')"
    fi
    
    print_success "Verification: PyTorch is working correctly"
}

################################################################################
# STEP 8: Set and Verify Environment Variables
################################################################################

setup_environment_variables() {
    print_step 8 "Setting Up Environment Variables for libtorch"
    
    SETUP_SCRIPT="${INSTALL_DIR}/setup_modkit_env.sh"
    
    # Generate the environment setup script at the installation directory
    echo "Generating environment setup script at: ${SETUP_SCRIPT}"
    cat > "${SETUP_SCRIPT}" << 'SETUP_EOF'
#!/bin/bash
################################################################################
# Modkit Environment Setup Script
################################################################################
#
# This script sets up the environment variables required to run modkit
# compiled with PyTorch (libtorch) and macOS GPU (MPS) support.
#
# IMPORTANT: This script must be SOURCED, not executed, so that the
# environment variables persist in your current shell session.
#
# Usage:
#   source setup_modkit_env.sh <installation_directory>
#
# Arguments:
#   installation_directory: Required. The path used during modkit installation
#                           (the same path passed to mac_compile_modkit.sh).
#
# Examples:
#   source setup_modkit_env.sh /path/to/install
#   source setup_modkit_env.sh ~/tools
#
# What this script does:
#   1. Deactivates any active conda environment to avoid conflicts
#   2. Activates the Python virtual environment (venv_modkit)
#   3. Detects the Python version inside the virtual environment
#   4. Exports LIBTORCH_USE_PYTORCH, LIBTORCH, DYLD_LIBRARY_PATH, LD_LIBRARY_PATH
#   5. Adds the modkit binary directory to PATH
#   6. Verifies that all required paths exist
#
# After sourcing, you can run modkit directly:
#   modkit --version
#   modkit pileup input.bam output.bed
#   modkit open-chromatin predict --device mps -i input.bam -o output.bed
#
################################################################################

# Detect if script is being sourced or executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced, not executed."
    echo "Usage: source ${0} <installation_directory>"
    exit 1
fi

# Check for mandatory argument
if [[ $# -lt 1 ]]; then
    echo "Usage: source setup_modkit_env.sh <installation_directory>"
    echo "  installation_directory: Required. Path used during modkit installation."
    return 1
fi

MODKIT_INSTALL_DIR="$1"

################################################################################
# Validate installation directory
################################################################################

if [[ ! -d "${MODKIT_INSTALL_DIR}" ]]; then
    echo "Error: Installation directory not found: ${MODKIT_INSTALL_DIR}"
    return 1
fi

MODKIT_VENV_DIR="${MODKIT_INSTALL_DIR}/venv_modkit"
MODKIT_REPO_DIR="${MODKIT_INSTALL_DIR}/modkit"
MODKIT_BINARY="${MODKIT_REPO_DIR}/target/release/modkit"

if [[ ! -d "${MODKIT_VENV_DIR}" ]]; then
    echo "Error: Virtual environment not found: ${MODKIT_VENV_DIR}"
    echo "Has modkit been installed to ${MODKIT_INSTALL_DIR}?"
    return 1
fi

if [[ ! -f "${MODKIT_BINARY}" ]]; then
    echo "Warning: modkit binary not found at: ${MODKIT_BINARY}"
    echo "The environment will be set up, but modkit may not be runnable."
fi

################################################################################
# Step 1: Deactivate conda if active
################################################################################

if command -v conda &> /dev/null; then
    if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
        echo "Deactivating conda environment '${CONDA_DEFAULT_ENV}' to avoid conflicts..."
        eval "$(conda shell.bash hook)"
        conda deactivate 2>/dev/null || true
    fi
fi

################################################################################
# Step 2: Activate the Python virtual environment
################################################################################

echo "Activating virtual environment: ${MODKIT_VENV_DIR}"
source "${MODKIT_VENV_DIR}/bin/activate"

################################################################################
# Step 3: Detect Python version and set environment variables
################################################################################

MODKIT_PYTHON_VER=$(python3 -c "import sys; print(f'python{sys.version_info.major}.{sys.version_info.minor}')")

export LIBTORCH_USE_PYTORCH=1
export LIBTORCH_BYPASS_VERSION_CHECK=1
export LIBTORCH="${MODKIT_VENV_DIR}/lib/${MODKIT_PYTHON_VER}/site-packages/torch"
export DYLD_LIBRARY_PATH="${LIBTORCH}/lib"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${DYLD_LIBRARY_PATH}"

################################################################################
# Step 4: Add modkit binary to PATH
################################################################################

export PATH="${MODKIT_REPO_DIR}/target/release:${PATH}"

################################################################################
# Step 5: Verify setup
################################################################################

echo ""
echo "Modkit environment configured:"
echo "  LIBTORCH_USE_PYTORCH          = ${LIBTORCH_USE_PYTORCH}"
echo "  LIBTORCH_BYPASS_VERSION_CHECK = ${LIBTORCH_BYPASS_VERSION_CHECK}"
echo "  LIBTORCH                      = ${LIBTORCH}"
echo "  DYLD_LIBRARY_PATH             = ${DYLD_LIBRARY_PATH}"
echo "  Python                        = $(which python3) (${MODKIT_PYTHON_VER})"

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
    echo ""
    echo "Environment variables are set, but modkit binary is missing."
    echo "You may need to compile modkit first."
fi

# Clean up local variables to avoid polluting the caller's namespace
unset MODKIT_INSTALL_DIR MODKIT_VENV_DIR MODKIT_REPO_DIR MODKIT_BINARY MODKIT_PYTHON_VER
SETUP_EOF

    chmod +x "${SETUP_SCRIPT}"
    print_success "Environment setup script generated at: ${SETUP_SCRIPT}"
    
    # Source the generated script to configure the current session
    echo "Sourcing environment setup script..."
    # Save variables that the setup script's cleanup will unset
    local SAVED_MODKIT_REPO_DIR="${MODKIT_REPO_DIR}"
    local SAVED_VENV_DIR="${VENV_DIR}"
    
    source "${SETUP_SCRIPT}" "${INSTALL_DIR}"
    
    # Restore variables needed by subsequent build steps
    MODKIT_REPO_DIR="${SAVED_MODKIT_REPO_DIR}"
    VENV_DIR="${SAVED_VENV_DIR}"
    
    # Additional verification of paths
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

################################################################################
# STEP 9: Build Modkit with macOS GPU Support (continued)
################################################################################

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
        echo "  1. Ensure all environment variables are set correctly"
        echo "  2. Check that PyTorch is installed: pip3 list | grep torch"
        echo "  3. Try cleaning and rebuilding: cargo clean && cargo build --release --features accelerate,tch"
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

################################################################################
# Save Installation Info
################################################################################

save_installation_info() {
    INFO_FILE="${INSTALL_DIR}/installation_info.txt"
    SETUP_SCRIPT="${INSTALL_DIR}/setup_modkit_env.sh"
    
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
Python Version: $(python3 --version)
PyTorch Version: $(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "Unknown")

Quick Start:
-----------
source "${SETUP_SCRIPT}" "${INSTALL_DIR}"
modkit --version

For detailed usage instructions, run:
modkit --help
EOF
    
    echo "Installation information saved to: ${INFO_FILE}"
}

################################################################################
# Main Installation Flow
################################################################################

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
    echo ""
    echo "The installation process includes 9 steps and may take 15-30 minutes."
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
    create_venv
    install_pytorch
    setup_environment_variables
    build_modkit
    
    # Save installation info
    save_installation_info
    
    # Calculate installation time
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    MINUTES=$((DURATION / 60))
    REMAINING_SECONDS=$((DURATION % 60))
    
    echo ""
    echo -e "${GREEN}Total installation time: ${MINUTES} minutes ${REMAINING_SECONDS} seconds${NC}"
    echo ""
    
}

################################################################################
# Script Entry Point
################################################################################

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
if [[ "${MACOS_MAJOR}" -lt 11 ]]; then
    print_error "macOS 11 (Big Sur) or later is required for Metal GPU support"
    echo "Detected macOS version: ${MACOS_VERSION} (major: ${MACOS_MAJOR})"
    echo "Please upgrade your macOS before running this script."
    exit 1
fi

# Run main installation
main

exit 0
