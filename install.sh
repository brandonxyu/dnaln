#!/bin/sh
# Installs the latest prebuilt dnaln release (macOS, or Linux x86-64 / ARM64):
#
#   curl -fsSL https://raw.githubusercontent.com/brandonxyu/dnaln/main/install.sh | sh
#
# Puts the `dnaln` program in ~/.local/bin (override with DNALN_INSTALL_DIR=...).
set -eu

REPO="${DNALN_REPO:-brandonxyu/dnaln}"
INSTALL_DIR="${DNALN_INSTALL_DIR:-$HOME/.local/bin}"

case "$(uname -s)-$(uname -m)" in
  Darwin-*) asset=dnaln-macos-universal ;;
  Linux-x86_64) asset=dnaln-linux-x86_64 ;;
  Linux-aarch64 | Linux-arm64) asset=dnaln-linux-arm64 ;;
  *)
    echo "No prebuilt dnaln for $(uname -s) $(uname -m). Build from source instead:" >&2
    echo "  git clone https://github.com/$REPO.git && cd dnaln && make && make install" >&2
    exit 1
    ;;
esac

url="${DNALN_BASE_URL:-https://github.com/$REPO/releases/latest/download}/$asset.tar.gz"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "Downloading $url"
curl -fsSL "$url" | tar xz -C "$tmp"
mkdir -p "$INSTALL_DIR"
cp "$tmp/$asset/dnaln" "$INSTALL_DIR/dnaln"
chmod 755 "$INSTALL_DIR/dnaln"
echo "Installed $("$INSTALL_DIR/dnaln" --version) to $INSTALL_DIR/dnaln"

case ":$PATH:" in
  *":$INSTALL_DIR:"*) echo "Try it: dnaln --help" ;;
  *)
    profile="$HOME/.profile"
    case "${SHELL:-}" in */zsh) profile="$HOME/.zshrc" ;; */bash) profile="$HOME/.bashrc" ;; esac
    echo
    echo "$INSTALL_DIR is not on your PATH yet. Run this once, then open a new terminal:"
    echo "  echo 'export PATH=\"$INSTALL_DIR:\$PATH\"' >> $profile"
    ;;
esac
