#!/bin/sh
# Installs the latest prebuilt dnaln release (macOS 11+, or Linux x86-64 / ARM64):
#
#   curl -fsSL https://raw.githubusercontent.com/brandonxyu/dnaln/main/install.sh | sh
#
# Puts the `dnaln` program in ~/.local/bin (override with DNALN_INSTALL_DIR=...).
set -eu

REPO="${DNALN_REPO:-brandonxyu/dnaln}"
INSTALL_DIR="${DNALN_INSTALL_DIR:-$HOME/.local/bin}"
BASE_URL="${DNALN_BASE_URL:-https://github.com/$REPO/releases/latest/download}"

fail() {
  echo "error: $*" >&2
  exit 1
}

case "$(uname -s)-$(uname -m)" in
  Darwin-*) asset=dnaln-macos-universal ;;
  Linux-x86_64) asset=dnaln-linux-x86_64 ;;
  Linux-aarch64 | Linux-arm64) asset=dnaln-linux-arm64 ;;
  *)
    fail "no prebuilt dnaln for $(uname -s) $(uname -m). Build from source instead:
  git clone https://github.com/$REPO.git && cd dnaln && make && make install"
    ;;
esac

if command -v sha256sum >/dev/null 2>&1; then
  sha256() { sha256sum "$1" | cut -d' ' -f1; }
elif command -v shasum >/dev/null 2>&1; then
  sha256() { shasum -a 256 "$1" | cut -d' ' -f1; }
else
  fail "need sha256sum or shasum to verify the download"
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "Downloading $BASE_URL/$asset.tar.gz"
curl -fsSL -o "$tmp/$asset.tar.gz" "$BASE_URL/$asset.tar.gz" ||
  fail "download failed. Check your connection and that a release exists at https://github.com/$REPO/releases"
curl -fsSL -o "$tmp/SHA256SUMS" "$BASE_URL/SHA256SUMS" || fail "could not download SHA256SUMS"

expected="$(grep " $asset.tar.gz\$" "$tmp/SHA256SUMS" | cut -d' ' -f1)"
[ -n "$expected" ] || fail "$asset.tar.gz is not listed in SHA256SUMS"
[ "$(sha256 "$tmp/$asset.tar.gz")" = "$expected" ] || fail "checksum mismatch for $asset.tar.gz; not installing"

tar xzf "$tmp/$asset.tar.gz" -C "$tmp"
"$tmp/$asset/dnaln" --version >/dev/null 2>&1 ||
  fail "the downloaded dnaln does not run on this system (Macs need macOS 11 or later)"
mkdir -p "$INSTALL_DIR"
cp "$tmp/$asset/dnaln" "$INSTALL_DIR/dnaln"
chmod 755 "$INSTALL_DIR/dnaln"
echo "Installed $("$INSTALL_DIR/dnaln" --version) to $INSTALL_DIR/dnaln (checksum verified)"

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
