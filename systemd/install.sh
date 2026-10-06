#!/usr/bin/env bash

set -e

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNIT_DIR="${HOME}/.config/systemd/user"

mkdir -p "$UNIT_DIR"

echo "Installing Weather Wall systemd units..."

cp "$SOURCE_DIR"/*.service "$UNIT_DIR/"
cp "$SOURCE_DIR"/*.timer "$UNIT_DIR/"

systemctl --user daemon-reload

echo
echo "Installed Weather Wall units:"
systemctl --user list-unit-files 'ww-*'