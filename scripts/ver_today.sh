#!/usr/bin/env bash

# Stamp today's date as the project version, and point the README's
# "commits since release" badge at the newest tag that actually exists.
#
# This used to shell out to `pcregrep`, which is no longer packaged by Homebrew
# (it ships `pcre2grep` now), so the loops silently iterated over nothing and the
# version stopped being updated. It also globbed `./*/*.py`, which stopped matching
# the package once it moved to `src/acpype/`. Both failures were invisible because
# the script always exited 0.

set -euo pipefail

today="$(date +%Y.%-m.%-d)"
last_tag="$(git describe --abbrev=0 --tags)"

# Matches a YYYY.M.D version. The word boundaries matter: without them this would
# also rewrite unrelated numbers such as the `2005.12.00` citation in the README.
version_re='[0-9]{4}\.[12]?[0-9]\.[123]?[0-9]'

# The only files carrying the project version.
version_files=(src/acpype/__init__.py pyproject.toml)

for afile in "${version_files[@]}"; do
    if ! grep -qwE "$version_re" "$afile"; then
        echo "No version string found in $afile -- has it moved?" >&2
        exit 1
    fi
    echo "Updated version to $today in file $afile"
    perl -pi -e "s/\b${version_re}\b/$today/g" "$afile"
    git add "$afile"
done

# The README badges must reference a release that exists, not the pending version.
if grep -qwE "$version_re" README.md; then
    echo "Updated last_tag to $last_tag in file README.md"
    perl -pi -e "s/\b${version_re}\b/$last_tag/g" README.md
    git add README.md
fi

# The version stamped above lives in uv.lock too; CI runs `uv sync --locked`, so the
# lock has to be refreshed in the same commit or every build fails on a stale lock.
if command -v uv >/dev/null 2>&1; then
    uv lock --quiet
    git add uv.lock
fi

exit 0
