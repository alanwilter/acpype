#!/bin/bash

today="$(date +%Y.%-m.%-d)"
last_tag=$(git describe --abbrev=0 --tags)

while IFS= read -r afile; do
	echo "Updated version to $today in file $afile"
	perl -pi -e "s/\b[0-9]{4}\.[12]?[0-9]\.[123]?[0-9]\b/$today/g" "$afile"
	git add "$afile"
done < <(pcregrep -lwr "\d{4}\.[12]?[0-9]\.[123]?[0-9]" ./*/*.py pyproject.toml)

while IFS= read -r afile; do
	echo "Updated last_tag to $last_tag in file $afile"
	perl -pi -e "s/\b[0-9]{4}\.[12]?[0-9]\.[123]?[0-9]\b/$last_tag/g" "$afile"
	git add "$afile"
done < <(pcregrep -lwr "\d{4}\.[12]?[0-9]\.[123]?[0-9]" README.md)

# The version stamped above lives in uv.lock too; CI runs `uv sync --locked`, so the
# lock has to be refreshed in the same commit or every build fails on a stale lock.
if command -v uv >/dev/null 2>&1; then
	uv lock --quiet
	git add uv.lock
fi

exit 0
