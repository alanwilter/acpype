#!/bin/bash

today="$(date +%Y.%-m.%-d)"
last_tag=$(git describe --abbrev=0 --tags)

while IFS= read -r afile; do
	echo "Updated version to $today in file $afile"
	perl -pi -e "s/\b[0-9]{4}\.[12]?[0-9]\.[12]?[0-9]\b/$today/g" "$afile"
	git add "$afile"
done < <(grep -lEr "\b\d{4}\.[12]?[0-9]\.[12]?[0-9]\b" -- **/*.py pyproject.toml)

while IFS= read -r afile; do
	echo "Updated last_tag to $last_tag in file $afile"
	perl -pi -e "s/\b[0-9]{4}\.[12]?[0-9]\.[12]?[0-9]\b/$last_tag/g" "$afile"
	git add "$afile"
done < <(grep -lEr "\b\d{4}\.[12]?[0-9]\.[12]?[0-9]\b" -- README.md)

exit 0
