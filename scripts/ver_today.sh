#!/bin/bash

today="$(date +%Y.%-m.%-d)"

while IFS= read -r afile; do
	echo "Updated version to $today in file $afile"
	perl -pi -e "s/\b[0-9]{4}\.[12]?[0-9]\.[12]?[0-9]\b/$today/g" "$afile"
	# git add "$afile"
done < <(grep -lEr "\b\d{4}\.[12]?[0-9]\.[12]?[0-9]\b" -- **/*.py README.md pyproject.toml)

exit 0
