@echo on
docker run --rm -i -t -v %cd%:/wdir -w /wdir acpype/acpype acpype %*
