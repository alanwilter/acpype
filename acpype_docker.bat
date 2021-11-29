@echo on
docker run -i -t -v %cd%:/wdir -w /wdir acpype/acpype acpype %*