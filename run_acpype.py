#!/usr/bin/env python3

import re
import sys

from acpype.cli import init_main

if __name__ == "__main__":
    sys.argv[0] = re.sub(r"(-script\.pyw?|\.exe)?$", "", sys.argv[0])
    sys.exit(init_main())
