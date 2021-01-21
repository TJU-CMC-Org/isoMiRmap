#!/usr/bin/env python
import sys

if sys.version_info[0] < 3 or (sys.version_info[0] >= 3 and sys.version_info[1] < 5):
    raise RuntimeError('Must be using Python 3.5 or later')

import SourceLibrary.IsoMiRmap_Library as START
if __name__ == "__main__":
    START.start_generation()
