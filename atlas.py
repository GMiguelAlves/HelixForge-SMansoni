#!/usr/bin/env python3
"""Start the frozen atlas without installing analysis dependencies."""
from scripts.atlas.runtime import main

if __name__ == "__main__":
    raise SystemExit(main())
