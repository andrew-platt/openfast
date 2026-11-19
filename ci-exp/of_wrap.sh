#!/bin/bash
# Run the real OpenFAST executable, optionally with a preloaded library.
if [ -n "$OF_PRELOAD" ]; then export LD_PRELOAD=$OF_PRELOAD; fi
exec "$OF_REAL" "$@"
