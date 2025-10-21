#!/bin/bash
# setup_k4env.sh
# Start a completely fresh shell with the desired Key4hep environment

echo "Sourcing Key4hep release 2025-05-29..."
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29

# Source local project setup if it exists
if [ -f setup.sh ]; then
  echo "Sourcing local project setup.sh..."
  source setup.sh
else
  echo "Warning: setup.sh not found in current directory!"
fi

echo "Fresh Key4hep environment ready!"
