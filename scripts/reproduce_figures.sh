#!/usr/bin/env bash
set +e

cd "$(dirname "$0")/.."

echo "========================================"
echo "REPRODUCE PUBLIC FIGURES"
echo "========================================"

python \
  src/make_publication_figures.py

RC=$?

echo
echo "Return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "Generated:"

    ls -lh \
      figures/

fi

echo
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
