#!/usr/bin/env bash
# Run one Chaa test case (or "all") and validate the result.
#   tests/run_case.sh <case-name>|all
# Env: CHAA_BIN (default bin/chaa), PY (default python3),
#      TEST_OUT (default test-output)
set -euo pipefail

HERE="$(cd "$(dirname "$0")/.." && pwd)"
CONF="$HERE/tests/cases.conf"
if [ -z "${CHAA_BIN:-}" ]; then
  if [ -x "$HERE/bin/chaa" ]; then CHAA_BIN="$HERE/bin/chaa"
  else CHAA_BIN="$HERE/build/bin/chaa"; fi
fi
: "${PY:=python3}"
: "${TEST_OUT:=$HERE/test-output}"

run_one() {
  local case="$1"
  local flags
  flags=$(awk -F'|' -v c="$case" \
    '$0 !~ /^[[:space:]]*(#|$)/ { name=$1; gsub(/[[:space:]]/, "", name);
      if (name == c) print $2 }' "$CONF")
  if [ -z "$flags" ]; then
    echo "unknown case: $case" >&2
    exit 2
  fi
  local out="$TEST_OUT/$case"
  rm -rf "$out"
  mkdir -p "$out"
  echo "=== running $case ==="
  # run from the repo root so relative --paramsFile paths resolve
  cd "$HERE"
  if [[ "$case" == restart-* ]]; then
    # graceful-stop/restart cases: (1) run with a pre-planted stop file
    # so the run saves a restart after step 1 and exits, (2) resume it,
    # (3) run uninterrupted into <case>-ref; the validator requires
    # machine-identical dumps.
    touch "$out/stop"
    # shellcheck disable=SC2086
    "$CHAA_BIN" $flags --outDir="$out" --logEvery=1000000000
    # shellcheck disable=SC2086
    "$CHAA_BIN" $flags --outDir="$out" --restart=true --logEvery=1000000000
    rm -rf "$out-ref"
    mkdir -p "$out-ref"
    # shellcheck disable=SC2086
    "$CHAA_BIN" $flags --outDir="$out-ref" --logEvery=1000000000
  else
    # shellcheck disable=SC2086
    "$CHAA_BIN" $flags --outDir="$out" --logEvery=1000000000
  fi
  "$PY" "$HERE/tests/validate/validate.py" "$case" "$out"
}

if [ "${1:-}" = "all" ]; then
  fails=0
  while read -r case; do
    if ! run_one "$case"; then
      echo "*** case $case FAILED"
      fails=$((fails+1))
    fi
  done < <(awk -F'|' '$0 !~ /^[[:space:]]*(#|$)/ { name=$1;
             gsub(/[[:space:]]/, "", name); print name }' "$CONF")
  echo "==============================="
  if [ "$fails" -gt 0 ]; then
    echo "$fails case(s) failed"
    exit 1
  fi
  echo "all cases passed"
else
  run_one "${1:?usage: run_case.sh <case>|all}"
fi
