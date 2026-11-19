#!/bin/bash
# usage: run_all.sh <mode> <rep>   (run from repo root; OF_REAL / OF_PRELOAD already exported)
mode=$1; rep=$2
export mode rep
export bd=$PWD/work/$mode
mkdir -p $bd results/$mode
CASES=$(grep -o 'of_regression_linear("[^"]*' reg_tests/CTestList.cmake | sed 's/.*"//' | awk '!s[$0]++')
run_one() {
  c=$1
  flag=$(grep -o "of_regression_linear(\"$c\" *\"[^\"]*\"" reg_tests/CTestList.cmake | head -1 | sed 's/"$//; s/.*"//')
  t0=$(date +%s)
  python reg_tests/executeOpenfastLinearRegressionCase.py $c $PWD/ci-exp/of_wrap.sh $PWD $bd 2 1.9 $flag > results/$mode/$c.log 2>&1
  rc=$?
  mkdir -p results/$mode/$c
  cp $bd/$c/*.lin results/$mode/$c/ 2>/dev/null
  h=$(cat results/$mode/$c/*.lin 2>/dev/null | grep -v -i -E 'predictions were generated|linked with' | sha256sum | cut -c1-10)
  echo "RESULT rep=$rep mode=$mode case=$c rc=$rc secs=$(( $(date +%s) - t0 )) linhash=$h"
}
export -f run_one
first=$(echo "$CASES" | head -1)
run_one $first
echo "$CASES" | tail -n +2 | xargs -P 4 -I{} bash -c 'run_one {}'
