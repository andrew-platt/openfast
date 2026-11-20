#!/usr/bin/env bash
# Enumerate every ctest test selected by each rtest-* job of automated-dev-tests.yml
# (filters copied verbatim from the workflow and .github/actions/tests-module-*).
# Output: TSV "<ci-job>\t<test-name>" on stdout. Run from the build directory.
set -euo pipefail
list() { # $1 = job tag, rest = ctest filter args
  local job=$1; shift
  ctest -N "$@" | sed -n 's/^ *Test *#[0-9]*: \(.*\)$/\1/p' | sed "s/^/${job}\t/"
}
list rtest-module-drivers  -R ad_      -LE python
list rtest-module-drivers  -R bd_
list rtest-module-drivers  -R hd_      -LE python
list rtest-module-drivers  -R ifw_     -LE python
list rtest-module-drivers  -R md_      -LE python
list rtest-module-drivers  -R seastate_ -LE python
list rtest-module-drivers  -R SD_
list rtest-interfaces      -L "cpp|python|fastlib" -LE "openfast_io"
list rtest-OF              -L openfast -LE "cpp|linear|python|fastlib|offshore"
list rtest-OF-offshore     -L openfast -L offshore -LE "cpp|linear|python|fastlib|multirotor"
list rtest-OF-multirotor   -L openfast -L multirotor -LE "cpp|linear|python|fastlib"
list rtest-OF-linearization -L linear
list rtest-OF-aeromap      -L aeromap -LE "cpp|linear|python"
list rtest-openfast_io     -L openfast_io
list rtest-FF              -L fastfarm
