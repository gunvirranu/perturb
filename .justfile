# This file is for the just command runner
# github.com/casey/just
#
# Nothing special, just a handy way to remember commands

default: build test

configure:
  cmake --preset=dev

build:
  cmake --build --preset=dev

test:
  ./build/dev/tests/test_perturb
