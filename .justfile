default: build test

configure:
  cmake --preset=dev

build:
  cmake --build --preset=dev

test:
  ./build/dev/tests/test_perturb
