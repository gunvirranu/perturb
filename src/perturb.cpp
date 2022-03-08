#include "perturb/perturb.hpp"

using namespace perturb;

Foo::Foo() : num(73) {
}

int Foo::gimme_number() const {
  return num;
}
