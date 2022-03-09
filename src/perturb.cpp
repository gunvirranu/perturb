#include "perturb/perturb.hpp"

namespace perturb {

Sgp4Error convert_sgp4_error_code(const int error_code) {
    if (error_code < 0 || error_code >= static_cast<int>(Sgp4Error::UNKNOWN)) {
        return Sgp4Error::UNKNOWN;
    }
    return static_cast<Sgp4Error>(error_code);
}

}  // namespace perturb
