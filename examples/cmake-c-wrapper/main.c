#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#include <perturb/perturb.h>

int main() {
    char * ISS_TLE_1 = "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996";
    char * ISS_TLE_2 = "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227";

    struct perturb_satellite sat;
    sat.p_satellite = malloc(PERTURB_C_SATELLITE_ALLOC_SIZE);

    const bool err = perturb_parse_tle_and_init_sat(ISS_TLE_1, ISS_TLE_2, sat);
    assert(err == false);

    struct perturb_state_vector sv = perturb_propagate_days(sat, 2.5);

    printf("Position [km]: { %f, %f, %f }\n", sv.position[0], sv.position[1], sv.position[2]);
    printf("Velocity [km/s]: { %f, %f, %f }\n", sv.velocity[0], sv.velocity[1], sv.velocity[2]);

    free(sat.p_satellite);

    /*
    // Let try simulating the orbit of the International Space Station
    // Got TLE from Celestrak sometime around 2022-03-12
    std::string ISS_TLE_1 = "1 25544U 98067A   22071.78032407  .00021395  00000-0  39008-3 0  9996";
    std::string ISS_TLE_2 = "2 25544  51.6424  94.0370 0004047 256.5103  89.8846 15.49386383330227";

    // Create and initialize a satellite object from the TLE
    auto sat = Satellite::from_tle(ISS_TLE_1, ISS_TLE_2);
    assert(sat.last_error() == Sgp4Error::NONE);
    assert(sat.epoch().to_datetime().day == 12);

    // Let's see what the ISS is doing on Pi Day
    const auto t = JulianDate(DateTime { 2022, 3, 14, 1, 59, 26.535 });
    const double delta_days = t - sat.epoch();
    assert(1 < delta_days && delta_days < 3);  // It's been ~2 days since the epoch

    // Calculate the position and velocity at the chosen time
    StateVector sv;
    const auto err = sat.propagate(t, sv);
    assert(err == Sgp4Error::NONE);
    const auto &pos = sv.position, &vel = sv.velocity;

    // Conclusion: The ISS is going pretty fast (~8 km/s)
    std::cout << "Position [km]: { " << pos[0] << ", " << pos[1] << ", " << pos[2] << " }\n";
    std::cout << "Velocity [km/s]: { " << vel[0] << ", " << vel[1] << ", " << vel[2] << " }\n";
    */
}
