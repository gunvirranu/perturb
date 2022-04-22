
Personally, I'm not a huge fan of auto-generated documentation; I usually just prefer reading the comments in the header. Since that's just me, this small site contains some library documentation that may be useful.

I know this main page is pretty sparse. There's no point in copying and pasting _all_ the _very cool_ features of `perturb` from the repo's README when you could just look there.

**Go to [the `perturb` repository's main page here](https://github.com/gunvirranu/perturb) for a bunch of useful info on features, installation, and examples.**

## Basic Usage

Everything in this library is in the <tt>@ref perturb</tt> namespace. Here's a quick intro to the typical "usage flow".

I won't cover the details of [SGP4][SGP4], but in brief, it's a very popular orbit propagator for Earth-centered spacecraft. Usually, the input orbit ephemeris is through a [TLE][TLE], such as the ones provided by [Celestrek][Celestrek]. These TLE inputs can be used to construct a `perturb::Satellite` object.

A specific point in time is represented as a `perturb::JulianDate`. You can either construct one from a specific date and time via `perturb::DateTime` or offset a number of days from the `perturb::Satellite::epoch()` of a satellite.

Passing in a `perturb::JulianDate` time point to the `perturb::Satellite::propagate` method of a satellite yields a `perturb::StateVector`, which contains a time-stamp, and a position and velocity `perturb::Vec3`. The vectors are just a `std::array<double, 3>`, measured in kilometres, and are represented in the [TEME][ECI-TEME] coordinate reference frame. The details of this frame can get a bit annoying, so this library does _not_ handle converting it to others. For handling Earth-centered reference frames such as TEME and transformations between them, you may be interested in the [`gelocus`][gelocus] library.

## Build Options

See the ["Build Options" section of the repository page][build-options] for more info.

## License

The `perturb` library is licensed under the [MIT license][mit-license].

<!-- Links -->
[SGP4]: https://en.wikipedia.org/wiki/Simplified_perturbations_models
[TLE]: https://en.wikipedia.org/wiki/Two-line_element_set
[Celestrek]: https://celestrak.com
[ECI-TEME]: https://en.wikipedia.org/wiki/Earth-centered_inertial
[gelocus]: https://github.com/gunvirranu/gelocus
[build-options]: https://github.com/gunvirranu/perturb#build-options
[mit-license]: https://github.com/gunvirranu/perturb#license
