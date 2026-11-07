# SPEIbase changelog

## SPEIbase v2.13.0

### Bug fixes

* PET conversion from mm/day to mm/month is now applied along the time
  dimension (`sweep(pet, 3, ndays, '*')`). Since v2.6, `pet * ndays`
  recycled the days-per-month vector along longitude instead of time, so
  each grid cell got one fixed multiplier (28, 29, 30 or 31) in every month.
  All SPEI datasets from v2.6 to v2.12.0 are affected; SPI is not. A
  corrected dataset release will follow (#19, reported by @arnabkdey).
* `spei.nc()` now stops with an error if the precipitation and PET files
  do not share the same time axis over the period being processed, because
  the days-per-month vector comes from the precipitation file and is applied
  to the PET data.

### Documentation

* README: added a note that the v2.5/v2.5.1 code (Zenodo release) does not
  convert PET from mm/day to mm/month, and that v2.6 to v2.12.0 applied the
  conversion incorrectly.

For earlier versions, see the version history in `README.md`.
