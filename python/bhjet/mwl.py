"""Helpers for preparing multi-wavelength flux points for Gammapy."""

from pathlib import Path

from astropy.table import Table
import astropy.units as u
import astropy.constants as const
import astropy.units as u
import yaml


DEFAULT_FRACTIONAL_ENERGY_HALF_WIDTH = 0.001
"""Default half-width for a nominally monochromatic flux point (0.1 percent)."""


def radio_yaml_to_dnde_csv(
    yaml_path,
    output_path,
    fractional_energy_half_width=DEFAULT_FRACTIONAL_ENERGY_HALF_WIDTH,
    instrument_energy_half_widths=None,
):
    """Convert YAML flux densities to a Gammapy-style CSV file.

    Every top-level YAML list is treated as a flux-point band. Each entry needs
    ``frequency_Hz``, ``flux_mjy``, and ``flux_err``. Other fields such as
    ``instrument``, ``date``, ``citation``, and ``note`` are retained in the
    input YAML but are not written to this numerical CSV.

    Radio flux density is converted using ``E**2 dN/dE = nu F_nu``. The output
    columns are ``e_min``, ``e_ref``, ``e_max``, ``dnde``, and ``dnde_err`` in
    eV, eV, eV, cm-2 s-1 eV-1, and cm-2 s-1 eV-1 respectively.

    ``e_min`` and ``e_max`` are nominal bounds for a monochromatic measurement:
    ``e_ref * (1 - fractional_energy_half_width)`` and
    ``e_ref * (1 + fractional_energy_half_width)``. Set
    ``instrument_energy_half_widths``, for example ``{"VLBA": 0.005}``, to
    use a different fractional half-width for an instrument.

    Parameters
    ----------
    yaml_path : str or pathlib.Path
        Input YAML file containing one or more flux-point lists.
    output_path : str or pathlib.Path
        CSV file to create. Parent directories are created when needed.
    Returns
    -------
    pathlib.Path
        The written CSV file.
    """
    yaml_path = Path(yaml_path)
    output_path = Path(output_path)
    instrument_energy_half_widths = instrument_energy_half_widths or {}

    _validate_fractional_half_width(
        fractional_energy_half_width,
        "fractional_energy_half_width",
    )

    for instrument, half_width in instrument_energy_half_widths.items():
        _validate_fractional_half_width(
            half_width,
            f"instrument_energy_half_widths[{instrument!r}]",
        )
    with yaml_path.open() as file:
        data = yaml.safe_load(file)

    if not isinstance(data, dict):
        raise ValueError(f"{yaml_path} must contain top-level flux-point lists.")

    rows = []
    for band, points in data.items():
        if not isinstance(points, list):
            raise ValueError(f"Top-level entry {band!r} must be a list.")

        for index, point in enumerate(points):
            if not isinstance(point, dict):
                raise ValueError(f"{band} entry {index} must be a mapping.")

            frequency_hz = _required_float(point, "frequency_Hz", f"{band} entry {index}")
            flux_mjy = _required_float(point, "flux_mjy", f"{band} entry {index}")
            flux_err_mjy = _required_float(point, "flux_err", f"{band} entry {index}")

            if frequency_hz <= 0:
                raise ValueError(f"{band} entry {index} has non-positive frequency_Hz.")
            if flux_mjy < 0 or flux_err_mjy < 0:
                raise ValueError(
                    f"{band} entry {index} has a negative flux_mjy or flux_err."
                )

            e_ref = (frequency_hz * u.Hz).to_value(
                u.eV,
                equivalencies=u.spectral(),
            )
            half_width = instrument_energy_half_widths.get(
                point.get("instrument"),
                fractional_energy_half_width,
            )

            dnde = _mjy_to_dnde(e_ref, flux_mjy)
            dnde_err = _mjy_to_dnde(e_ref, flux_err_mjy)

            rows.append(
                (
                    e_ref * (1.0 - half_width),
                    e_ref,
                    e_ref * (1.0 + half_width),
                    dnde,
                    dnde_err,
                )
            )

    rows.sort(key=lambda row: row[1])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as file:
        file.write("e_min e_ref e_max dnde dnde_err\n")

        for row in rows:
            file.write(" ".join(f"{value:.16e}" for value in row) + "\n")

    return output_path


def _mjy_to_dnde(energy_eV, flux_mjy):
    """Convert a flux density in mJy to dN/dE in cm-2 s-1 eV-1."""
    flux_density_cgs = flux_mjy * 1e-26
    return flux_density_cgs / (const.h.cgs.value * energy_eV)


def _required_float(point, name, index):
    try:
        return float(point[name])
    except KeyError as error:
        raise ValueError(f"radio entry {index} is missing {name!r}.") from error
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"radio entry {index} has a non-numeric {name!r}."
        ) from error


def _validate_fractional_half_width(value, name):
    if not 0.0 < float(value) < 1.0:
        raise ValueError(f"{name} must be greater than 0 and less than 1.")


def read_dnde_table(path):
    table = Table.read(path, format="ascii.csv", delimiter=" ")
    table["e_min"].unit = u.eV
    table["e_ref"].unit = u.eV
    table["e_max"].unit = u.eV
    table["dnde"].unit = u.Unit("cm-2 s-1 eV-1")
    table["dnde_err"].unit = u.Unit("cm-2 s-1 eV-1")
    return table


def read_dnde_flux_points(path):
    """Read each narrow-bin CSV row as an independent Gammapy flux point.

    A single ``FluxPoints`` object requires a continuous energy axis. This is
    unsuitable for sparse multi-wavelength measurements because it turns the
    gaps between bands into artificial wide bins. Returning one object per row
    preserves the ``e_min`` and ``e_max`` values written for that measurement.
    """
    from gammapy.estimators import FluxPoints

    table = read_dnde_table(path)
    return [
        FluxPoints.from_table(table[index : index + 1], sed_type="dnde")
        for index in range(len(table))
    ]
