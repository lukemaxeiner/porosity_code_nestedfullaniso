#!/usr/bin/env python3
"""Plot a histogram of porosity values from a legacy ASCII VTK file."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile


def read_vtk_porosity(vtk_path: Path) -> np.ndarray:
    """Read the point scalars named ``porosity`` from a legacy ASCII VTK file."""
    lines = vtk_path.read_text(encoding="utf-8").splitlines()

    scalar_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.strip().upper().startswith("SCALARS POROSITY")
        ),
        None,
    )
    if scalar_index is None:
        raise ValueError("No SCALARS porosity array was found in the VTK file.")

    lookup_index = next(
        (
            index
            for index in range(scalar_index + 1, len(lines))
            if lines[index].strip().upper().startswith("LOOKUP_TABLE")
        ),
        None,
    )
    if lookup_index is None:
        raise ValueError("No LOOKUP_TABLE entry was found after the porosity scalar array.")

    values: list[float] = []
    for line in lines[lookup_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        try:
            values.extend(float(token) for token in stripped.split())
        except ValueError as error:
            raise ValueError(
                f"Unexpected non-numeric data in the porosity array near: {stripped[:80]!r}"
            ) from error

    if not values:
        raise ValueError("The porosity scalar array contains no values.")

    return np.asarray(values, dtype=float)


def describe_values(values: np.ndarray) -> dict[str, float]:
    """Calculate summary statistics for one realization."""
    outside = (values < 0.0) | (values > 1.0)
    return {
        "count": float(values.size),
        "mean": float(np.mean(values)),
        "variance": float(np.var(values)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "below_zero": float(np.count_nonzero(values < 0.0)),
        "above_one": float(np.count_nonzero(values > 1.0)),
        "outside_fraction": float(np.count_nonzero(outside) / values.size),
    }


def parse_args() -> argparse.Namespace:
    project_folder = Path(__file__).resolve().parent.parent
    default_vtk = (
        project_folder
        / "GENERATED_POROSITY_FIELDS"
        / "S4_1.0487_2201_NI8crop_gaussian_realization_001.vtk"
    )

    parser = argparse.ArgumentParser(
        description="Plot a histogram of porosity values from a legacy ASCII VTK file."
    )
    parser.add_argument(
        "--vtk",
        type=Path,
        default=default_vtk,
        help="VTK file to visualize (default: existing realization 001)",
    )
    parser.add_argument("--bins", type=int, default=50, help="Histogram bin count")
    parser.add_argument(
        "--physical-range",
        action="store_true",
        help="Display only values from 0 to 1 (source values are unchanged)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "porosity_histogram.png",
        help="Where to save the histogram image",
    )
    parser.add_argument("--no-show", action="store_true", help="Do not open an interactive window")
    parser.add_argument(
        "--measured",
        type=Path,
        help="Measured blockwise porosity TIFF to overlay",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.bins < 1:
        raise SystemExit("--bins must be at least 1")

    vtk_path = args.vtk.expanduser().resolve()
    if not vtk_path.is_file():
        raise SystemExit(f"VTK file not found: {vtk_path}")

    values = read_vtk_porosity(vtk_path)
    statistics = describe_values(values)

    plot_values = values[np.isfinite(values)]
    if args.physical_range:
        plot_values = plot_values[(plot_values >= 0.0) & (plot_values <= 1.0)]

    measured_values = None
    measured_path = None
    measured_statistics = None
    if args.measured:
        measured_path = args.measured.expanduser().resolve()
        if not measured_path.is_file():
            raise SystemExit(f"Measured porosity TIFF not found: {measured_path}")

        measured_values = tifffile.imread(measured_path).astype(np.float64).ravel()
        measured_values = measured_values[np.isfinite(measured_values)]
        if measured_values.size == 0:
            raise SystemExit("The measured porosity TIFF contains no finite values.")
        measured_statistics = describe_values(measured_values)

        if args.physical_range:
            measured_values = measured_values[
                (measured_values >= 0.0) & (measured_values <= 1.0)
            ]
            if measured_values.size == 0:
                raise SystemExit("No measured porosity values lie in the range 0 to 1.")

    fig, axis = plt.subplots(figsize=(9, 6))

    if measured_values is not None:
        if args.physical_range:
            histogram_bins = np.linspace(0.0, 1.0, args.bins + 1)
        else:
            combined_values = np.concatenate((plot_values, measured_values))
            histogram_bins = np.histogram_bin_edges(combined_values, bins=args.bins)

        axis.hist(
            measured_values,
            bins=histogram_bins,
            density=True,
            alpha=0.35,
            edgecolor="black",
            color="tab:blue",
            label="Measured",
        )
        axis.hist(
            plot_values,
            bins=histogram_bins,
            density=True,
            histtype="step",
            linewidth=2.0,
            color="tab:orange",
            label="Reconstructed",
        )
        axis.set_ylabel("Probability density")
        axis.legend()
    else:
        axis.hist(plot_values, bins=args.bins, edgecolor="black", alpha=0.75)
        axis.set_ylabel("Voxel count")

    axis.set_xlabel("Porosity")
    axis.set_title(f"Porosity histogram: {vtk_path.name}")
    if args.physical_range:
        axis.set_xlim(0.0, 1.0)

    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()

    args.output.expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.expanduser().resolve(), dpi=150)

    print(f"VTK file: {vtk_path}")
    print(f"Values: {statistics['count']:.0f}")
    print(f"Mean: {statistics['mean']:.6f}")
    print(f"Variance: {statistics['variance']:.6f}")
    print(f"Minimum: {statistics['minimum']:.6f}")
    print(f"Maximum: {statistics['maximum']:.6f}")
    print(f"Below 0: {statistics['below_zero']:.0f}")
    print(f"Above 1: {statistics['above_one']:.0f}")
    print(f"Outside 0-1 fraction: {statistics['outside_fraction']:.6f}")

    if measured_path is not None:
        print(f"Measured TIFF: {measured_path}")
        print(f"Measured values: {measured_statistics['count']:.0f}")
        print(f"Measured mean: {measured_statistics['mean']:.6f}")
        print(f"Measured variance: {measured_statistics['variance']:.6f}")
        print(f"Measured minimum: {measured_statistics['minimum']:.6f}")
        print(f"Measured maximum: {measured_statistics['maximum']:.6f}")

    print(f"Saved histogram: {args.output.expanduser().resolve()}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
