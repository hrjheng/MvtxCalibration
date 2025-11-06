#!/usr/bin/env python3
import os
import argparse
import sys

sys.path.insert(0, "./packages")

import numpy as np
import matplotlib.pyplot as plt

# global factor to translate threshold values from ADC to e- (electrons)
ADC_TO_ELECTRONS = 10.0  # This is a placeholder value, adjust as needed

stave_map = {
    "0": ["L0_03", "L0_04", "L2_01", "L2_02", "L0_05", "L2_03", "L2_04", "L2_15"],
    "1": ["L1_00", "L1_01", "L2_06", "L2_07", "L1_02", "L1_03", "L2_08", "L2_09"],
    "2": ["L0_00", "L0_01", "L1_04", "L2_10", "L0_02", "L1_05", "L1_06", "L1_07"],
    "3": ["L0_09", "L0_10", "L2_05", "L2_11", "L0_11", "L2_12", "L2_13", "L2_14"],
    "4": ["L0_06", "L0_07", "L1_12", "L2_00", "L0_08", "L1_13", "L1_14", "L1_15"],
    "5": ["L1_08", "L1_09", "L2_16", "L2_17", "L1_10", "L1_11", "L2_18", "L2_19"],
}


def read_stave_data(staveid, data_directory):
    feeids = sorted(
        [(int(staveid[1]) << 12) | int(staveid[3:]) | (gbt << 8) for gbt in range(3)]
    )
    data_thrs = np.zeros((3 * len(feeids), 512, 1024), dtype=np.float32)
    data_rmss = np.zeros((3 * len(feeids), 512, 1024), dtype=np.float32)
    for feeid in feeids:
        gbt_channel = feeids.index(feeid)
        file_thrs = "ls " + data_directory + f"/*FEEID{feeid}_thrs.dat"
        # file_rmss = "ls " + data_directory + f"/*FEEID{feeid}_rmss.dat"
        file_rmss = "ls " + data_directory + f"/*FEEID{feeid}_tnoise.dat"
        # get filename
        thrmap_file = os.popen(file_thrs).read().strip().split("\n")[0]
        rmsmap_file = os.popen(file_rmss).read().strip().split("\n")[0]
        if not thrmap_file or not rmsmap_file:
            print(
                f"Warning: either threshold or rms file could be found for stave {staveid} (feeid {feeid}). Skipping."
            )
            continue
        # read the data from the threshold file
        fin = open(thrmap_file, "r")
        thrdata = np.fromfile(fin, dtype=np.float32)
        for j in range(512):
            for k in range(1024 * 3):
                i = k // 1024
                data_thrs[gbt_channel * 3 + i][j][k - i * 1024] = thrdata[
                    j * 3 * 1024 + k
                ]
        fin.close()
        # read the data from the rms file
        fin = open(rmsmap_file, "r")
        rmsdata = np.fromfile(fin, dtype=np.float32)
        for j in range(512):
            for k in range(1024 * 3):
                i = k // 1024
                data_rmss[gbt_channel * 3 + i][j][k - i * 1024] = rmsdata[
                    j * 3 * 1024 + k
                ]
        fin.close()

    return data_thrs, data_rmss


def get_stave_ids(felix_number):
    """Get stave IDs for a given Felix number."""
    if felix_number < 0 or felix_number > 5:
        raise ValueError("Felix number must be between 0 and 5.")
    return stave_map[str(felix_number)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Threshold analysis.")
    parser.add_argument(
        "-f",
        "--felix",
        help="Felix number (0-5)",
        type=int,
        choices=range(6),
        required=True,
    )
    parser.add_argument("-d", help="Data directory", required=True)
    parser.add_argument("-o", help="Output file prefix", required=False, default="./")
    parser.add_argument(
        "-nsteps", help="Threshold step", type=int, default=50, required=False
    )
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.makedirs(args.o)

    staveids = get_stave_ids(args.felix)

    # The threshold map
    thrmap_fig, thrmap_ax = plt.subplots(
        len(staveids),
        9,
        figsize=(20, 1.5 * len(staveids)),
        sharex=True,
        sharey="row",
        gridspec_kw={"hspace": 0, "wspace": 0},
    )
    # The rms map
    rmsmap_fig, rmsmap_ax = plt.subplots(
        len(staveids),
        9,
        figsize=(20, 1.5 * len(staveids)),
        sharex=True,
        sharey="row",
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    staveid_idx = 0
    staveid_idx_max = len(staveids)
    print(f"Processing {staveid_idx_max} staves: {staveids}")
    if staveid_idx_max == 1:
        ax_for_cb = thrmap_ax[0]
        ax_for_cb_rms = rmsmap_ax[0]
    else:
        ax_for_cb = thrmap_ax[0][0]
        ax_for_cb_rms = rmsmap_ax[0][0]

    all_avg = {}  # Dictionary to store average thresholds for each stave
    all_std = {}  # Dictionary to store standard deviation of thresholds for each stave
    all_rms = {}  # Dictionary to store rms values for each stave
    all_rms_std = (
        {}
    )  # Dictionary to store standard deviation of rms values for each stave
    for stave in staveids:
        stave_data, stave_data_rms = read_stave_data(stave, args.d)
        print(
            f"Processing stave {stave} with threshold data shape {stave_data.shape} and rms data shape {stave_data_rms.shape}."
        )
        if (
            stave_data.ndim != 3
            or stave_data.shape[0] != 9
            or stave_data.shape[1] != 512
            or stave_data.shape[2] != 1024
            or stave_data_rms.ndim != 3
            or stave_data_rms.shape[0] != 9
            or stave_data_rms.shape[1] != 512
            or stave_data_rms.shape[2] != 1024
        ):
            print(
                f"Error: data for stave {stave} has an unexpected shape {stave_data.shape} for thresholds and/or {stave_data_rms.shape} for rms. Expected (9, 512, 1024)."
            )
            exit(1)

        # initialize arrays for average and standard deviation
        avg_chip = np.zeros(
            (9,), dtype=float
        )  # Array to store average thresholds for each chip
        std_chip = np.zeros(
            (9,), dtype=float
        )  # Array to store standard deviation of thresholds for each chip
        avgrms_chip = np.zeros(
            (9,), dtype=float
        )  # Array to store average rms values for each chip
        stdrms_chip = np.zeros(
            (9,), dtype=float
        )  # Array to store standard deviation of rms values for each chip
        nparr_avg_stave = np.zeros(
            (1,), dtype=float
        )  # Array to store all non-nan entries for the stave
        nparr_rms_stave = np.zeros(
            (1,), dtype=float
        )  # Array to store all non-nan rms entries for the stave

        # loop over chips and calculate average and standard deviation of non-zero and non-nan entries. each chip has 512 rows and 1024 columns so a total of 512*1024 entries
        for chip in range(0, stave_data.shape[0]):
            # get non-zero and non-nan entries for each chip
            non_zero_entries = stave_data[chip][stave_data[chip] != 0]
            non_nan_entries = non_zero_entries[~np.isnan(non_zero_entries)]
            print(
                f"Chip {chip} has {non_nan_entries.size} non-zero and non-nan entries."
            )
            avg_chip[chip] = np.mean(non_nan_entries) * ADC_TO_ELECTRONS
            std_chip[chip] = np.std(non_nan_entries) * ADC_TO_ELECTRONS
            print(
                f"Chip {chip} average: {avg_chip[chip]:.2f}, std: {std_chip[chip]:.2f}"
            )
            # append the non_nan_entries to avg_stave
            nparr_avg_stave = np.append(nparr_avg_stave, non_nan_entries)
            # get non-zero and non-nan entries for rms
            non_zero_entries_rms = stave_data_rms[chip][stave_data_rms[chip] != 0]
            non_nan_entries_rms = non_zero_entries_rms[~np.isnan(non_zero_entries_rms)]
            print(
                f"Chip {chip} has {non_nan_entries_rms.size} non-zero and non-nan rms entries."
            )
            avgrms_chip[chip] = np.mean(non_nan_entries_rms) * ADC_TO_ELECTRONS
            stdrms_chip[chip] = np.std(non_nan_entries_rms) * ADC_TO_ELECTRONS
            print(
                f"Chip {chip} rms average: {avgrms_chip[chip]:.2f}, std: {stdrms_chip[chip]:.2f}"
            )
            # append the non_nan_entries to rms_stave
            nparr_rms_stave = np.append(nparr_rms_stave, non_nan_entries_rms)

        avg_stave = np.mean(nparr_avg_stave) * ADC_TO_ELECTRONS
        std_stave = np.std(nparr_avg_stave) * ADC_TO_ELECTRONS
        rms_stave = np.mean(nparr_rms_stave) * ADC_TO_ELECTRONS
        rms_std_stave = np.std(nparr_rms_stave) * ADC_TO_ELECTRONS

        for chip in range(0, stave_data.shape[0]):
            if staveid_idx_max == 1:
                ax_thrmap = thrmap_ax[chip]
                ax_rmsmap = rmsmap_ax[chip]
            else:
                ax_thrmap = thrmap_ax[staveid_idx][chip]
                ax_rmsmap = rmsmap_ax[staveid_idx][chip]

            if chip == 0:
                ax_thrmap.set_ylabel(f"{stave}", rotation=90, labelpad=20, fontsize=14)
                ax_rmsmap.set_ylabel(f"{stave}", rotation=90, labelpad=20, fontsize=14)
            else:
                ax_thrmap.yaxis.set_visible(False)
                ax_rmsmap.yaxis.set_visible(False)

            if staveid_idx == staveid_idx_max - 1:
                ax_thrmap.set_xlabel(
                    f"Chip {chip}", rotation=0, labelpad=20, fontsize=14
                )
                ax_thrmap.xaxis.set_ticks([500, 1000])
                ax_rmsmap.set_xlabel(
                    f"Chip {chip}", rotation=0, labelpad=20, fontsize=14
                )
                ax_rmsmap.xaxis.set_ticks([500, 1000])
            else:
                ax_thrmap.xaxis.set_visible(False)
                ax_rmsmap.xaxis.set_visible(False)

            chip_thrmap = stave_data[chip] * ADC_TO_ELECTRONS
            ax_thrmap.invert_yaxis()
            cmap = plt.colormaps["Reds"].copy()
            ax_thrmap.imshow(chip_thrmap, cmap=cmap, aspect="auto", origin="lower")
            ax_thrmap.text(
                0.5,
                0.5,
                rf"{avg_chip[chip]:.2f} $\pm$ {std_chip[chip]:.2f}",
                horizontalalignment="center",
                transform=ax_thrmap.transAxes,
                fontsize=10,
            )

            chip_rmsmap = stave_data_rms[chip] * ADC_TO_ELECTRONS
            ax_rmsmap.invert_yaxis()
            cmap_rms = plt.colormaps["Reds"].copy()
            ax_rmsmap.imshow(chip_rmsmap, cmap=cmap_rms, aspect="auto", origin="lower")
            ax_rmsmap.text(
                0.5,
                0.5,
                rf"{avgrms_chip[chip]:.2f} $\pm$ {stdrms_chip[chip]:.2f}",
                horizontalalignment="center",
                transform=ax_rmsmap.transAxes,
                fontsize=10,
            )

        all_avg[staveid_idx] = avg_stave
        all_std[staveid_idx] = std_stave
        print(f"Stave {stave} average: {avg_stave:.2f}, std: {std_stave:.2f}")
        all_rms[staveid_idx] = rms_stave
        all_rms_std[staveid_idx] = rms_std_stave
        staveid_idx += 1

    thrmap_fig.colorbar(
        ax_for_cb.images[0],
        format="%.1f",
        ax=thrmap_ax,
        label="Threshold (e$^-$)",
        orientation="vertical",
    )
    thrmap_fig.savefig(f"{args.o}/mvtx{args.felix}.png", dpi=300, bbox_inches="tight")
    thrmap_fig.savefig(f"{args.o}/mvtx{args.felix}.pdf", bbox_inches="tight")
    plt.close(thrmap_fig)

    rmsmap_fig.colorbar(
        ax_for_cb_rms.images[0],
        format="%.1f",
        ax=rmsmap_ax,
        label="RMS (e$^-$)",
        orientation="vertical",
    )
    rmsmap_fig.savefig(
        f"{args.o}/mvtx{args.felix}_rms.png", dpi=300, bbox_inches="tight"
    )
    rmsmap_fig.savefig(f"{args.o}/mvtx{args.felix}_rms.pdf", bbox_inches="tight")
    plt.close(rmsmap_fig)

    # print the map
    print("Threshold map (avg):", all_avg)
    print("Threshold map (std):", all_std)
    print("RMS map (avg):", all_rms)
    print("RMS map (std):", all_rms_std)

    thrs_file = f"{args.o}/mvtx{args.felix}_thresholds.txt"
    with open(thrs_file, "w") as f:
        for stave in staveids:
            stave_idx = staveids.index(stave)
            mean_thr = all_avg[stave_idx]
            std_thr = all_std[stave_idx]
            f.write(f"{stave}: {mean_thr:.2f} e$^-$ $\\pm$ {std_thr:.2f} e$^-$\n")

    print(f"Thresholds saved to {thrs_file}")

    # rms_file = f"{args.o}/mvtx{args.felix}_rms.txt"
    rms_file = f"{args.o}/mvtx{args.felix}_tnoise.txt"
    with open(rms_file, "w") as f:
        for stave in staveids:
            stave_idx = staveids.index(stave)
            mean_rms = all_rms[stave_idx]
            std_rms = all_rms_std[stave_idx]
            f.write(f"{stave}: {mean_rms:.2f} e$^-$ $\\pm$ {std_rms:.2f} e$^-$\n")

    # print(f"RMS values saved to {rms_file}")
    print(f"Temporal noise values saved to {rms_file}")

    sys.exit(0)
# End of script
