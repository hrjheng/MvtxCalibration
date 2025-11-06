import os

FEEIDs = [
    [
        [3, 4, 259, 260, 515, 516, 8193, 8194, 8449, 8450, 8705, 8706],
        [5, 261, 517, 8195, 8196, 8207, 8451, 8452, 8463, 8707, 8708, 8719],
    ],
    [
        [4096, 4097, 4352, 4353, 4608, 4609, 8198, 8199, 8454, 8455, 8710, 8711],
        [4098, 4099, 4354, 4355, 4610, 4611, 8200, 8201, 8456, 8457, 8712, 8713],
    ],
    [
        [0, 1, 256, 257, 512, 513, 4100, 4356, 4612, 8202, 8458, 8714],
        [2, 258, 514, 4101, 4102, 4103, 4357, 4358, 4359, 4613, 4614, 4615],
    ],
    [
        [9, 10, 265, 266, 521, 522, 8197, 8203, 8453, 8459, 8709, 8715],
        [11, 267, 523, 8204, 8205, 8206, 8460, 8461, 8462, 8716, 8717, 8718],
    ],
    [
        [6, 7, 262, 263, 518, 519, 4108, 4364, 4620, 8192, 8448, 8704],
        [8, 264, 520, 4109, 4110, 4111, 4365, 4366, 4367, 4621, 4622, 4623],
    ],
    [
        [4104, 4105, 4360, 4361, 4616, 4617, 8208, 8209, 8464, 8465, 8720, 8721],
        [4106, 4107, 4362, 4363, 4618, 4619, 8210, 8211, 8466, 8467, 8722, 8723],
    ],
]

# Setup
runs = [76388, 76389, 76390]
Nsegments = 1
runtype = "calib"
filenames = []
queuetext = ""

for flx in range(6):
    for ep in range(2):
        argument_p = 2000 + flx * 10 + (ep + 1)
        for feeid in FEEIDs[flx][ep]:
            for run in runs:
                for seg in range(Nsegments):
                    filename = "/sphenix/lustre01/sphnxpro/physics/MVTX/{}/{}_mvtx{}-{:08d}-{:04d}.evt".format(
                        runtype, runtype, flx, run, seg
                    )

                    if not os.path.isfile(filename):
                        filename = "/sphenix/lustre01/sphnxpro/fromhpss/physics_2025/MVTX/{}/{}_mvtx{}-{:08d}-{:04d}.evt".format(
                            runtype, runtype, flx, run, seg
                        )

                    if not os.path.isfile(filename):
                        print("File {} does not exist in either location".format(filename))
                        continue

                    prefix = "{}_{:08d}_{:04d}_FEEID{}".format(
                        runtype.capitalize(), run, seg, feeid
                    )
                    queuetext += "{} {} {} {}\n".format(
                        argument_p, filename, feeid, prefix
                    )

with open("./queue.list", "w") as f:
    f.write(queuetext)

f.close()
