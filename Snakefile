inputs = ["src/scripts/pipeline.py"]
for i in range(15):
    inputs.extend(
        [
            f"src/data/fiducial_10_10_{i}.h5",
            f"src/data/fiducial_11_10_{i}.h5",
            f"src/data/fiducial_11_11_{i}.h5",
            f"src/data/fiducial_12_10_12_{i}.h5",
            f"src/data/q3_10_10_{i}.h5",
            f"src/data/q3_11_10_{i}.h5",
            f"src/data/q3_11_11_{i}.h5",
            f"src/data/q3_12_10_12_{i}.h5",
            f"src/data/alpha25_10_10_{i}.h5",
            f"src/data/alpha25_11_10_{i}.h5",
            f"src/data/alpha25_11_11_{i}.h5",
            f"src/data/alpha25_12_10_12_{i}.h5",
            f"src/data/alpha5_10_10_{i}.h5",
            f"src/data/alpha5_11_10_{i}.h5",
            f"src/data/alpha5_11_11_{i}.h5",
            f"src/data/alpha5_12_10_12_{i}.h5"
        ]
    )

rule results:
    input:
        inputs
    cache:
        True
    output:
        "src/data/results.hdf"
    script:
        "src/scripts/pipeline.py"