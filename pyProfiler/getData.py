import json

NUM_THREADS = 8

runs = json.load(open("ra.json"))

infos = []

for run in runs:

    speedup = run["timeSeq"] / run["timePar"]
    efficiency = speedup / NUM_THREADS
    info = {
        "R": run["R"],
        "C": run["C"],
        "A": run["A"],
        "Speedup": f'{speedup:.2f}',
        "Efficiency": f'{efficiency:.2f}',
        "Sequential Time": f'{run["timeSeq"]:.6f}s',
        "Parallel Time": f'{run["timePar"]:.6f}s',
    }
    infos.append(info)

json.dump(infos, open("data.json", "w"), indent=4)
