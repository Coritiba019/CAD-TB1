import subprocess
import os
import json
MAX_LIMIT = 1000000

R = 1
C = 1
A = 1

runs = []


def run(R, C, A):
    with open("in", 'w') as f:
        f.write("%d\n%d\n%d\n%d\n" % (R, C, A, 7))
    os.system("./sseq < in > outSeq")
    os.system("./spar < in > outPar")
    diffCode = subprocess.call(["diff", "outSeq", "outPar"])
    with open("timeSeq", 'r') as f:
        timeSeq = float(f.read())
    with open("timePar", 'r') as f:
        timePar = float(f.read())
    return {
        "R": R,
        "C": C,
        "A": A,
        "timeSeq": timeSeq,
        "timePar": timePar,
        "error": diffCode != 0
    }


runs = []

while R < MAX_LIMIT:
    C = 1
    while C < MAX_LIMIT:
        A = 1
        while A < MAX_LIMIT:
            if(R*C*A < MAX_LIMIT):
                runs.append(run(R, C, A))
            A *= 100
        C *= 100
    R *= 100

with open("runs.json", 'w') as f:
    f.write(json.dumps(runs, indent=4))
