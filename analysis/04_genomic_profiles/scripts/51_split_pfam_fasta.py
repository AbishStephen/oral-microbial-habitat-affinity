from pathlib import Path
import heapq

BASE=Path("__PROJECT_ROOT__/study13")
SCR=Path("__WORK_ROOT__/HOMD_convergence_expanded13_work")

inp=BASE/"work/proteins/HOMD_1121_all_proteins.faa"
out=SCR/"pfam/chunks"
out.mkdir(parents=True,exist_ok=True)

N=48
fh=[open(out/f"chunk_{i:02d}.faa","w") for i in range(N)]
heap=[(0,i) for i in range(N)]
heapq.heapify(heap)

name=None
seq=[]
nseq=0
naa=0

def emit(name,seq):
    global nseq,naa
    if name is None:
        return
    s="".join(seq)
    total,i=heapq.heappop(heap)
    fh[i].write(name)
    fh[i].write(s)
    heapq.heappush(heap,(total+len(s.replace("\n","")),i))
    nseq+=1
    naa+=len(s.replace("\n",""))

with open(inp) as f:
    for line in f:
        if line.startswith(">"):
            emit(name,seq)
            name=line
            seq=[]
        else:
            seq.append(line)
    emit(name,seq)

for x in fh:
    x.close()

print("Proteins:",nseq)
print("Residues:",naa)
for total,i in sorted(heap,key=lambda x:x[1]):
    print(f"chunk_{i:02d}",total)
