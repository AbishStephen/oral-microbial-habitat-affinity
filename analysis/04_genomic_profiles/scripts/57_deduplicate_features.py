import csv,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
B=Path('__PROJECT_ROOT__/study13'); H=B/'handoff'; V=B/'08_validation'
with open(H/'X_primary_transferable.tsv') as f:
    r=csv.reader(f,delimiter='\t'); hd=next(r); rows=list(r)
meta=[x[:5] for x in rows]; names=hd[5:]
X=np.array([[float(v) for v in x[5:]] for x in rows],dtype=np.float64)
groups=defaultdict(list)
for j,n in enumerate(names):
    # round only to make text-generated floating representations stable; values themselves retained
    key=np.round(X[:,j],12).tobytes()
    groups[key].append((j,n))
# choose representative preferring biologically interpretable abstraction, then source
priority={'KPATH':0,'KMOD':1,'RHEA':2,'EC':3,'PFAM':4,'COG':5,'X5':6}
keep=[]; eq=[]
for gid,(key,members) in enumerate(groups.items(),1):
    members=sorted(members,key=lambda z:(priority.get(z[1].split(':',1)[0],99),z[1]))
    rep=members[0]; keep.append(rep[0])
    for j,n in members:
        eq.append([f'EQ{gid:05d}',rep[1],n,'YES' if j==rep[0] else 'NO'])
keep=sorted(keep)
with open(H/'feature_equivalence_groups.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['equivalence_group','representative_feature','feature','is_representative']); w.writerows(eq)
with open(H/'X_primary_transferable_DEDUP.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(hd[:5]+[names[j] for j in keep])
    for i,m in enumerate(meta): w.writerow(m+[X[i,j] for j in keep])
sizes=Counter(len(v) for v in groups.values())
layer_removed=Counter()
for v in groups.values():
    if len(v)>1:
        ordered=sorted(v,key=lambda z:(priority.get(z[1].split(':',1)[0],99),z[1]))
        for j,n in ordered[1:]: layer_removed[n.split(':',1)[0]]+=1
print('Original features:',len(names))
print('Unique feature vectors:',len(keep))
print('Exact duplicate columns removed:',len(names)-len(keep))
print('Equivalence group size distribution:',sorted(sizes.items()))
print('Removed by layer:',dict(layer_removed))
# largest redundant groups
largest=sorted((v for v in groups.values() if len(v)>1),key=len,reverse=True)[:20]
print('Largest equivalence groups:')
for v in largest:
    print(len(v),' | '.join(n for _,n in v[:12]))
