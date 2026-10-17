#!/usr/bin/env python3
"""Exploratory RegPrecise HxlR-family motif scan of selected oral intergenic DNA,
with dinucleotide-preserving shuffles and separation by independent motif regulog.
Small, phylogenetically correlated DNA sample = NOT a promoter validation.
"""
from pathlib import Path
import pandas as pd,numpy as np,collections,random,json
from numpy.lib.stride_tricks import sliding_window_view
B=Path('__PROJECT_ROOT__/study13')
O=B/'16_HxlR_reference_phylogeny_20261014'
old=B/'12_HxlR_SnoaL_regulatory_genomics_20261013'
redox=B/'15_HxlR_regulator_target_ligand_20261014'
dna_files={
 'HxlR-SnoaL2':old/'HxlR_SnoaL2_intergenic_all_available.fna',
 'quinone_reductase':redox/'QorR_like_quinone_reductase_intergenic_1perHMT.fna',
 'CoA_disulfide_reductase':redox/'CoA_disulfide_reductase_intergenic_1perHMT.fna',
 'NADH_oxidase':redox/'NADH_oxidase_intergenic_1perHMT.fna',
}
def fasta(path):
 d={};n=None;buf=[]
 for line in open(path):
  if line.startswith('>'):
   if n:d[n]=''.join(buf).upper()
   n=line[1:].strip();buf=[]
  else:buf.append(line.strip())
 if n:d[n]=''.join(buf).upper()
 return d
def rc(seq):return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
m=pd.read_csv(O/'RegPrecise_HxlR11_regulog_annotated_binding_sites.tsv',sep='\t')
lookup=np.full(256,255,dtype=np.uint8)
for j,c in enumerate('ACGT'):lookup[ord(c)]=j
models={}
for (regulog,name),z in m.groupby(['regulog_id','regulator']):
 sites=list(z.sequence)
 if not sites or len({len(x) for x in sites})!=1:continue
 a=np.array([[('ACGT').index(ch) for ch in s] for s in sites])
 count=np.array([(a==j).sum(axis=0) for j in range(4)],dtype=float)
 freq=(count+.5)/(len(sites)+2.)
 weights=np.log2(freq/.25).T
 models[str(regulog)]=(name,weights,len(sites))
print('MODELS',[(k,v[0],v[2]) for k,v in models.items()],flush=True)
def score_one(seq,weights):
 k=len(weights);seq=seq.upper()
 if len(seq)<k or any(x not in 'ACGT' for x in seq):return np.nan,None,None,None
 scores=[]
 for strand,st in [('+',seq),('-',rc(seq))]:
  a=lookup[np.frombuffer(st.encode(),dtype=np.uint8)]
  windows=sliding_window_view(a,k)
  v=weights[np.arange(k)[None,:],windows].sum(axis=1)
  j=int(np.argmax(v))
  scores.append((float(v[j]),int(j),strand,st[j:j+k]))
 return max(scores,key=lambda v:v[0])
def shuffle_dinuc(s,rng):
 # Eulerian traversal using each adjacent nucleotide edge exactly once.
 dd=collections.defaultdict(list)
 for a,b in zip(s[:-1],s[1:]):dd[a].append(b)
 for k in dd:rng.shuffle(dd[k])
 stack=[s[0]];walk=[]
 while stack:
  v=stack[-1]
  if dd[v]:stack.append(dd[v].pop())
  else:walk.append(stack.pop())
 ans=''.join(reversed(walk))
 if len(ans)!=len(s):raise RuntimeError('Dinucleotide shuffler failed')
 return ans
seqs=[]
for module,path in dna_files.items():
 if not path.exists():continue
 d=fasta(path)
 # exact dedup per module only: avoid repeated copies dominating
 unique={s:n for n,s in d.items() if s and all(b in 'ACGT' for b in s)}
 for s,title in unique.items():seqs.append((module,title,s))
print('DNA',len(seqs),collections.Counter(i[0] for i in seqs),flush=True)
out=[];rng=random.Random(20261014);N=99
for i,(module,title,seq) in enumerate(seqs,1):
 shuf=[shuffle_dinuc(seq,rng) for _ in range(N)]
 for mid,(name,weights,nrefs) in models.items():
  observed=score_one(seq,weights)
  if not isinstance(observed,tuple):continue
  null=np.array([score_one(x,weights)[0] for x in shuf],dtype=float)
  p=(sum(null>=observed[0])+1)/(N+1)
  out.append({'module':module,'sample':title,'intergenic_bp':len(seq),'regulog_id':mid,'regulator':name,
   'binding_sites_in_reference':nrefs,'best_pwm_score':observed[0],'match_start_0based':observed[1],
   'match_strand':observed[2],'matched_DNA_word':observed[3],
   'within_sequence_dinucleotide_shuffle_p':p,'null_mean_pwm_max':round(float(null.mean()),3),'null_95pct_pwm_max':round(float(np.quantile(null,.95)),3)})
 if i%10==0:print('SCANNED',i,len(seqs),flush=True)
df=pd.DataFrame(out)
df.to_csv(O/'HxlR_11_regulog_motif_scans_shuffle99_exploratory.tsv',sep='\t',index=False)
summary=df.groupby(['module','regulator','regulog_id']).agg(n_tested=('sample','size'),n_unadjusted_p_le_0_05=('within_sequence_dinucleotide_shuffle_p',lambda x:int((x<=.05).sum())),
 median_p=('within_sequence_dinucleotide_shuffle_p','median'),best_p=('within_sequence_dinucleotide_shuffle_p','min')).reset_index().sort_values(['module','best_p'])
summary.to_csv(O/'HxlR_motif_exploratory_permodule_summary.tsv',sep='\t',index=False)
print('SUMMARY\n'+summary.to_string(index=False),flush=True)
(O/'HxlR_regulatory_motif_scan_caveat.json').write_text(json.dumps({'status':'exploratory','shuffles_per_region':N,
'model_count':len(models),'region_count_after_sequence_dedup':len(seqs),
'caution':'No multiple-hypothesis-corrected binding motifs, no direct DNA binding demonstrated. Site families are computationally inferred, and numerous sequences are phylogenetically correlated. Do NOT assign chemicals or call promoter from matches.'},indent=2))
