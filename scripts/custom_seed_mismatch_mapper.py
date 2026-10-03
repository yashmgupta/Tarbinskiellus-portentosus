from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np, json, csv, time

def run_mapping(root=None, k=19, max_mm=12, max_candidates=80):
    root=Path(root or Path.cwd())
    if not (root/'raw_data').exists() and (root.parent/'raw_data').exists():
        root=root.parent
    raw=root/'raw_data'; out=root/'derived_data'; out.mkdir(exist_ok=True)
    datasets=[
      ('BAN108',raw/'TSN20260330-010-00157_TP-mito_20260331-BAN108_D06_D06_i29.R1.fq','fastq',raw/'TSN20260330-010-00157_TP-mito_20260331-BAN108_D06_D06_i29.R2.fq','fastq'),
      ('BAN72',raw/'TSN20260330-010-00157_TP-mito_20260403-BAN72_A05_A05_i27.R1.fasta','fasta',raw/'TSN20260330-010-00157_TP-mito_20260403-BAN72_A05_A05_i27.R2.fasta','fasta')]
    rec=SeqIO.read(raw/'sequence.gb','genbank'); ref=str(rec.seq).upper(); L=len(ref); ref2=ref+ref
    idx=defaultdict(list)
    for i in range(L):
        mer=ref2[i:i+k]
        if len(mer)==k and 'N' not in mer: idx[mer].append(i)
    def revcomp(s): return str(Seq(s).reverse_complement())
    def candidates(q):
        n=len(q); seed_pos=sorted(set([0,n//5,2*n//5,3*n//5,4*n//5,max(0,n-k)])); votes=Counter()
        for p in seed_pos:
            mer=q[p:p+k]
            if len(mer)<k or 'N' in mer: continue
            for hit in idx.get(mer,[]): votes[(hit-p)%L]+=1
        return [s for s,v in votes.most_common(max_candidates)]
    def hamming(q,start):
        target=ref2[start:start+len(q)] if start+len(q)<=len(ref2) else (ref2+ref)[start:start+len(q)]
        return sum(a!=b and a!='N' for a,b in zip(q,target))
    def map_read(seq):
        seq=seq.upper(); best=None; ties=0
        for strand,q in [('+',seq),('-',revcomp(seq))]:
            for st in candidates(q):
                mm=hamming(q,st)
                if best is None or mm<best[0]: best=(mm,st,strand,len(q)); ties=1
                elif mm==best[0]: ties+=1
        if best is None or best[0]>max_mm: return None
        return best+(ties,)
    def add_cov(cov,st,n):
        end=st+n
        if end<=L: cov[st:end]+=1
        else: cov[st:L]+=1; cov[0:end-L]+=1
    all_cov=np.zeros(L,dtype=np.int64); stats={}; maprows=[]
    for dname,r1,fmt1,r2,fmt2 in datasets:
        cov=np.zeros(L,dtype=np.int64); counts=Counter(); mm_hist=Counter(); tie_hist=Counter(); n_pairs=0; t=time.time()
        for a,b in zip(SeqIO.parse(r1,fmt1),SeqIO.parse(r2,fmt2)):
            n_pairs+=1
            for mate,rr in [(1,a),(2,b)]:
                counts['total_reads']+=1; m=map_read(str(rr.seq))
                if m is None: counts['unmapped']+=1; continue
                mm,st,strand,n,ties=m; counts['mapped']+=1; mm_hist[mm]+=1; tie_hist[ties]+=1
                if ties>1: counts['multimapped_best_tie']+=1
                add_cov(cov,st,n); add_cov(all_cov,st,n); maprows.append((dname,n_pairs,mate,st+1,strand,n,mm,ties))
        stats[dname]={'pairs':n_pairs,**counts,'mean_depth':float(cov.mean()),'median_depth':float(np.median(cov)),'min_depth':int(cov.min()),'max_depth':int(cov.max()),'covered_1x':int((cov>0).sum()),'covered_10x':int((cov>=10).sum()),'mm_hist':dict(mm_hist),'tie_hist':dict(tie_hist),'seconds':time.time()-t}
        np.savetxt(out/f'{dname}_coverage.tsv',np.column_stack((np.arange(1,L+1),cov)),fmt='%d',delimiter='\t',header='position\tdepth',comments='')
    stats['combined']={'total_reads':sum(stats[d]['total_reads'] for d,*_ in datasets),'mapped':sum(stats[d]['mapped'] for d,*_ in datasets),'unmapped':sum(stats[d]['unmapped'] for d,*_ in datasets),'mean_depth':float(all_cov.mean()),'median_depth':float(np.median(all_cov)),'min_depth':int(all_cov.min()),'max_depth':int(all_cov.max()),'covered_1x':int((all_cov>0).sum()),'covered_10x':int((all_cov>=10).sum()),'covered_100x':int((all_cov>=100).sum())}
    np.savetxt(out/'combined_coverage.tsv',np.column_stack((np.arange(1,L+1),all_cov)),fmt='%d',delimiter='\t',header='position\tdepth',comments='')
    with open(out/'mapping_stats.json','w') as f: json.dump(stats,f,indent=2)
    with open(out/'read_mappings.tsv','w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['dataset','pair_no','mate','start_1based','strand','read_length','mismatches','best_hit_ties']); w.writerows(maprows)
    return stats

if __name__=='__main__': print(json.dumps(run_mapping(),indent=2))
