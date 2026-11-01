#!/usr/bin/env python3
"""Deterministic signed/positive/CSR smoke and regression tests.

Optional: OLD_BENCH=/path/to/original/benchmark_four checks exact checksums
against the pre-refactor implementations (compiled with the same compiler flags).
"""
import json,os,random,struct,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OLD=os.getenv('OLD_BENCH','')
def execute(cmd):
 p=subprocess.run([str(v) for v in cmd],capture_output=True,text=True)
 if p.returncode:raise AssertionError(str(cmd)+'\n'+p.stderr)
def fixture(folder,name,n,edges,groups):
 p=folder/name;p.mkdir()
 edges=sorted(edges);g=p/'graph.bin';g.write_bytes(struct.pack('<qq',n,len(edges))+b''.join(struct.pack('<ii',*e) for e in edges))
 E=set(edges);blocks={b:[u for u in range(n) if groups[u]==b] for b in sorted(set(groups))}
 for mode in ('CoRe','Mags'):
  se=[];cp=[];cm=[];active=set()
  for a,AA in blocks.items():
   for b,BB in blocks.items():
    if b<a:continue
    poss=[(u,v) for u in AA for v in BB if u<v]
    m=sum(x in E for x in poss)
    if not poss:continue
    use=(m==len(poss) and m>1) if mode=='CoRe' else (2*m>len(poss)+(a!=b))
    if use:se.append((a,b));active.add((a,b))
  for u in range(n):
   for v in range(u+1,n):
    ab=tuple(sorted((groups[u],groups[v])))
    if ab in active:
     if (u,v) not in E:cm.append((u,v))
    elif (u,v) in E:cp.append((u,v))
  for k,values in [('partition',groups),('se',se),('cp',cp),('cm',cm)]:
   fmt='<i' if k=='partition' else '<ii'
   (p/f'{mode}.{k}.bin').write_bytes(b''.join(struct.pack(fmt,*v) if isinstance(v,tuple) else struct.pack(fmt,v) for v in values))
 return p,g
def main():
 os.environ['QUERY_SEED']='20261101'  # stabilize default-source fixture checks
 for t in ('bfs','dfs','pagerank','coloring'):
  assert (ROOT/'bin'/t).is_file(),'build first: bash querys/scripts/build.sh'
 random_gen=random.Random(111)
 groups=[i//3 for i in range(18)]
 mixed={(u,v) for u in range(18) for v in range(u+1,18)
        if random_gen.random()<(0.83 if groups[u]==groups[v] else 0.77 if abs(groups[u]-groups[v])==1 else 0.09)}
 examples=[('mixed',18,mixed,groups),('clique_140',140,{(u,v) for u in range(140) for v in range(u+1,140)},[u//20 for u in range(140)])]
 cases=0
 with tempfile.TemporaryDirectory(prefix='core_query_test_') as tmp:
  for name,n,E,p in examples:
   base,g=fixture(Path(tmp),name,n,E,p)
   results={}
   for mode in ('CSR','CoRe','Mags'):
    prefix=g if mode=='CSR' else base/mode
    oldrows={}
    if OLD:
     oldfile=base/f'{mode}_original.jsonl';execute([OLD,mode,g,prefix,oldfile,'1','0','0'])
     oldrows={(x['task'],x['kernel']):x for x in map(json.loads,oldfile.read_text().splitlines())}
    for task,exe in [('BFS','bfs'),('DFS','dfs'),('PR','pagerank'),('Coloring','coloring')]:
     out=base/f'{mode}_{task}.jsonl'
     execute([ROOT/'bin'/exe,mode,g,prefix,out,'1','0']+(['0'] if OLD else []))
     rows=list(map(json.loads,out.read_text().splitlines()))
     assert rows and all(r['status']=='ok' for r in rows),rows
     if OLD:
      for r in rows:
       prior=oldrows[(task,r['kernel'])]
       assert r['checksum']==prior['checksum'],(name,mode,task,r['kernel'])
     results[(mode,task)]=rows;cases+=len(rows)
   # One dataset and seed must produce the exact same root in every mode.
   for task in ('BFS','DFS'):
    chosen={results[(mode,task)][0]['source'] for mode in ('CSR','CoRe','Mags')}
    assert len(chosen)==1,(name,task,chosen)
    for mode in ('CSR','CoRe','Mags'):
     for row in results[(mode,task)]:
      assert row['source_method']==('explicit' if OLD else 'seeded_random'),row
      if not OLD: assert row['source_seed']==20261101,row
   if not OLD:
    assert results[('CSR','BFS')][0]['source']==results[('CSR','DFS')][0]['source']
    # Explicitly selecting the generated source must reproduce every DFS/BFS kernel.
    for mode in ('CSR','CoRe','Mags'):
     prefix=g if mode=='CSR' else base/mode
     for task,exe in (('BFS','bfs'),('DFS','dfs')):
      source=results[(mode,task)][0]['source']
      explicit=base/f'{mode}_{task}_explicit.jsonl'
      execute([ROOT/'bin'/exe,mode,g,prefix,explicit,'1','0',str(source)])
      orig=results[(mode,task)]
      chosen=list(map(json.loads,explicit.read_text().splitlines()))
      assert [x['checksum'] for x in orig]==[x['checksum'] for x in chosen],(name,mode,task)
      assert all(x['source_method']=='explicit' for x in chosen)
   for task in ('BFS','Coloring'):
    expected=results[('CSR',task)][0]['checksum']
    for mode in ('CoRe','Mags'):
     for row in results[(mode,task)]:
      assert row['checksum']==expected,(name,mode,task,row['kernel'])
   if name=='clique_140':
    # Confirm root selection responds to seed and truly changes DFS traversal.
    seeds=[]
    for seed in (1,2,3,4,5):
     env=dict(os.environ,QUERY_SEED=str(seed))
     for mode in ('CSR','CoRe','Mags'):
      prefix=g if mode=='CSR' else base/mode
      for task,exe in (('BFS','bfs'),('DFS','dfs')):
       out=base/f'{mode}_{task}_seed{seed}.jsonl'
       run=subprocess.run([str(ROOT/'bin'/exe),mode,str(g),str(prefix),str(out),'1','0','-1'],
                          capture_output=True,text=True,env=env)
       assert run.returncode==0,run.stderr
       rr=list(map(json.loads,out.read_text().splitlines()))
       assert all(x['status']=='ok' and x['source_seed']==seed for x in rr)
       if mode=='CSR' and task=='DFS':seeds.append((seed,rr[0]['source'],rr[0]['checksum']))
       if mode=='CSR' and task=='BFS':csr_bfs_source=rr[0]['source']
       if mode=='CSR' and task=='DFS':assert rr[0]['source']==csr_bfs_source
    assert len({source for seed,source,checksum in seeds})>1,seeds
    assert len({checksum for seed,source,checksum in seeds})>1,seeds
   assert (base/'Mags.cm.bin').stat().st_size>0 or name=='clique_140','test must exercise negative corrections'
   print('[PASS]',name,'all four queries / three representations',flush=True)
 print(f'PASS: {cases} retained kernel executions; random-root parity and old-binary compatibility={bool(OLD)}')
if __name__=='__main__':main()
