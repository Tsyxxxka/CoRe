#!/usr/bin/env python3
"""Run only the paper's retained BFS/DFS/PageRank/Coloring implementations.

No graph loading or representation layout time enters query_seconds_*.
Results are saved per task/mode/dataset with explicit missing/failure markers.
"""
import csv,json,os,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
PROJECT=Path(os.getenv('MAIN_ROOT',str(ROOT.parent))).resolve()
GRAPH_DIR=Path(os.getenv('GRAPH_DIR',str(PROJECT/'data'))).resolve()
CORE_DIR=Path(os.getenv('CORE_SUMMARY_DIR',str(PROJECT/'results/final/summaries'))).resolve()
MAGS_DIR=Path(os.getenv('MAGS_SUMMARY_DIR',str(PROJECT/'outputs'))).resolve()
OUT=Path(os.getenv('OUT_DIR',str(ROOT/'output'))).resolve()
DATASETS=os.getenv('DATASETS','caida enron brightkite email-euall slashdot dblp amazon0601 cnr-2000 youtube skitter in-2004 eu-2005 eswiki-2013 livejournal hollywood-2011 indochina-2004 uk-2005 it-2004').split()
MODES=os.getenv('MODES','CoRe Mags CSR').split()
TASKS=os.getenv('TASKS','BFS DFS PR Coloring').split()
KERNELS={
 'CSR':{'BFS':['csr'],'DFS':['parent_spine'],'PR':['two_pass'],'Coloring':['first_fit']},
 'CoRe':{'BFS':['expand','sparse_block'],'DFS':['expand','spine'],'PR':['dq5_auto'],'Coloring':['expand','block_bit']},
 'Mags':{'BFS':['expand','runtime_hybrid'],'DFS':['expand_frame'],'PR':['dq5_auto'],'Coloring':['expand','block_u32']},
}
EXE={'BFS':'bfs','DFS':'dfs','PR':'pagerank','Coloring':'coloring'}
FIELDS=['dataset','representation','status','task','kernel','n','m','source','source_method','source_seed','pr_iters','warmups','repeats','read_seconds','graph_metadata_seconds','representation_load_seconds','layout_seconds','verification_seconds','query_seconds_mean','query_seconds_median','query_seconds_runs','checksum','elapsed_wall_seconds','message']
REPEATS=int(os.getenv('REPEATS','3'))
WARMUPS=int(os.getenv('WARMUPS','1'))
SOURCE=int(os.getenv('SOURCE','-1'))
QUERY_SEED=int(os.getenv('QUERY_SEED','20261101'))
TIME_LIMIT=int(os.getenv('TIME_LIMIT_SECONDS','86400'))
FORCE=os.getenv('FORCE','0')=='1'
def prefix(ds,mode):
 if mode=='CoRe':return CORE_DIR/(ds+'_core')
 if mode=='Mags':return MAGS_DIR/(ds+'_mags')
 return GRAPH_DIR/(ds+'.bin')
def placeholder(ds,mode,task,kernel,status,message=''):
 return dict(dataset=ds,representation=mode,task=task,kernel=kernel,status=status,message=message)
def parse(path,ds,mode,task):
 rows=[]
 for line in path.read_text().splitlines():
  if not line.strip():continue
  d=json.loads(line)
  if d['task']!=task or d['representation']!=mode or d['kernel'] not in KERNELS[mode][task]:raise ValueError('unexpected result identity')
  if d['status']=='ok':
   if d['repeats']!=REPEATS or d['warmups']!=WARMUPS:
    raise ValueError('stale run parameters')
   if SOURCE>=0:
    if d.get('source')!=SOURCE or d.get('source_method')!='explicit':
     raise ValueError('stale explicit source')
   elif task in ('BFS','DFS'):
    if d.get('source_method')!='seeded_random' or d.get('source_seed')!=QUERY_SEED:
     raise ValueError('stale random source seed')
  d['dataset']=ds
  if isinstance(d.get('query_seconds_runs'),list):d['query_seconds_runs']=json.dumps(d['query_seconds_runs'],separators=(',',':'))
  rows.append(d)
 found=[x['kernel'] for x in rows]
 if len(found)!=len(set(found)) or set(found)!=set(KERNELS[mode][task]):raise ValueError('incomplete results')
 return rows
def main():
 if QUERY_SEED<0 or QUERY_SEED>=2**64:raise SystemExit('QUERY_SEED must be an unsigned 64-bit integer')
 for mode in MODES:
  if mode not in KERNELS:raise SystemExit('invalid mode: '+mode)
 for task in TASKS:
  if task not in EXE:raise SystemExit('invalid task: '+task)
 for task in TASKS:
  if not (ROOT/'bin'/EXE[task]).is_file():raise SystemExit('build first: bash querys/scripts/build.sh')
 OUT.mkdir(parents=True,exist_ok=True);rows=[]
 for ds in DATASETS:
  graph=GRAPH_DIR/(ds+'.bin')
  for mode in MODES:
   pre=prefix(ds,mode)
   required=[graph]+([Path(str(pre)+'.'+s+'.bin') for s in ('partition','se','cp','cm')] if mode!='CSR' else [])
   missing=[str(f) for f in required if not f.is_file()]
   for task in TASKS:
    raw=OUT/'raw'/mode/task/(ds+'.jsonl');raw.parent.mkdir(parents=True,exist_ok=True)
    if missing:
     rows.extend(placeholder(ds,mode,task,k,'missing','; '.join(missing)) for k in KERNELS[mode][task]);continue
    if raw.exists() and not FORCE:
     try:
      valid=parse(raw,ds,mode,task)
      if all(r['status']=='ok' for r in valid):rows.extend(valid);print('[SKIP]',ds,mode,task,flush=True);continue
     except (ValueError,KeyError,TypeError,FileNotFoundError):pass
    raw.unlink(missing_ok=True)
    cmd=[str(ROOT/'bin'/EXE[task]),mode,str(graph),str(pre),str(raw),str(REPEATS),str(WARMUPS),str(SOURCE)]
    log=OUT/'logs'/mode/task/(ds+'.log');log.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();status='failed';message=''
    try:
     with log.open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=TIME_LIMIT,check=False)
     status='ok' if p.returncode==0 else 'failed';message=f'exit={p.returncode}' if p.returncode else ''
    except subprocess.TimeoutExpired:status='timeout';message=f'exceeded {TIME_LIMIT} seconds'
    except Exception as e:message=repr(e)
    try:
     results=parse(raw,ds,mode,task)
     rows.extend(results)
     if not all(r['status']=='ok' for r in results):status='partial'
    except Exception as e:
     message+='; '+str(e)
     by_kernel={}
     if raw.exists():
      for line in raw.read_text().splitlines():
       try:
        d=json.loads(line)
        if d.get('kernel') in KERNELS[mode][task] and d.get('task')==task:
         d['dataset']=ds
         if isinstance(d.get('query_seconds_runs'),list):d['query_seconds_runs']=json.dumps(d['query_seconds_runs'])
         by_kernel[d['kernel']]=d
       except ValueError:pass
     rows.extend(by_kernel.get(k,placeholder(ds,mode,task,k,status,message)) for k in KERNELS[mode][task])
    print(f'[{status}] {ds} {mode} {task} {time.monotonic()-start:.2f}s',flush=True)
    with (OUT/'results.csv').open('w',newline='') as f:
     w=csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore');w.writeheader();w.writerows(rows)
 with (OUT/'results.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,extrasaction='ignore');w.writeheader();w.writerows(rows)
 print('results:',OUT/'results.csv', 'rows:',len(rows))
if __name__=='__main__':main()
