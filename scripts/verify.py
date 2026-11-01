#!/usr/bin/env python3
"""Independent full-edge reconstruction and optimal fixed-partition cost audit.
Does not call either C++ builder and does not use their block-count caches.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def verify(graph,prefix,mode):
    t=time.perf_counter();prefix=str(prefix)
    with open(graph,'rb') as f:
        n,m=map(int,np.fromfile(f,'<i8',2));edges=np.fromfile(f,'<i4').reshape(-1,2)
    assert len(edges)==m
    p=np.fromfile(prefix+'.partition.bin','<i4')
    assert len(p)==n and np.all(p>=0) and np.all(p<n)
    k=int(p.max())+1;sz=np.bincount(p,minlength=k)
    assert np.all(sz>0),'expected compact partition labels'
    def load(s): return np.fromfile(prefix+s,'<i4').reshape(-1,2)
    se,cp,cm=load('.se.bin'),load('.cp.bin'),load('.cm.bin')
    assert not len(se) or (se.min()>=0 and se.max()<k and np.all(se[:,0]<=se[:,1]))
    assert len(np.unique(se,axis=0))==len(se)
    def codes(a):
        assert not len(a) or (a.min()>=0 and a.max()<n and np.all(a[:,0]<a[:,1]))
        b=a.astype(np.uint64);return (b[:,0]<<np.uint64(32))|b[:,1]
    original=codes(edges);assert len(np.unique(original))==m
    order=np.argsort(p,kind='stable');offs=np.r_[0,np.cumsum(sz)];pieces=[]
    for a,b in se:
        aa=order[offs[a]:offs[a+1]].astype(np.uint64);bb=order[offs[b]:offs[b+1]].astype(np.uint64)
        if a==b:
            i,j=np.triu_indices(len(aa),1);x=aa[i];y=aa[j]
        else:
            x=np.repeat(aa,len(bb));y=np.tile(bb,len(aa))
        pieces.append((np.minimum(x,y)<<np.uint64(32))|np.maximum(x,y))
    expanded=np.concatenate(pieces) if pieces else np.empty(0,np.uint64)
    cpc,cmc=codes(cp),codes(cm)
    allgen=np.concatenate([expanded,cpc]);assert len(np.unique(allgen))==len(allgen)
    assert len(np.unique(cmc))==len(cmc)
    assert not len(np.setdiff1d(cmc,expanded)),'negative correction outside expansion'
    restored=np.setdiff1d(allgen,cmc)
    assert np.array_equal(restored,np.sort(original)),'reconstructed graph differs from input'
    bc=np.sort(p[edges],axis=1).astype(np.uint64)
    keys,ec=np.unique((bc[:,0]<<np.uint64(32))|bc[:,1],return_counts=True)
    a=(keys>>np.uint64(32)).astype(np.int64);b=(keys&np.uint64(0xffffffff)).astype(np.int64)
    cap=np.where(a==b,sz[a]*(sz[a]-1)//2,sz[a]*sz[b])
    assert np.all(ec<=cap)
    signed=mode=='signed';expected=np.minimum(ec,1+cap-ec) if signed else np.where(ec==cap,1,ec)
    cost=len(se)+len(cp)+len(cm);assert int(expected.sum())==cost
    native=np.minimum(2*ec,np.where(a==b,1,2)+2*(cap-ec)) if signed else np.where(ec==cap,np.where(a==b,1,2),2*ec)
    if not signed: assert len(cm)==0
    metapath=Path(prefix+'.json')
    if metapath.exists():
        j=json.loads(metapath.read_text());assert j['cost']==cost and j['negative_corrections']==len(cm)
        assert j['superedges']==len(se) and j['positive_corrections']==len(cp)
        assert j['native_opt_cost']==int(native.sum())
    return dict(prefix=prefix,mode=mode,n=n,m=m,blocks=k,decoded_edges=len(restored),cost=cost,
        relative_size=cost/m if m else 0,superedges=len(se),positive_corrections=len(cp),negative_corrections=len(cm),
        verified_lossless=True,verified_no_negative=(len(cm)==0),native_opt_cost=int(native.sum()),
        seconds=time.perf_counter()-t,hashes={ext:sha(prefix+ext) for ext in ['.partition.bin','.se.bin','.cp.bin','.cm.bin']})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('graph');p.add_argument('prefix');p.add_argument('mode',choices=['positive','signed']);p.add_argument('--out');a=p.parse_args()
    result=verify(a.graph,a.prefix,a.mode);text=json.dumps(result,indent=2);print(text)
    if a.out:Path(a.out).write_text(text)
