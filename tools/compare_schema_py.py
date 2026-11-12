# Compare every dataset and attribute of the Julia and Python outputs, up to edge order.
# Usage: compare_schema_py.py JULIA_DIR PYTHON_DIR
import h5py, numpy as np, glob, os, sys
jd, pd_ = sys.argv[1], sys.argv[2]
def canon_edges(f):
    e=f['vertex_to_vertex'][()]; c2e=f['cell_to_edge'][()]; em=f['boundary/edge_markers'][()] if 'boundary/edge_markers' in f else None
    pair=[tuple(sorted(r)) for r in e]
    cell_edges=[sorted(pair[i] for i in row) for row in c2e]
    emk=dict(zip(pair,em)) if em is not None else None
    return pair, cell_edges, emk
tot=0; problems={}
for y in sorted(glob.glob(pd_+'/*.yaml')):
    n=os.path.basename(y)[:-5]; jp=f'{jd}/{n}.h5'; pp=f'{pd_}/{n}.h5'
    if not (os.path.exists(jp) and os.path.exists(pp)): continue
    a=h5py.File(jp); b=h5py.File(pp); tot+=1
    bad=[]
    ka=set(); kb=set()
    a.visititems(lambda k,o: ka.add(k) if isinstance(o,h5py.Dataset) else None)
    b.visititems(lambda k,o: kb.add(k) if isinstance(o,h5py.Dataset) else None)
    if ka!=kb: bad.append(f'dataset names differ: {sorted(ka^kb)}')
    for k in sorted(ka&kb):
        if k in ('vertex_to_vertex','cell_to_edge','boundary/edge_markers'): continue
        x,z=a[k][()],b[k][()]
        if x.shape!=z.shape or x.dtype!=z.dtype: bad.append(f'{k}: shape/dtype {x.shape}{x.dtype} vs {z.shape}{z.dtype}'); continue
        if x.dtype.kind=='f':
            if not np.allclose(x,z,rtol=1e-11,atol=1e-12): bad.append(f'{k}: values differ (max {np.abs(x-z).max():.2e})')
        elif not np.array_equal(x,z): bad.append(f'{k}: values differ')
    if 'vertex_to_vertex' in ka and 'vertex_to_vertex' in kb:
        pa,ca,ma=canon_edges(a); pb,cb,mb=canon_edges(b)
        if sorted(pa)!=sorted(pb): bad.append('edge set differs')
        if ca!=cb: bad.append('cell edge sets differ')
        if ma!=mb: bad.append('edge markers differ')
        la=[tuple(r) for r in a['vertex_to_vertex'][()]]; lb=[tuple(r) for r in b['vertex_to_vertex'][()]]
        if la!=lb: bad.append('edge ORDER or direction differs (set equal)' if sorted(pa)==sorted(pb) else 'edge list differs')
    for k in sorted(set(a.attrs)|set(b.attrs)):
        va,vb=a.attrs.get(k),b.attrs.get(k)
        va=va.decode() if isinstance(va,bytes) else va; vb=vb.decode() if isinstance(vb,bytes) else vb
        if not np.array_equal(va,vb): bad.append(f'attr {k}: {va!r} vs {vb!r}')
    for g in ('periodicity',):
        if g in a:
            for k in set(a[g].attrs)|set(b[g].attrs):
                if a[g].attrs.get(k)!=b[g].attrs.get(k): bad.append(f'{g} attr {k} differs')
    for m in bad: problems.setdefault(m.split(':')[0].split(' (')[0] if not m.startswith('attr') else m,[]).append(n)
print(tot,'configs compared')
for m,ns in problems.items(): print(f'  {m}  -> {len(ns)} configs: {", ".join(ns[:6])}{" ..." if len(ns)>6 else ""}')
if not problems: print('  all datasets and attributes agree')
