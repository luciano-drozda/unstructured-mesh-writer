# Print the dihedral statistics of every tetrahedral .h5 file in one folder.
# Usage: quality_stats.py FOLDER
# Columns: cell count, smallest and largest dihedral angle (degrees), the
# number of cells below 10 degrees, below 5 degrees and above 170 degrees,
# and the percentage of cells below 10 degrees. It needs h5py and numpy.
import h5py, numpy as np, sys, glob, itertools
def dih(P,C):
    a,b,c=(P[C[:,k]]-P[C[:,0]] for k in (1,2,3))
    S=np.zeros((len(C),4,3)); S[:,1]=0.5*np.cross(c,b); S[:,2]=0.5*np.cross(a,c); S[:,3]=0.5*np.cross(b,a); S[:,0]=-S[:,1:].sum(1)
    m=np.linalg.norm(S,axis=2); out=[]
    for k,l in itertools.combinations(range(4),2):
        out.append(np.degrees(np.arccos(np.clip(-(S[:,k]*S[:,l]).sum(1)/m[:,k]/m[:,l],-1,1))))
    out=np.stack(out,1); return out.min(1),out.max(1)
print(f"{'file':18s} {'cells':>7s} {'min':>6s} {'max':>7s} {'<10deg':>7s} {'<5deg':>6s} {'>170deg':>8s} {'%<10':>6s}")
for f in sorted(glob.glob(sys.argv[1]+'/*.h5')):
    h=h5py.File(f); P=h['points'][()]; C=h['cell_to_vertex'][()]
    mn,mx=dih(P,C)
    print(f"{f.split('/')[-1][:-3]:18s} {len(C):7d} {mn.min():6.2f} {mx.max():7.2f} {(mn<10).sum():7d} {(mn<5).sum():6d} {(mx>170).sum():8d} {100*(mn<10).mean():6.2f}")
