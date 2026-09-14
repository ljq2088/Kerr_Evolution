"""Rank audit for local jump matrices; no sourced RHS or flux fit."""
import json
from pathlib import Path
import numpy as np
from paper_jump_basis import projected_matrix


def scaled_spectrum(matrix):
    rownorm=np.linalg.norm(matrix,axis=1)
    keep=rownorm>1e-13*np.max(rownorm)
    reduced=matrix[keep]/rownorm[keep,None]
    colnorm=np.linalg.norm(reduced,axis=0)
    if np.any(colnorm==0):raise ValueError('Unconstrained zero column')
    scaled=reduced/colnorm
    singular=np.linalg.svd(scaled,compute_uv=False)
    rank=int(np.sum(singular>singular[0]*1e-10))
    return dict(equations=matrix.shape[0],nonzero_equations=int(sum(keep)),unknowns=matrix.shape[1],rank=rank,
      singular_values=singular.tolist(),condition_number=float(singular[0]/singular[-1]) if rank==matrix.shape[1] else None)


def run(m,L,q,a):
    A,degrees,labels=projected_matrix(20.,a,m,L,L,q)
    common=np.zeros(A.shape[:-1],bool);common[:,:,0]=common[:,:,1]=True
    corrected=common.copy()
    if abs(m)==1:
        corrected[:,0,1]=False;corrected[:,0,2]=True
    return dict(m=m,ellmax=L,quadrature=q,a=a,standard=scaled_spectrum(A[common]),
       dipole_corrected=scaled_spectrum(A[corrected]),labels=labels),A[corrected]


def main():
    rows=[];matrices={}
    for a,m,L,q in ((0.,1,1,18),(.8771530275949366,1,1,18),(.8771530275949366,1,1,24),
                    (.8771530275949366,2,2,18),(.8771530275949366,1,4,18),(.8771530275949366,1,4,24)):
        row,matrix=run(m,L,q,a);rows.append(row);matrices[(a,m,L,q)]=matrix
        print(json.dumps(row),flush=True)
    for L in (1,4):
        a=.8771530275949366;lo=matrices[(a,1,L,18)];hi=matrices[(a,1,L,24)]
        rows.append(dict(ellmax=L,quadrature_change='18 to 24',relative_matrix_max_difference=float(np.max(abs(hi-lo))/np.max(abs(hi)))))
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/paper_jump_rank_audit.json'
    out.write_text(json.dumps(dict(status='homogeneous_jump_matrix_prototype_not_full_reconstruction',rows=rows,
        limitations=['No known Weyl/trace source RHS assembled',
                     'No flux or complete point-source solution computed',
                     'Only three tetrad components implemented for this matrix prototype',
                     'Double-precision angular eigenfunctions and local jets; arbitrary precision remains to be implemented']),indent=2)+'\n')


if __name__=='__main__':main()
