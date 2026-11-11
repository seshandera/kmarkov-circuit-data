"""
N-scaling of the charge-1 PSWAP operator-entanglement plateau (Referee 2, pt 4).
Cheap because the orbit is only 2N^2. Uses the k=1 sym E[R]=4 rule (small run
lengths -> small layer-matrix cache); the plateau is rule-independent. Half-ring
cut with one boundary next to the initial sigma_z at site 0.

Goal: is the OSEE plateau O(1), O(log N), or O(N)? Compare to log2(2N^2) bound
and to volume law ~N/2.
"""
import numpy as np, scipy.sparse as sp, time

THETA = np.pi/300
INIT_SITE = 0

def pair_T(phi):
    c, s = np.cos(phi), np.sin(phi)
    U = np.array([[1,0,0,0],[0,c,1j*s,0],[0,1j*s,c,0],[0,0,0,1]], dtype=complex)
    pa=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.array([[1,0],[0,-1]])]
    T={}
    for p0 in range(4):
        for p1 in range(4):
            P=np.kron(pa[p0],pa[p1]); Pe=U.conj().T@P@U; out=[]
            for q0 in range(4):
                for q1 in range(4):
                    Q=np.kron(pa[q0],pa[q1]); co=np.trace(Q.conj().T@Pe)/4
                    if abs(co)>1e-12: out.append((q0*4+q1,co))
            T[p0*4+p1]=out
    return T

def apply_pair(state,a,b,T):
    res={}
    for s,c0 in state.items():
        for o,c2 in T[s[a]*4+s[b]]:
            q0,q1=divmod(o,4); ls=list(s); ls[a]=q0; ls[b]=q1
            t=tuple(ls); res[t]=res.get(t,0)+c0*c2
    return {k:v for k,v in res.items() if abs(v)>1e-12}

def full_layer(state,offset,T,N):
    cur=state
    for j in range(N//2):
        cur=apply_pair(cur,(offset+2*j)%N,(offset+2*j+1)%N,T)
    return cur

def build_basis(N):
    Tref=pair_T(np.pi/4); z0=tuple([3 if i==INIT_SITE else 0 for i in range(N)])
    orbit={z0}; changed=True
    while changed:
        changed=False
        for s in list(orbit):
            for off in(0,1):
                for st2 in full_layer({s:1.0+0j},off,Tref,N):
                    if st2 not in orbit: orbit.add(st2); changed=True
    basis=sorted(orbit); return basis,{s:i for i,s in enumerate(basis)}

def build_layer_matrix(basis,idx,offset,phi,N):
    T=pair_T(phi); K=len(basis); rows,cols,data=[],[],[]
    for j,s in enumerate(basis):
        for st2,c in full_layer({s:1.0+0j},offset,T,N).items():
            rows.append(idx[st2]); cols.append(j); data.append(c)
    return sp.csr_matrix((data,(rows,cols)),shape=(K,K),dtype=complex)

def bipartition_indices(basis,N,Lsites):
    Ls=sorted(Lsites); Rs=sorted(set(range(N))-set(Lsites))
    lkey,rkey={},{}; row=np.empty(len(basis),int); col=np.empty(len(basis),int)
    for i,s in enumerate(basis):
        a=tuple(s[j] for j in Ls); b=tuple(s[j] for j in Rs)
        if a not in lkey: lkey[a]=len(lkey)
        if b not in rkey: rkey[b]=len(rkey)
        row[i]=lkey[a]; col[i]=rkey[b]
    return row,col,len(lkey),len(rkey)

def osee(v,row,col,nrow,ncol):
    v=v/np.sqrt((np.abs(v)**2).sum())
    M=np.zeros((nrow,ncol),complex); M[row,col]=v
    sv=np.linalg.svd(M,compute_uv=False); p=sv**2; p=p[p>1e-14]
    return float(-(p*np.log2(p)).sum())

if __name__=="__main__":
    Ns=[20,30,40,50,60]
    N_STEPS=3000; SAMPLE=150
    print(f"{'N':>4} {'K=2N^2':>7} {'halves':>7} {'OSEE_plateau':>12} {'log2(2N^2)':>11} {'N/2':>5} {'steps':>6} {'t(s)':>6}")
    out={}
    for N in Ns:
        t0=time.time()
        basis,idx=build_basis(N); K=len(basis)
        row,col,nrow,ncol=bipartition_indices(basis,N,set(range(0,N//2)))
        cache={}; rng=np.random.default_rng(7)
        v=np.zeros(K,complex); v[idx[tuple([3 if i==INIT_SITE else 0 for i in range(N)])]]=1.0
        S=[]; Ls=[]
        for t in range(1,N_STEPS+1):
            RA=int(rng.geometric(0.25)); RB=int(rng.geometric(0.25))  # k=1 sym E[R]=4
            for (off,R) in ((0,RA),(1,RB)):
                key=(off,R)
                if key not in cache: cache[key]=build_layer_matrix(basis,idx,off,R*THETA,N)
                v=cache[key]@v
            if t%SAMPLE==0: S.append(osee(v,row,col,nrow,ncol)); Ls.append(t)
        S=np.array(S); plateau=S[len(S)//2:].mean()
        out[N]=(S,np.array(Ls))
        print(f"{N:4d} {2*N*N:7d} {nrow:7d} {plateau:12.3f} {np.log2(2*N*N):11.2f} {N/2:5.0f} "
              f"{N_STEPS:6d} {time.time()-t0:6.0f}")
    np.savez("operator_entanglement/osee_nscaling.npz",
             Ns=Ns, **{f"S_{N}":out[N][0] for N in Ns}, **{f"L_{N}":out[N][1] for N in Ns})
    print("saved operator_entanglement/osee_nscaling.npz")
