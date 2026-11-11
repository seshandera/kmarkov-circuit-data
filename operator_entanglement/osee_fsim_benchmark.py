"""
fSim benchmark: does adding a |11> phase (density-density interaction) move the
operator-entanglement (OSEE) saturation from O(1) (free-fermion PSWAP) to
volume-law (interacting/chaotic)?  Referee 2, point 4/5.

Gate: fSim(theta,phi) = [[1,0,0,0],[0,c,is,0],[0,is,c,0],[0,0,0,e^{i phi}]].
phi=0  -> charge-1 PSWAP (free-fermion / matchgate): sigma_z stays a fermion
          bilinear, orbit = 2N^2, OSEE = O(1).
phi!=0 -> interacting (still U(1)/number-conserving): bilinears couple to
          higher operators, orbit grows past 2N^2, OSEE -> volume-law.

Dense-dict operator evolution (orbit not assumed), deterministic brickwork
(even then odd layer per step), theta=pi/4. Half-ring cut next to sigma_z at
site 0.  Small N only (orbit grows).
"""
import numpy as np, time

THETA = np.pi/4
INIT_SITE = 0

def pair_T(theta, phi):
    c, s = np.cos(theta), np.sin(theta)
    U = np.array([[1,0,0,0],[0,c,1j*s,0],[0,1j*s,c,0],[0,0,0,np.exp(1j*phi)]], dtype=complex)
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

def osee(state,N,Lsites):
    Ls=sorted(Lsites); Rs=sorted(set(range(N))-set(Lsites))
    lkey,rkey={},{}; rows=[]; cols=[]; vals=[]
    tot=np.sqrt(sum(abs(v)**2 for v in state.values()))
    for s,v in state.items():
        a=tuple(s[j] for j in Ls); b=tuple(s[j] for j in Rs)
        if a not in lkey: lkey[a]=len(lkey)
        if b not in rkey: rkey[b]=len(rkey)
        rows.append(lkey[a]); cols.append(rkey[b]); vals.append(v/tot)
    M=np.zeros((len(lkey),len(rkey)),complex)
    for r,c,val in zip(rows,cols,vals): M[r,c]=val
    sv=np.linalg.svd(M,compute_uv=False); p=sv**2; p=p[p>1e-14]
    return float(-(p*np.log2(p)).sum())

def run(N, phi, n_layers=80, sample=4):
    T=pair_T(THETA,phi)
    z0=tuple([3 if i==INIT_SITE else 0 for i in range(N)])
    state={z0:1.0+0j}
    Lsites=set(range(0,N//2))
    S=[]; sizes=[]
    for t in range(1,n_layers+1):
        state=full_layer(state, 0 if t%2==1 else 1, T, N)
        if t%sample==0:
            S.append(osee(state,N,Lsites)); sizes.append(len(state))
    return np.array(S), max(sizes)

if __name__=="__main__":
    print(f"theta=pi/4, deterministic brickwork, half-ring cut next to site 0\n")
    print(f"{'N':>3} {'phi':>6} {'orbit_max':>10} {'2N^2':>6} {'OSEE_sat':>9} {'volume ~N/2':>11}")
    for N in [6,8,10]:
        for phi,lab in [(0.0,"0"),(np.pi/2,"pi/2")]:
            t0=time.time()
            S,omax=run(N,phi, n_layers=max(60,8*N))
            sat=S[len(S)//2:].mean()
            print(f"{N:3d} {lab:>6} {omax:10d} {2*N*N:6d} {sat:9.3f} {N/2:11.1f}   ({time.time()-t0:.0f}s)")
    print("\nExpect: phi=0 -> orbit=2N^2, OSEE flat ~O(1); phi=pi/2 -> orbit>>2N^2, OSEE grows with N.")
