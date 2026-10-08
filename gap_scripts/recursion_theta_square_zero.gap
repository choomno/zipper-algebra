# Determine the minimal polynomial of the lasso twist Theta_i on the slider X_i
# in Zip_4 with ONLY the 2-strand untwisting + zippers (NO full-twist untwisting).
# Calibrate on i=3 (expect to recover the algebra law for X_3), then read off i=4.
# Reported to the user; NOT for the writeup.
LoadPackage("gbnp");
SetInfoLevel(InfoGBNP,1); SetInfoLevel(InfoGBNPTime,0);
SetRecursionTrapInterval(200000);

p:=2^31-1; F:=GF(p); q:=5*One(F); t:=7*One(F);
A:=FreeAssociativeAlgebraWithOne(F,"x1","x2","x3","y1","y2","y3");
gg:=GeneratorsOfAlgebra(A); e:=gg[1];
x1:=gg[2];x2:=gg[3];x3:=gg[4];y1:=gg[5];y2:=gg[6];y3:=gg[7];
xs:=[x1,x2,x3]; ys:=[y1,y2,y3];

rels:=[];
for i in [1..3] do Add(rels, xs[i]*ys[i]-e); Add(rels, ys[i]*xs[i]-e); od;
# far commutativity (1,3)
Add(rels, x1*x3-x3*x1); Add(rels,y1*y3-y3*y1); Add(rels,y1*x3-x3*y1); Add(rels,x1*y3-y3*x1);
# adjacent braid
Add(rels, x1*x2*x1-x2*x1*x2); Add(rels, x2*x3*x2-x3*x2*x3);
Add(rels, y1*y2*y1-y2*y1*y2); Add(rels, y2*y3*y2-y3*y2*y3);
# mixed length-3 (i=1,2)
for i in [1,2] do
  Add(rels, xs[i]*xs[i+1]*ys[i]-ys[i+1]*xs[i]*xs[i+1]);
  Add(rels, ys[i]*ys[i+1]*xs[i]-xs[i+1]*ys[i]*ys[i+1]);
  Add(rels, xs[i]*ys[i+1]*ys[i]-ys[i+1]*ys[i]*xs[i+1]);
  Add(rels, ys[i]*xs[i+1]*xs[i]-xs[i+1]*xs[i]*ys[i+1]);
od;

S2:=List([1..3], j-> q*ys[j]+(1-q)*e-xs[j]);
S31:=(q^2*y1*y2 - x1*x2)*S2[1];
S41:=(q^3*y1*y2*y3 - x1*x2*x3)*S31;

# zipper 2-strand j=1,2
for j in [1,2] do
  Add(rels, (q*ys[j+1]+(1-q)*e-xs[j+1] - q*ys[j]*ys[j+1]+xs[j]*xs[j+1])*S2[j]);
od;
# zipper 3-strand j=1
Add(rels, (q^2*y2*y3 - x2*x3 - q^2*y1*y2*y3 + x1*x2*x3)*S31);
# 2-strand untwisting (sigma-t)^2 [SQUARE-ZERO] ONLY
for j in [1..3] do Add(rels, (xs[j]*xs[j]-2*t*xs[j]+t^2*e)*S2[j]); od;

Print("Computing Grobner basis of Zip_4 (2-strand untwist + zippers, NO full twist)...\n");
t0:=Runtime();
GB:=SGrobner(GP2NPList(rels));
Print("GB done in ", Runtime()-t0, " ms; |GB|=", Length(GB), "\n");
diminfo := DimQA(GB,6);
Print("dim Zip_4^partial = ", diminfo, "\n");

Theta3:=x1*x2*x2*x1;
Theta4:=x1*x2*x3*x3*x2*x1;

NFofGP := function(elt) return StrongNormalFormNP(GP2NPList([elt])[1], GB); end;
IsZeroNF := function(nf) return nf=[] or (IsList(nf) and Length(nf)>=1 and nf[1]=[]); end;

Print("X3 reduces to 0? ", IsZeroNF(NFofGP(S31)), "\n");
Print("X4 reduces to 0? ", IsZeroNF(NFofGP(S41)), "\n");

# minimal polynomial of Theta on X via normal forms + linear dependence
MinPoly := function(Theta, X, maxk)
  local nfs, monset, k, cur, nf, mons, coeffs, i2, basis, vecs, v, idx, m, mat, dep, rk, prevrk;
  nfs:=[]; cur:=X;
  for k in [0..maxk] do
    nf:=NFofGP(cur);
    Add(nfs, nf);
    cur:=cur*Theta;   # Theta^{k+1} X (as GAP elt; reduce at NF time)
  od;
  # collect monomials
  monset:=[];
  for nf in nfs do if not IsZeroNF(nf) then for m in nf[1] do AddSet(monset,m); od; fi; od;
  Print("  (#monomials in NF span = ", Length(monset), ")\n");
  # build coefficient vectors
  vecs:=[];
  for nf in nfs do
    v:=ListWithIdenticalEntries(Length(monset), Zero(F));
    if not IsZeroNF(nf) then
      for i2 in [1..Length(nf[1])] do
        idx:=Position(monset, nf[1][i2]); v[idx]:=nf[2][i2];
      od;
    fi;
    Add(vecs,v);
  od;
  # first k with vecs[1..k+1] dependent: rank stops increasing
  prevrk:=0;
  for k in [1..Length(vecs)] do
    rk:=RankMat(vecs{[1..k]});
    if rk < k then
      # dependence: vecs[1..k] dependent, min poly degree k-1
      dep:=NullspaceMat(TransposedMat(vecs{[1..k]}));  # combos giving 0
      return rec(degree:=k-1, relation:=dep);
    fi;
  od;
  return rec(degree:=fail, note:="no dependence up to maxk");
end;

if not IsZeroNF(NFofGP(S31)) then
  Print("=== Theta_3 on X_3 ===\n");
  r3:=MinPoly(Theta3, S31, 8);
  Print("Theta_3 min-poly degree on X_3 = ", r3.degree, "\n");
  if r3.degree<>fail then Print("  relation (coeffs of sum a_k Theta^k X=0): ", r3.relation, "\n"); fi;
fi;
if not IsZeroNF(NFofGP(S41)) then
  Print("=== Theta_4 on X_4 ===\n");
  r4:=MinPoly(Theta4, S41, 10);
  Print("Theta_4 min-poly degree on X_4 = ", r4.degree, "\n");
  if r4.degree<>fail then Print("  relation: ", r4.relation, "\n"); fi;
fi;
Print("q=",q," q^2=",q^2," q^3=",q^3," q^4=",q^4," t=",t," -t=",-t,"\n");
QUIT_GAP();
