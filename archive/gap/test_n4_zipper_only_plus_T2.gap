# Test the n=4 claim:  Z_n^zip / <T_2> ?= BMW_n at n = 4.
#
# Concretely: compute dim Q(q,t)[B_4] / <R_22, R_23, R_24, R_25>, where
#   R_22 = level-2 zipper at j=1
#   R_23 = level-2 zipper at j=2
#   R_24 = level-3 zipper at j=1   (the only level-3 zipper at n=4)
#   R_25 = level-2 untwisting at j=1
# plus braid + invertibility for B_4.
#
# If dim = 105 = dim BMW_4, then X_3 in <Z_2, T_2, Z_3> (and X_4 in it too),
# verifying my earlier claim that the level-3 zipper R_24 substitutes for
# R_28 in the cofactor.  If dim > 105, the claim is false at n=4.

LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetInfoLevel(InfoGBNPTime, 1);
SetRecursionTrapInterval(50000);

n := 4;
pStr := "2^31-1";
p := EvalString(pStr);
R := GF(p);
q_int := 5;
t_int := 7;
q := q_int * One(R);
t := t_int * One(R);

generator_names := Concatenation(
    List([1..n-1], i -> Concatenation("x", String(i))),
    List([1..n-1], i -> Concatenation("y", String(i)))
);
F := FreeAssociativeAlgebraWithOne(R, generator_names);
gens := GeneratorsOfAlgebra(F);
one := gens[1];
xs := gens{[2..n]};
ys := gens{[n+1..2*n-1]};

# Sliders S[i][j]
S := List([1..n], i -> []);
for j in [1..n-1] do
  S[2][j] := q*ys[j] + (1-q)*one - xs[j];
od;
for i in [3..n] do
  for j in [1..n-i+1] do
    yProd := ys[j];
    for k in [j+1..j+i-2] do yProd := yProd * ys[k]; od;
    xProd := xs[j];
    for k in [j+1..j+i-2] do xProd := xProd * xs[k]; od;
    S[i][j] := (q^(i-1)*yProd - xProd) * S[i-1][j];
  od;
od;
X3 := S[3][1];
X3np := GP2NPList([X3])[1];

# Build relations
inv_rels := [];
for i in [1..n-1] do
  Add(inv_rels, xs[i]*ys[i] - one);
  Add(inv_rels, ys[i]*xs[i] - one);
od;
braid_rels := [];
for i in [1..n-2] do
  for j in [i+1..n-1] do
    if j - i >= 2 then
      Add(braid_rels, xs[i]*xs[j] - xs[j]*xs[i]);
      Add(braid_rels, ys[i]*ys[j] - ys[j]*ys[i]);
      Add(braid_rels, ys[i]*xs[j] - xs[j]*ys[i]);
      Add(braid_rels, xs[i]*ys[j] - ys[j]*xs[i]);
    else
      Add(braid_rels, xs[i]*xs[j]*xs[i] - xs[j]*xs[i]*xs[j]);
      Add(braid_rels, ys[i]*ys[j]*ys[i] - ys[j]*ys[i]*ys[j]);
    fi;
  od;
od;
for i in [1..n-2] do
  Add(braid_rels, xs[i]*xs[i+1]*ys[i]   - ys[i+1]*xs[i]*xs[i+1]);
  Add(braid_rels, ys[i]*ys[i+1]*xs[i]   - xs[i+1]*ys[i]*ys[i+1]);
  Add(braid_rels, xs[i]*ys[i+1]*ys[i]   - ys[i+1]*ys[i]*xs[i+1]);
  Add(braid_rels, ys[i]*xs[i+1]*xs[i]   - xs[i+1]*xs[i]*ys[i+1]);
od;

# Zipper relations: level-2 (j=1,2) and level-3 (j=1)
zipper_rels := [];
# Level-2 zippers, j in [1..n-2] = [1,2]
for j in [1..n-2] do
  Add(zipper_rels,
      (q*ys[j+1] + (1-q)*one - xs[j+1]
       - q*ys[j]*ys[j+1] + xs[j]*xs[j+1]) * S[2][j]);
od;
# Level-3 zipper, j in [1..n-i] with i=3 = [1]
for i in [3..n-1] do
  for j in [1..n-i] do
    yProd_after := ys[j+1];
    for k in [j+2..j+i-1] do yProd_after := yProd_after * ys[k]; od;
    yProd_full := ys[j] * yProd_after;
    xProd_after := xs[j+1];
    for k in [j+2..j+i-1] do xProd_after := xProd_after * xs[k]; od;
    xProd_full := xs[j] * xProd_after;
    Add(zipper_rels,
        (q^(i-1)*yProd_after - xProd_after
         - q^(i-1)*yProd_full + xProd_full) * S[i][j]);
  od;
od;

# Level-2 untwist at j=1 only (others recovered by braid conjugation)
untwist_rels := [(xs[1] - t*one) * S[2][1]];

Print("Relation counts:\n");
Print("  inv          : ", Length(inv_rels), "\n");
Print("  braid        : ", Length(braid_rels), "\n");
Print("  zipper       : ", Length(zipper_rels), " (level-2 x ", n-2,
      ", level-3 x ", n-3, ")\n");
Print("  T_2 only     : ", Length(untwist_rels), "\n");

# Pre-reduce by invertibility
inv_GB := GP2NPList(inv_rels);
zipper_np := Filtered(List(GP2NPList(zipper_rels),
                           pp -> StrongNormalFormNP(pp, inv_GB)),
                      pp -> pp <> []);
untwist_np := Filtered(List(GP2NPList(untwist_rels),
                            pp -> StrongNormalFormNP(pp, inv_GB)),
                       pp -> pp <> []);
braid_np := GP2NPList(braid_rels);

I := Concatenation(inv_GB, braid_np, zipper_np, untwist_np);
Print("Total inputs to SGrobner: ", Length(I), "\n\n");

Print("Running SGrobner ... (this may take a while at n=4)\n");
t_start := Runtime();
B := SGrobner(I);
t_elapsed := Runtime() - t_start;
Print("Done in ", t_elapsed, " ms.\n");
Print("Grobner basis size: ", Length(B), "\n");

nGen := 2*(n-1);
L := LMonsNP(B);
growth := DetermineGrowthQA(L, nGen, false);
Print("Growth (Gel'fand-Kirillov dimension proxy): ", growth, "\n");

if growth = 0 then
  dimQA := DimQA(B, nGen);
  Print("Quotient dimension DimQA(B, ", nGen, ") = ", dimQA, "\n");
  Print("Compare: dim BMW_", n, " = (2n-1)!! = ", Product([1,3..2*n-1]), "\n");
fi;

# Reduce X_3 mod GB
X3_red := StrongNormalFormNP(X3np, B);
X3_is_zero := (X3_red = []) or (Length(X3_red) >= 1 and X3_red[1] = []);
if X3_is_zero then
  Print("X_3 reduces to 0 modulo GB  ==>  X_3 in <Z_2, Z_3, T_2> at n=4.\n");
else
  Print("X_3 mod GB has ", Length(X3_red[1]), " nonzero monomials.\n");
fi;
