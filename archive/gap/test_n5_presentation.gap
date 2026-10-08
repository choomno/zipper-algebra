# Test the conjecture at n = 5.
#
# Relation set used here (n=5):
#   inv (2*(n-1)=8 relations)
#   braid + far-comm + length-3 untwistings
#   R_22 = Z_2 (single relation, only at j=1)
#   R_25 = T_2 (single relation)
#   R_28 = T_3 (single relation)
#
# Expected: dim = (2n-1)!! = 945 at n=5, and X_3 := S[3][1] reduces to 0.

LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetInfoLevel(InfoGBNPTime, 1);
SetRecursionTrapInterval(50000);

n := 5;
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

# The three specific relations from the conjecture (still only at j=1).
R22 := (q*ys[2] + (1-q)*one - xs[2] - q*ys[1]*ys[2] + xs[1]*xs[2]) * S[2][1];
R25 := (xs[1] - t*one) * S[2][1];
R28 := (xs[2]*xs[1] - q^2*ys[2]*ys[1]) * S[3][1];

extra_rels := [R22, R25, R28];

Print("Relation counts:\n");
Print("  inv          : ", Length(inv_rels), "\n");
Print("  braid        : ", Length(braid_rels), "\n");
Print("  R_22,R_25,R_28: ", Length(extra_rels), "\n");

inv_GB := GP2NPList(inv_rels);
extras_np := Filtered(List(GP2NPList(extra_rels),
                           pp -> StrongNormalFormNP(pp, inv_GB)),
                      pp -> pp <> []);
braid_np := GP2NPList(braid_rels);

I := Concatenation(inv_GB, braid_np, extras_np);
Print("Total inputs to SGrobner: ", Length(I), "\n\n");

Print("Running SGrobner ... (this may take a while at n=5)\n");
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

X3_red := StrongNormalFormNP(X3np, B);
X3_is_zero := (X3_red = []) or (IsList(X3_red) and Length(X3_red) >= 1 and X3_red[1] = []);
if X3_is_zero then
  Print("X_3 reduces to 0 modulo GB.\n");
else
  Print("X_3 mod GB is NOT zero.\n");
fi;
