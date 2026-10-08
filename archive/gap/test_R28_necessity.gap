# Test whether the level-3 untwisting R_28 is needed for X_3 = 0.
#
# Hypothesis under test (from user):
#   In the algebra obtained from Q(q,t)[B_3] by adding only the zipper
#   relation R_22 and the level-2 untwisting R_25 (= sigma_1 X_2 = t X_2)
#   -- with no level-3 untwisting and no second level-2 untwist (which is
#   implied by R_25 + braid via conjugation) -- the slider X_3 = S[3][1]
#   already vanishes, and the dimension is 15 = dim BMW_3.
#
# We compute three Grobner bases at n=3 to compare:
#   (A) full Z_3 with all relations (baseline; dim should be 15)
#   (B) drop the level-3 untwisting (only level-2 untwistings + zipper)
#   (C) drop both level-3 untwisting AND zipper (only level-2 untwistings)
#   (D) drop both level-3 untwisting AND zipper AND second level-2 untwist
#       (only R_25 = sigma_1 X_2 = t X_2)
#
# For each: report dim, whether X_3 reduces to 0 modulo the GB.

LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 0);
SetInfoLevel(InfoGBNPTime, 0);
SetRecursionTrapInterval(50000);

n := 3;
pStr := "2^31-1";
p := EvalString(pStr);
R := GF(p);
q_int := 5;
t_int := 7;
q := q_int * One(R);
t := t_int * One(R);

# Free algebra setup
generator_names := Concatenation(
    List([1..n-1], i -> Concatenation("x", String(i))),
    List([1..n-1], i -> Concatenation("y", String(i)))
);
F := FreeAssociativeAlgebraWithOne(R, generator_names);
generators := GeneratorsOfAlgebra(F);
one := generators[1];
xs := generators{[2..n]};;
ys := generators{[n+1..2*n-1]};;

# Build the sliders S[i][j]
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

# The target X_3 = S[3][1]
X3 := S[3][1];
X3_np := GP2NPList([X3])[1];
Print("X_3 has ", Length(X3_np[1]), " monomials over GF(", pStr, ") at q=", q_int, " t=", t_int, "\n\n");

# Build relation groups
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

# Zipper relations (only j=1 at n=3)
zipper_rels := [];
for j in [1..n-2] do
  Add(zipper_rels,
      (q*ys[j+1] + (1-q)*one - xs[j+1]
       - q*ys[j]*ys[j+1] + xs[j]*xs[j+1]) * S[2][j]);
od;

# Level-2 untwistings (j=1..n-1: j=1, j=2 at n=3)
l2_untwist_rels := [];
for j in [1..n-1] do
  Add(l2_untwist_rels, (xs[j] - t*one) * S[2][j]);
od;

# Level-3 untwisting (j=1 only at n=3)
l3_untwist_rels := [];
i := 3;
for j in [1..n-i+1] do
  xProd_rev := xs[j+i-2];
  yProd_rev := ys[j+i-2];
  for k in [2..i-1] do
    xProd_rev := xProd_rev * xs[j+i-1-k];
    yProd_rev := yProd_rev * ys[j+i-1-k];
  od;
  Add(l3_untwist_rels, (xProd_rev - q^(i-1)*yProd_rev) * S[i][j]);
od;

# Level-3 ZIPPER (manually constructed for n=3; would be R_24 at n=4, j=1).
# The polynomial expression is well-defined in B_3, even though formally
# zipper.gap only generates this relation for n >= 4.
# Formula: (q^(i-1)*yProd_after - xProd_after
#          - q^(i-1)*yProd_full + xProd_full) * S[i][j], with i=3, j=1
# yProd_after = y_2, xProd_after = x_2
# yProd_full = y_1 * y_2, xProd_full = x_1 * x_2
l3_zipper_rels := [];
i := 3;
for j in [1..Maximum(n-i, 0)] do
  yProd_after := ys[j+1];
  for k in [j+2..j+i-1] do yProd_after := yProd_after * ys[k]; od;
  yProd_full := ys[j] * yProd_after;
  xProd_after := xs[j+1];
  for k in [j+2..j+i-1] do xProd_after := xProd_after * xs[k]; od;
  xProd_full := xs[j] * xProd_after;
  Add(l3_zipper_rels,
      (q^(i-1)*yProd_after - xProd_after
       - q^(i-1)*yProd_full + xProd_full) * S[i][j]);
od;
# Manually add level-3 zipper for n=3 (computation requires only B_3 vars).
if n = 3 then
  Add(l3_zipper_rels,
      (q^2*ys[2] - xs[2] - q^2*ys[1]*ys[2] + xs[1]*xs[2]) * S[3][1]);
fi;

# Print sizes
Print("Relation group sizes:\n");
Print("  inv         : ", Length(inv_rels), "\n");
Print("  braid       : ", Length(braid_rels), "\n");
Print("  zipper      : ", Length(zipper_rels), "\n");
Print("  l2_untwist  : ", Length(l2_untwist_rels), "\n");
Print("  l3_untwist  : ", Length(l3_untwist_rels), "\n\n");

# Pre-reduce zipper / untwisting by invertibility
inv_GB := GP2NPList(inv_rels);

simplify_mod_inv := function(rels)
  local np, p;
  np := GP2NPList(rels);
  np := List(np, p -> StrongNormalFormNP(p, inv_GB));
  return Filtered(np, p -> p <> []);
end;

zipper_np    := simplify_mod_inv(zipper_rels);
l2_np        := simplify_mod_inv(l2_untwist_rels);
l3_np        := simplify_mod_inv(l3_untwist_rels);
l3_zip_np    := simplify_mod_inv(l3_zipper_rels);
braid_np     := GP2NPList(braid_rels);

# Common base: inv + braid
base_np := Concatenation(inv_GB, braid_np);

# Define the test cases
test_cases := [
  rec(name := "(A) FULL: inv + braid + l2-zip + l2-untw + l3-untw",
      rels := Concatenation(base_np, zipper_np, l2_np, l3_np)),
  rec(name := "(B) DROP l3-untw: inv + braid + l2-zip + l2-untw",
      rels := Concatenation(base_np, zipper_np, l2_np)),
  rec(name := "(C) ONLY l2-untw: inv + braid + l2-untw",
      rels := Concatenation(base_np, l2_np)),
  rec(name := "(D) ONLY R_25: inv + braid + (x1-t)*S[2][1]",
      rels := Concatenation(base_np, [l2_np[1]])),
  rec(name := "(E) USER: inv + braid + l2-zip + R_25 + l3-zip",
      rels := Concatenation(base_np, zipper_np, [l2_np[1]], l3_zip_np)),
  rec(name := "(F) Compare: inv + braid + l2-zip + R_25 + l3-untw",
      rels := Concatenation(base_np, zipper_np, [l2_np[1]], l3_np))
];

# Run each test
nGen := 2*(n-1);
results := [];
for tc in test_cases do
  Print("\n----------------------------------------\n");
  Print("Test ", tc.name, "\n");
  Print("Total input relations: ", Length(tc.rels), "\n");
  Print("Running SGrobner ...\n");
  B := SGrobner(tc.rels);
  Print("Grobner basis size: ", Length(B), "\n");
  L := LMonsNP(B);
  growth := DetermineGrowthQA(L, nGen, false);
  Print("Growth: ", growth, "\n");
  dimQA := fail;
  if growth = 0 then
    dimQA := DimQA(B, nGen);
    Print("Dimension: ", dimQA, "\n");
  fi;
  # Reduce X3 mod GB
  X3_red := StrongNormalFormNP(X3_np, B);
  X3_is_zero := (X3_red = []) or (Length(X3_red) >= 1 and X3_red[1] = []);
  if X3_is_zero then
    Print("X_3 mod GB = 0  ===>  X_3 is in the ideal.\n");
  else
    Print("X_3 mod GB has ", Length(X3_red[1]), " monomials  ===>  X_3 NOT in ideal.\n");
  fi;
  Add(results, rec(name := tc.name, gb_size := Length(B),
                   growth := growth, dim := dimQA,
                   X3_zero := X3_is_zero));
od;

Print("\n\n", "=" )*( 70, "\n");
Print("SUMMARY at n=", n, " over GF(", pStr, "), q=", q_int, " t=", t_int, "\n");
Print("======================================================================\n");
for r in results do
  Print(r.name, "\n");
  Print("  GB size = ", r.gb_size,
        "  growth = ", r.growth,
        "  dim = ", r.dim,
        "  X_3 = 0? ", r.X3_zero, "\n");
od;
Print("Compare: dim BMW_", n, " = (2n-1)!! = ", Product([1, 3..2*n-1]), "\n");
