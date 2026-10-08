LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 2);
SetInfoLevel(InfoGBNPTime, 2);
SetRecursionTrapInterval(50000);


###############################################################################
# USER PARAMETERS
###############################################################################
#
# Edit the values in this block, then run
#   gap -o <workspace> -K <workspace> zipper.gap
# from the braid-GAP/ directory.
#
# Number of strands.
n := 4;

# The "true" coefficient ring is the function field Q(q,t), but rational
# arithmetic is too expensive; specializing to GF(p) at chosen integer values is
# the computational substitute.
#
# The coefficient ring is GF(p) for the prime p below; q and t are the residue
# classes of q_int and t_int. Run at multiple primes/parameters to validate, or
# for lifting coefficients to Q(q,t).
pStr := "2^31-1";
q_int := 5;
t_int := 7;

# Type of untwisting relations to add. The value is one of:
#   "bigelow"        - Bigelow's conjectured form: 2-strand (xs[j] - t)X_2,
#                      i-strand (reverse_x - q^(i-1)*reverse_y)X_i.
#   "indeterminate"  - Same 2-strand relation as bigelow (eigenvalue t), but
#                      the i-strand eigenvalue is s^(i-1) where
#                      s = s_int*One(R) is an independent parameter. The
#                      s_int parameter below is used only in this case.
#   "silly_question" - Bigelow's i-strand relation, but the 2-strand relation
#                      is weakened to (xs[j]^2 - t^2)X_2 (eigenvalue +/- t).
#   "none"           - Omit the untwisting relations entirely. Quotient may
#                      potentially be infinite-dimensional in this case.
untwisting_type := "silly_question";
s_int := 11;

# Below is a parameter that determines whether to compute the full trace
# (for instance, for the purpose of expressing each basis element as a linear
# combination of the input relations) or just a bare Grobner basis.
# Tracing is dramatically heavier in time, RAM, and disk. Set either:
#   traced := true
#   traced := false
traced := false;


###############################################################################
# END USER PARAMETERS, START AUXILIARY (read, ignoreable)
###############################################################################

# Validate untwisting_type and throw error if invalid.
if not (untwisting_type in
        [ "bigelow", "indeterminate", "silly_question", "none" ]) then
  Error("invalid untwisting_type: \"", untwisting_type,
        "\" (expected \"bigelow\", \"indeterminate\", \"silly question\", ",
        "or \"none\")");
fi;

# Set control flow for traced/untraced computation
# and strings for output file name
if traced then
  mode_str := "traced";
  out_str := "trace";
else
  mode_str := "untraced";
  out_str := "untraced";
fi;

# Generate parameter suffix used in output filenames. For "indeterminate"
# untwisting we have the additional s_int parameter; for "bigelow" and "none"
# only (q,t) appear.
if untwisting_type = "indeterminate" then
  params_str := Concatenation(String(q_int), "-", String(t_int), "-",
                              String(s_int));
else
  params_str := Concatenation(String(q_int), "-", String(t_int));
fi;

# logfile name encodes (mode, untwisting_type, n, p, params_str)
LogTo(Concatenation("logs_and_traces/logfile_zipper_BKL_", untwisting_type,
                    "-", mode_str, "-", String(n), "-(", pStr, ")-",
                    params_str, ".txt"));


###############################################################################
# END AUXILIARY, START ALGEBRA PRESENTATION
###############################################################################


# Field parameters
p := EvalString(pStr);
R := GF(p);
q := q_int * One(R);
t := t_int * One(R);
s := s_int * One(R);


# Define generator names for the BKL (Birman-Ko-Lee / band) presentation.
# Each band a_{ts}, 1 <= s < t <= n, gets a generator; ys are inverses.
# Generator order is lex by (t,s): (2,1), (3,1), (3,2), (4,1), (4,2), (4,3), ...
# x's first, then y's in the same order. Total 2 * Binomial(n,2) generators.
n_pairs := n*(n-1)/2;
generator_names := [];
for t in [2..n] do
  for s in [1..t-1] do
    Add(generator_names, Concatenation("a_", String(t), "_", String(s)));
  od;
od;
for t in [2..n] do
  for s in [1..t-1] do
    Add(generator_names, Concatenation("y_", String(t), "_", String(s)));
  od;
od;
F := FreeAssociativeAlgebraWithOne(R, generator_names);
generators := GeneratorsOfAlgebra(F);
one := generators[1];

# Build 2D arrays xs[t][s] = a_{ts}, ys[t][s] = a_{ts}^{-1}. The flat list
# generator_names runs lex by (t,s), so we walk it in the same order.
xs := List([1..n], i -> []);
ys := List([1..n], i -> []);
gen_idx := 2;   # skip identity at position 1
for t in [2..n] do
  for s in [1..t-1] do
    xs[t][s] := generators[gen_idx];
    gen_idx := gen_idx + 1;
  od;
od;
for t in [2..n] do
  for s in [1..t-1] do
    ys[t][s] := generators[gen_idx];
    gen_idx := gen_idx + 1;
  od;
od;


# build_relations: build the four relation groups in any coefficient ring,
# in the BKL (band-generator) presentation. Used for GF(p) below and over
# Q(q,t,s) for human-readable printing. Returns rec(inv, braid, zipper,
# untwisting). xs_, ys_ are 2D arrays indexed by (t,s) with t > s.
#
# The relations are grouped so that we can reduce the high-degree zipper/untwist
# relations modulo the invertibility ones before handing them to SGrobnerTrace,
# and so that the display can label them.
build_relations := function(one_, xs_, ys_, q_, t_, s_, untw)
  local n_, inv, braid, zipper, untwisting, S,
        t1, s1, t2, s2, r, j, k, ii,
        ax, ay, chain_y, chain_x, chain_y_rev, chain_x_rev;
  n_ := Length(xs_);
  inv := []; braid := []; zipper := []; untwisting := [];


  # Invertibility for each band generator.
  for t1 in [2..n_] do
    for s1 in [1..t1-1] do
      Add(inv, xs_[t1][s1]*ys_[t1][s1] - one_);
      Add(inv, ys_[t1][s1]*xs_[t1][s1] - one_);
    od;
  od;


  # BKL commutativity relations: a_{ts} a_{rq} = a_{rq} a_{ts} whenever the
  # intervals [s,t] and [q,r] are disjoint, or one is strictly nested inside
  # the other. y- and mixed-sign versions included to save reduction work.
  for t1 in [2..n_] do
    for s1 in [1..t1-1] do
      for t2 in [t1..n_] do
        for s2 in [1..t2-1] do
          if t1 < t2 or (t1 = t2 and s1 < s2) then
            if (t1 < s2) or (t2 < s1)            # disjoint
               or (s1 < s2 and t2 < t1)          # (q,r) nested in (s,t)
               or (s2 < s1 and t1 < t2) then     # (s,t) nested in (q,r)
              Add(braid, xs_[t1][s1]*xs_[t2][s2] - xs_[t2][s2]*xs_[t1][s1]);
              Add(braid, ys_[t1][s1]*ys_[t2][s2] - ys_[t2][s2]*ys_[t1][s1]);
              Add(braid, xs_[t1][s1]*ys_[t2][s2] - ys_[t2][s2]*xs_[t1][s1]);
              Add(braid, ys_[t1][s1]*xs_[t2][s2] - xs_[t2][s2]*ys_[t1][s1]);
            fi;
          fi;
        od;
      od;
    od;
  od;


  # BKL triangle relations: for t > s > r,
  #   a_{ts} a_{sr} = a_{tr} a_{ts}  (defines a_{tr} in terms of a_{ts}, a_{sr})
  #   a_{tr} a_{ts} = a_{sr} a_{tr}
  # which together also imply a_{ts} a_{sr} = a_{sr} a_{tr}. We feed the first
  # two plus their y-versions to SGrobner.
  for t1 in [3..n_] do
    for s1 in [2..t1-1] do
      for r in [1..s1-1] do
        Add(braid, xs_[t1][s1]*xs_[s1][r]  - xs_[t1][r]*xs_[t1][s1]);
        Add(braid, xs_[t1][r]*xs_[t1][s1]  - xs_[s1][r]*xs_[t1][r]);
        Add(braid, ys_[t1][s1]*ys_[s1][r]  - ys_[t1][r]*ys_[t1][s1]);
        Add(braid, ys_[t1][r]*ys_[t1][s1]  - ys_[s1][r]*ys_[t1][r]);
      od;
    od;
  od;


  # Views of just the adjacent band generators a_{i+1,i} = sigma_i. Below,
  # sliders/zippers/untwistings use ONLY these (consecutive strands), so
  # they are equivalent to the Artin versions in zipper.gap; the BKL braid
  # structure above provides extra reduction paths for SGrobner.
  ax := List([1..n_-1], k -> xs_[k+1][k]);
  ay := List([1..n_-1], k -> ys_[k+1][k]);


  # Sliders S[number of strands][anchor], consecutive-strand only.
  S := List([1..n_], i -> []);
  for j in [1..n_-1] do
    S[2][j] := q_*ay[j] + (1-q_)*one_ - ax[j];
  od;
  for k in [3..n_] do
    for j in [1..n_-k+1] do
      chain_y := ay[j];
      for ii in [j+1..j+k-2] do chain_y := chain_y * ay[ii]; od;
      chain_x := ax[j];
      for ii in [j+1..j+k-2] do chain_x := chain_x * ax[ii]; od;
      S[k][j] := (q_^(k-1)*chain_y - chain_x) * S[k-1][j];
    od;
  od;


  # Zipper relations (consecutive only).
  for j in [1..n_-2] do
    Add(zipper,
        (q_*ay[j+1] + (1-q_)*one_ - ax[j+1]
         - q_*ay[j]*ay[j+1] + ax[j]*ax[j+1]) * S[2][j]);
  od;
  for k in [3..n_-1] do
    for j in [1..n_-k] do
      chain_y := ay[j+1];
      for ii in [j+2..j+k-1] do chain_y := chain_y * ay[ii]; od;
      chain_x := ax[j+1];
      for ii in [j+2..j+k-1] do chain_x := chain_x * ax[ii]; od;
      Add(zipper,
          (q_^(k-1)*chain_y - chain_x
           - q_^(k-1)*ay[j]*chain_y + ax[j]*chain_x) * S[k][j]);
    od;
  od;


  # Untwisting relations (consecutive only).
  if untw <> "none" then

    for j in [1..n_-1] do
      if untw = "silly_question" then
        Add(untwisting, (ax[j]^2 - t_^2*one_) * S[2][j]);
      else
        Add(untwisting, (ax[j] - t_*one_) * S[2][j]);
      fi;
    od;

    for k in [3..n_] do
      for j in [1..n_-k+1] do
        chain_x_rev := ax[j+k-2];
        chain_y_rev := ay[j+k-2];
        for ii in [2..k-1] do
          chain_x_rev := chain_x_rev * ax[j+k-1-ii];
          chain_y_rev := chain_y_rev * ay[j+k-1-ii];
        od;
        if untw = "indeterminate" then
          Add(untwisting,
              (chain_x_rev - s_^(k-1)*chain_y_rev) * S[k][j]);
        else  # bigelow or silly_question
          Add(untwisting,
              (chain_x_rev - q_^(k-1)*chain_y_rev) * S[k][j]);
        fi;
      od;
    od;

  fi;

  return rec(inv := inv, braid := braid,
             zipper := zipper, untwisting := untwisting);
end;


# GF(p)-algebra presentation
rels := build_relations(one, xs, ys, q, t, s, untwisting_type);
inv_relations        := rels.inv;
braid_relations      := rels.braid;
zipper_relations     := rels.zipper;
untwisting_relations := rels.untwisting;

# Relations preprocessing: pre-reduce zipper / untwisting relations modulo
# invertibility. The invertibility relations already form a Grobner basis on
# their own (their leading monomials don't overlap; all S-polys reduce to zero),
# so we can use them directly in StrongNormalFormNP without an SGrobner call.
# Reducing kills x_i*y_i and y_i*x_i instances, reducing polynomial degree.
# We drop any that reduced to zero (would mean the relation was already implied
# by invertibility alone. This shouldn't happen, but we do it conservatively).
inv_GB := GP2NPList(inv_relations);
zipper_simplified := List(GP2NPList(zipper_relations),
                          p -> StrongNormalFormNP(p, inv_GB));
zipper_simplified := Filtered(zipper_simplified, p -> p <> []);
untwisting_simplified := List(GP2NPList(untwisting_relations),
                              p -> StrongNormalFormNP(p, inv_GB));
untwisting_simplified := Filtered(untwisting_simplified, p -> p <> []);


###############################################################################
# END ALGEBRA PRESENTATION, START AUXILIARY II (read, also ignoreable)
###############################################################################


# Print a legend so the log opens with all the run parameters and the
# convention used for the algebra's variable names.
Print("=== zipper_BKL.gap run (BKL / band-generator presentation) ===\n");
Print("  n  = ", n, " strands\n");
Print("  Coefficient ring : GF(", pStr, ")\n");
Print("  q  = ", q_int, "\n");
Print("  t  = ", t_int, "\n");
if untwisting_type = "indeterminate" then
  Print("  s  = ", s_int, "\n");
fi;
Print("  Untwisting type  : ", untwisting_type, "\n");
Print("  Mode             : ", mode_str, "\n");
Print("\n");
Print("  Variable conventions (BKL band generators):\n");
Print("    a_t_s  =  a_{ts}    (band crossing strand t over strand s, t>s)\n");
Print("    y_t_s  =  a_{ts}^{-1}\n");
Print("    S[I]   =  k-strand slider on ordered strand tuple I (|I|=k>=2)\n");
Print("\n");


# Display the input we're about to feed SGrobner, formatted symbolically.
# We rebuild the same relations over the function field Q(q,t[,s]) and use a
# small NP printer below to format them.
if untwisting_type = "indeterminate" then
  symR := FunctionField(Rationals, ["q", "t", "s"]);
  sym_inds := IndeterminatesOfFunctionField(symR);
  sym_q := sym_inds[1]; sym_t := sym_inds[2]; sym_s := sym_inds[3];
else
  symR := FunctionField(Rationals, ["q", "t"]);
  sym_inds := IndeterminatesOfFunctionField(symR);
  sym_q := sym_inds[1]; sym_t := sym_inds[2]; sym_s := Zero(symR);
fi;
symF := FreeAssociativeAlgebraWithOne(symR, generator_names);
sym_gens := GeneratorsOfAlgebra(symF);
sym_one := sym_gens[1];
sym_xs := List([1..n], i -> []);
sym_ys := List([1..n], i -> []);
gen_idx := 2;
for t1 in [2..n] do
  for s1 in [1..t1-1] do
    sym_xs[t1][s1] := sym_gens[gen_idx];
    gen_idx := gen_idx + 1;
  od;
od;
for t1 in [2..n] do
  for s1 in [1..t1-1] do
    sym_ys[t1][s1] := sym_gens[gen_idx];
    gen_idx := gen_idx + 1;
  od;
od;
sym_rels := build_relations(sym_one, sym_xs, sym_ys,
                            sym_q, sym_t, sym_s, untwisting_type);
sym_inv_GB := GP2NPList(sym_rels.inv);
sym_zipper_simpl := List(GP2NPList(sym_rels.zipper),
                         p -> StrongNormalFormNP(p, sym_inv_GB));
sym_zipper_simpl := Filtered(sym_zipper_simpl, p -> p <> []);
sym_untwisting_simpl := List(GP2NPList(sym_rels.untwisting),
                             p -> StrongNormalFormNP(p, sym_inv_GB));
sym_untwisting_simpl := Filtered(sym_untwisting_simpl, p -> p <> []);

# pretty_coeff_str: render a function-field element. Wrap in parens only if
# the printed form has an inner +/- (i.e. it's a multi-term polynomial).
pretty_coeff_str := function(c)
  local s, i;
  s := String(c);
  for i in [2..Length(s)] do
    if s[i] = '+' or s[i] = '-' then
      return Concatenation("(", s, ")");
    fi;
  od;
  return s;
end;

# pretty_monomial_str: render a monomial (list of generator indices) using
# the names in generator_names (a_t_s / y_t_s) with run-length compression
# for repeated factors.
pretty_monomial_str := function(mon)
  local parts, count, name, i;
  if Length(mon) = 0 then return ""; fi;
  parts := [];
  i := 1;
  while i <= Length(mon) do
    count := 1;
    while i + count <= Length(mon) and mon[i + count] = mon[i] do
      count := count + 1;
    od;
    name := generator_names[mon[i]];
    if count > 1 then
      Add(parts, Concatenation(name, "^", String(count)));
    else
      Add(parts, name);
    fi;
    i := i + count;
  od;
  return JoinStringsWithSeparator(parts, "*");
end;

# pretty_print_np: print a symbolic NP polynomial as "term1 + term2 - term3 ...".
# Folds the leading sign of each coefficient into the term separator so we get
# "- q*x1" rather than "+ -q*x1".
pretty_print_np := function(np)
  local mons, coefs, i, cs, ms, term, first;
  if np = [] then Print("0"); return; fi;
  mons := np[1]; coefs := np[2];
  first := true;
  for i in [1..Length(mons)] do
    cs := pretty_coeff_str(coefs[i]);
    ms := pretty_monomial_str(mons[i]);
    if ms = "" then
      term := cs;
    elif cs = "1" then
      term := ms;
    elif cs = "-1" then
      term := Concatenation("-", ms);
    else
      term := Concatenation(cs, "*", ms);
    fi;
    if first then
      Print(term);
      first := false;
    elif Length(term) >= 1 and term[1] = '-' then
      Print(" - ", term{[2..Length(term)]});
    else
      Print(" + ", term);
    fi;
  od;
end;

print_labeled := function(label, np_list)
  local r;
  if Length(np_list) = 0 then return; fi;
  Print("  [", label, "]\n");
  for r in np_list do
    Print("    ");
    pretty_print_np(r);
    Print("\n");
  od;
end;

Print("Input relations (in q, t, s, (1-q)):\n");
print_labeled("invertibility", GP2NPList(sym_rels.inv));
print_labeled("braid",         GP2NPList(sym_rels.braid));
print_labeled("zipper",        GP2NPList(sym_rels.zipper));
print_labeled("untwisting",    GP2NPList(sym_rels.untwisting));
Print("Simplified zipper / untwisting relations ",
      "(post-reduction mod invertibility):\n");
print_labeled("zipper",     sym_zipper_simpl);
print_labeled("untwisting", sym_untwisting_simpl);

###############################################################################
# END AUXILIARY II, START GROBNER CALCULATION
###############################################################################

I := Concatenation(inv_GB,
                   GP2NPList(braid_relations),
                   zipper_simplified,
                   untwisting_simplified);
Print("Total input relations: ", Length(I), "\n\n");

# Run SGrobnerTrace or SGrobner depending on the user parameter `traced`.
if traced then
  Print("Traced option selected.\n",
        "Beginning SGrobnerTrace calculation.\n");
  GBT := SGrobnerTrace(I);
  B := List(GBT, r -> r.pol);   # plain basis for downstream calls
else
  Print("Untraced option selected.\n",
        "Beginning SGrobner calculation (no trace).\n");
  B := SGrobner(I);
fi;

# NOTE: Add a human-readable print of the grobner basis here.
Print("Groebner basis size: ", Length(B), "\n\n");

# Save the result to a parameter-tagged file.
outFile := Concatenation("logs_and_traces/grobner_zipper_BKL_",
                          untwisting_type, "-", out_str, "-",
                          String(n), "-(", pStr, ")-",
                          params_str, ".gap");
if traced then
  PrintTo(outFile,
    "# Saved traced Grobner basis and parameters from zipper.gap.\n",
    "# GBT is the traced basis (list of rec(pol, trace)); trace tuples\n",
    "# [left, inputIdx, right, coeff] reference positions in I.\n",
    "n := ", n, ";\n",
    "pStr := \"", pStr, "\";\n",
    "p := EvalString(pStr);\n",
    "R := GF(p);\n",
    "q := ", Int(q), "*One(R);\n",
    "t := ", Int(t), "*One(R);\n",
    "I := ", I, ";\n",
    "GBT := ", GBT, ";\n");
else
  PrintTo(outFile,
    "# Saved (untraced) Grobner basis and parameters from zipper.gap.\n",
    "# B is the basis as a list of NP polynomials.\n",
    "n := ", n, ";\n",
    "pStr := \"", pStr, "\";\n",
    "p := EvalString(pStr);\n",
    "R := GF(p);\n",
    "q := ", Int(q), "*One(R);\n",
    "t := ", Int(t), "*One(R);\n",
    "I := ", I, ";\n",
    "B := ", B, ";\n");
fi;

###############################################################################
# END GROBNER CALCULATION, START DIMENSION ANALYSIS
###############################################################################

# Determination of dimension of quotient algebra
# First we determine the rate of growth of homogeneous subspaces
# Only leading monomials are required for determining dimension.
L := LMonsNP(B);   # leading monomials
nGen := n*(n-1);   # number of generators (a's + y's = 2 * Binomial(n,2))
growth := DetermineGrowthQA(L, nGen, true);
# Print quotient algebra betti number growth
if growth=0 then
  Print("The quotient algebra by the ideal is finite-dimensional.\n",
  "Its Gel'fand-Kirillov dimension is: ", growth, "\n");
elif not IsString(growth) then
  Print("The quotient algebra is infinite-dimensional of polynomial growth.\n",
  "Its Gel'fand-Kirillov dimension is: ", growth, "\n");
else 
  Print("The quotient algebra is infinite-dimensional of exponential growth.\n",
  "Its Gel'fand-Kirillov dimension is infinite.\n");
fi;

# HilbertSeriesQA(Lm, t, d): values of the Hilbert series up to degree d
# This is a sanity check for comparison to BMW, and is also just kind of neat.
dMax := n*(n-1);
hilb := HilbertSeriesQA(L, nGen, dMax);
hilb_deg := Length(hilb);
dims := List([1..hilb_deg], i -> Sum(hilb{[1..i]}));

# Print Hilbert series result
Print("Hilbert series (degree 0,...,", hilb_deg-1, "): ", hilb, "\n");
Print("Cumulative dims: ", dims, "\n");

if growth=0 then
  Print("Computing quotient dimension.\n");
  dimQA := DimQA(B, nGen);
  Print("Quotient algebra dimension DimQA(B, ",nGen,") = ", dimQA, "\n");
  Print("Compare: dim BMW_", n, " = (2n-1)!! = ",
      Product([1, 3..2*n-1]), "\n");
fi;