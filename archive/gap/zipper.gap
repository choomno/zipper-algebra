LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 2);
SetInfoLevel(InfoGBNPTime, 2);
SetRecursionTrapInterval(50000);


###############################################################################
# USER PARAMETERS
###############################################################################
#
# Each parameter reads from an environment variable. If the variable is unset
# or empty, the fallback default below is used. Run with (e.g.):
#   N=4 Q=5 T=7 UNTWIST=bigelow TRACED=true gap -o 42G -K 42G zipper.gap
# Override only the values you want to change; the rest fall back to defaults.

ReadParam := function(name, default_)
  if IsBound(GAPInfo.SystemEnvironment.(name))
     and GAPInfo.SystemEnvironment.(name) <> "" then
    return GAPInfo.SystemEnvironment.(name);
  fi;
  return default_;
end;

# Number of strands.
n := Int(ReadParam("N", "3"));

# The "true" coefficient ring is the function field Q(q,t), but rational
# arithmetic is too expensive; specializing to GF(p) at chosen integer values is
# the computational substitute.
#
# The coefficient ring is GF(p) for the prime p below; q and t are the residue
# classes of q_int and t_int. Run at multiple primes/parameters to validate, or
# for lifting coefficients to Q(q,t).
pStr  :=     ReadParam("PSTR", "2^31-1");
q_int := Int(ReadParam("Q",    "5"));
t_int := Int(ReadParam("T",    "7"));

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
untwisting_type :=     ReadParam("UNTWIST", "bigelow");
s_int           := Int(ReadParam("S",       "11"));

# Whether to compute the full trace (express each basis element as a linear
# combination of the input relations) or just a bare Grobner basis. Tracing
# is dramatically heavier in time, RAM, and disk. Set TRACED=true / false
# (case-sensitive).
traced := (ReadParam("TRACED", "false") = "true");


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
LogTo(Concatenation("logs_and_traces/gap_traces/logfile_zipper_", untwisting_type, "-",
                    mode_str, "-", String(n), "-(", pStr, ")-", params_str,
                    ".txt"));


###############################################################################
# END AUXILIARY, START ALGEBRA PRESENTATION
###############################################################################


# Field parameters
p := EvalString(pStr);
R := GF(p);
q := q_int * One(R);
t := t_int * One(R);
s := s_int * One(R);


# Define generator names for a free associative algebra on 2*(n-1) variables
# xs are generators, ys their inverses
# Generator order: x1, x2, ..., x_{n-1}, y1, y2, ..., y_{n-1}.
generator_names := Concatenation(
    List([1..n-1], i -> Concatenation("x", String(i))),
    List([1..n-1], i -> Concatenation("y", String(i)))
);
F := FreeAssociativeAlgebraWithOne(R, generator_names);
generators := GeneratorsOfAlgebra(F);
one := generators[1];
xs := generators{[2..n]};;
ys := generators{[n+1..2*n-1]};;


# build_relations: build the four relation groups in any coefficient ring.
# Used for the GF(p) computation below and a second time over Q(q,t,s) for
# printing the relations in human-readable format. 
# Returns rec(inv, braid, zipper, untwisting) (GAP record data type).
#
# The relations are grouped so that we can reduce the high-degree zipper/untwist
# relations modulo the invertibility ones before handing them to SGrobnerTrace,
# and so that the display can label them.
build_relations := function(one_, xs_, ys_, q_, t_, s_, untw)
  local n_, inv, braid, zipper, untwisting, S, i, j, k,
        yProd, xProd, yProd_after, xProd_after, yProd_full, xProd_full,
        xProd_reverse, yProd_reverse;
  n_ := Length(xs_) + 1;
  inv := [];          # invertibility
  braid := [];        # braid + far-commutativity (no x*y pairs)
  zipper := [];       # zipper relations (contains x*y pairs)
  untwisting := [];   # untwisting relations (contains x*y pairs)


  # Add invertibility relations
  for i in [1..n_-1] do
    Add(inv, xs_[i]*ys_[i] - one_);
    Add(inv, ys_[i]*xs_[i] - one_);
  od;


  # Add far commutativity and braid relations
  # We include the y-braid and mixed versions explicitly: all derivable from the
  # x-braid relations + invertibility, but giving them upfront saves SGrobner from
  # rediscovering them via S-polynomial reductions.
  for i in [1..n_-2] do
    for j in [i+1..n_-1] do
      if j - i >= 2 then
        Add(braid, xs_[i]*xs_[j] - xs_[j]*xs_[i]);
        Add(braid, ys_[i]*ys_[j] - ys_[j]*ys_[i]);
        Add(braid, ys_[i]*xs_[j] - xs_[j]*ys_[i]);
        Add(braid, xs_[i]*ys_[j] - ys_[j]*xs_[i]);
      else  # j = i+1
        Add(braid, xs_[i]*xs_[j]*xs_[i] - xs_[j]*xs_[i]*xs_[j]);
        Add(braid, ys_[i]*ys_[j]*ys_[i] - ys_[j]*ys_[i]*ys_[j]);
      fi;
    od;
  od;


  # Adds additional length-3 braid identities for adjacent strands.
  # These are derivable from the x-braid + invertibility relations above, but
  # without including them they show up as Grobner basis elements via long
  # S-polynomial chains.
  for i in [1..n_-2] do
    Add(braid, xs_[i]*xs_[i+1]*ys_[i]   - ys_[i+1]*xs_[i]*xs_[i+1]);
    Add(braid, ys_[i]*ys_[i+1]*xs_[i]   - xs_[i+1]*ys_[i]*ys_[i+1]);
    Add(braid, xs_[i]*ys_[i+1]*ys_[i]   - ys_[i+1]*ys_[i]*xs_[i+1]);
    Add(braid, ys_[i]*xs_[i+1]*xs_[i]   - xs_[i+1]*xs_[i]*ys_[i+1]);
  od;


  # Define sliders S[number of strands][anchor]
  # The i-strand sliders have anchor index j ranging over [1..n-i+1].
  S := List([1..n_], i -> []);
  for j in [1..n_-1] do
    S[2][j] := q_*ys_[j] + (1-q_)*one_ - xs_[j];
  od;
  for i in [3..n_] do
    for j in [1..n_-i+1] do
      yProd := ys_[j];
      for k in [j+1..j+i-2] do yProd := yProd * ys_[k]; od;
      xProd := xs_[j];
      for k in [j+1..j+i-2] do xProd := xProd * xs_[k]; od;
      S[i][j] := (q_^(i-1)*yProd - xProd) * S[i-1][j];
    od;
  od;


  # Add zipper relations
  # Note that zipper and untwisting relations for offset sliders (S[i][j],j > 1)
  # can be derived from the S[i][1] ones via conjugation, but again, including
  # them saves from rediscovery via long reductions.
  # Zipper relations for 2-strand sliders
  for j in [1..n_-2] do
    Add(zipper,
        (q_*ys_[j+1] + (1-q_)*one_ - xs_[j+1]
         - q_*ys_[j]*ys_[j+1] + xs_[j]*xs_[j+1]) * S[2][j]);
  od;

  # Zipper relations for 3 or more strand sliders
  # Since the zipper relation only applies if there is free strand to the right of
  # a slider, j ranges in [1..n-i]
  for i in [3..n_-1] do
    for j in [1..n_-i] do
      yProd_after := ys_[j+1];
      for k in [j+2..j+i-1] do yProd_after := yProd_after * ys_[k]; od;
      yProd_full := ys_[j] * yProd_after;
      xProd_after := xs_[j+1];
      for k in [j+2..j+i-1] do xProd_after := xProd_after * xs_[k]; od;
      xProd_full := xs_[j] * xProd_after;
      Add(zipper,
          (q_^(i-1)*yProd_after - xProd_after
           - q_^(i-1)*yProd_full + xProd_full) * S[i][j]);
    od;
  od;


  # Add untwisting relations -- behavior controlled by untw.
  if untw <> "none" then

    # Untwisting relations for 2-strand sliders. "silly_question" weakens the
    # eigenvalue to +/- t (degree 2 in xs[j]); the others fix it at t.
    for j in [1..n_-1] do
      if untw = "silly_question" then
        Add(untwisting, (xs_[j]^2 - t_^2*one_) * S[2][j]);
      else
        Add(untwisting, (xs_[j] - t_*one_) * S[2][j]);
      fi;
    od;

    # Untwisting relations for 3 or more strand sliders.
    # An i-strand slider anchored at strand index j requires strands
    # j, j+1, ..., j+i-2, so j ranges [1..n-i+1].
    # "indeterminate" uses s^(i-1); the others use q^(i-1).
    for i in [3..n_] do
      for j in [1..n_-i+1] do
        xProd_reverse := xs_[j+i-2];
        yProd_reverse := ys_[j+i-2];
        for k in [2..i-1] do
          xProd_reverse := xProd_reverse * xs_[j+i-1-k];
          yProd_reverse := yProd_reverse * ys_[j+i-1-k];
        od;
        if untw = "indeterminate" then
          Add(untwisting, (xProd_reverse - s_^(i-1)*yProd_reverse) * S[i][j]);
        else  # bigelow or silly question
          Add(untwisting, (xProd_reverse - q_^(i-1)*yProd_reverse) * S[i][j]);
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
Print("=== zipper.gap run ===\n");
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
Print("  Variable conventions:\n");
Print("    x_i  =  sigma_i           (positive crossing on strand i)\n");
Print("    y_i  =  sigma_i^{-1}      (negative crossing; y_i = x_i^{-1})\n");
Print("    S[i][j]  =  X_{i,j}       (i-strand slider anchored at strand j)\n");
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
sym_xs := sym_gens{[2..n]};
sym_ys := sym_gens{[n+1..2*n-1]};
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
# x_i / y_i names with run-length compression for repeated factors.
pretty_monomial_str := function(mon)
  local parts, idx, count, name, i;
  if Length(mon) = 0 then return ""; fi;
  parts := [];
  i := 1;
  while i <= Length(mon) do
    count := 1;
    while i + count <= Length(mon) and mon[i + count] = mon[i] do
      count := count + 1;
    od;
    idx := mon[i];
    if idx <= n - 1 then
      name := Concatenation("x", String(idx));
    else
      name := Concatenation("y", String(idx - (n - 1)));
    fi;
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

Print("Grobner basis size: ", Length(B), "\n\n");

# Save the result to a parameter-tagged file.
outFile := Concatenation("logs_and_traces/gap_traces/grobner_zipper_",
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
nGen := 2*(n-1);   # number of generators of the free algebra (xs and ys)
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