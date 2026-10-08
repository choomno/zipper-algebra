LoadPackage("gbnp");
SetInfoLevel(InfoGBNP, 1);
SetRecursionTrapInterval(50000);


###############################################################################
# USER PARAMETERS
###############################################################################
#
# Edit the values in this block, then run
#   gap -o <workspace> -K <workspace> trace_replay.gap
# from the braid-GAP/ directory.
#
# This script extracts a witness for the target element X_n = S[n][1] from a
# saved traced Grobner basis (the "anchor" run), then re-derives the witness
# coefficients at a list of additional parameter points by linear solve. The
# verification at each sample point is "rebuild = target", confirming the
# linear solver found a valid decomposition. (Whether the solver picks the
# SAME branch across samples -- needed for CRT/interpolation -- is checked by
# whether the recovered c-vectors interpolate consistently downstream.)


# Path to the anchor saved-traced-basis file. Loading it brings into scope:
#   n, pStr, p, R, q, t, I, GBT
# Override via ANCHOR_FILE environment variable, or edit the default below.
if IsBound(GAPInfo.SystemEnvironment.ANCHOR_FILE)
   and GAPInfo.SystemEnvironment.ANCHOR_FILE <> "" then
  basisFile := GAPInfo.SystemEnvironment.ANCHOR_FILE;
else
  basisFile := "logs_and_traces/gap_traces/grobner_zipper_bigelow-trace-3-(2^31-1)-5-7.gap";
fi;

# Sample parameter points at which to replay the trace. Each entry is
#   [p_val, q_int, t_int]
# with p_val a prime fitting comfortably in 31 bits, and q_int, t_int small
# nonzero integers (preferably distinct small primes, avoiding 0 and 1).
# The first entry should be the anchor parameters themselves (a sanity check
# that the linear solve recovers a valid witness at the original point).
samples := [
  [ 2^31 - 1,   5,  7 ],   # anchor (matches the loaded basis)
  [ 2^31 - 1,   3, 11 ],   # different (q,t), same prime
  [ 2^31 - 19,  5,  7 ],   # different prime, same (q,t)
  [ 2^31 - 19,  3, 11 ],   # different prime AND different (q,t)
  [ 2^31 - 61,  5,  7 ],   # another prime, same (q,t)
];

# What to replay at each sample.
#   "single" - just the trace of the target X_n = S[n][1] (the conjectured
#              identity we're trying to verify). Cheaper.
#   "full"   - the trace of every basis element in GBT, lifting each B[k].pol
#              from the anchor field to the sample field by preserving its
#              integer-representative coefficients. Useful if you want to
#              recover symbolic coefficients for the entire basis (not just
#              X_n) via downstream CRT/interpolation. Roughly N times more
#              work where N = Length(GBT).
replay_mode := "single";

###############################################################################
# END USER PARAMETERS
###############################################################################


# === Load anchor traced basis ===
Read(basisFile);
Print("Loaded basis from ", basisFile, "\n");
Print("  n=", n, ", |GBT|=", Length(GBT),
      ", GF(", pStr, "), q=", Int(q), ", t=", Int(t), "\n");

# Log everything that follows to a parameter-tagged logfile.
LogTo(Concatenation("logs_and_traces/gap_traces/logfile_replay-",
                    String(n), "-(", pStr, ")-",
                    String(Int(q)), "-", String(Int(t)), ".txt"));


# === Helper: build the input list and target X_n at chosen parameters ===
# Mirrors zipper.gap's relation construction. Returns rec(I, target, R).
# Take p_val (a prime) plus q_int, t_int (lifted to GF(p_val) via *One).
BuildInputsAndTarget := function(p_val, q_int, t_int)
  local R_loc, q_loc, t_loc, gens, F_loc, gen_list, one_loc, xs_loc, ys_loc,
        inv_rels, basic_rels, complex_rels, S_loc, i, j, k,
        yProd, xProd, yProd_after, yProd_full, xProd_after, xProd_full,
        xProd_reverse, yProd_reverse,
        inv_GB, complex_np, complex_simp, I_loc, target_loc;

  R_loc := GF(p_val);
  q_loc := q_int * One(R_loc);
  t_loc := t_int * One(R_loc);

  gens := Concatenation(
    List([1..n-1], a -> Concatenation("x", String(a))),
    List([1..n-1], a -> Concatenation("y", String(a))));
  F_loc := FreeAssociativeAlgebraWithOne(R_loc, gens);
  gen_list := GeneratorsOfAlgebra(F_loc);
  one_loc := gen_list[1];
  xs_loc := gen_list{[2..n]};
  ys_loc := gen_list{[n+1..2*n-1]};

  # invertibility
  inv_rels := [];
  for i in [1..n-1] do
    Add(inv_rels, xs_loc[i]*ys_loc[i] - one_loc);
    Add(inv_rels, ys_loc[i]*xs_loc[i] - one_loc);
  od;

  # braid + far-comm + length-3 mixed-sign identities
  basic_rels := [];
  for i in [1..n-2] do
    for j in [i+1..n-1] do
      if j - i >= 2 then
        Add(basic_rels, xs_loc[i]*xs_loc[j] - xs_loc[j]*xs_loc[i]);
        Add(basic_rels, ys_loc[i]*ys_loc[j] - ys_loc[j]*ys_loc[i]);
        Add(basic_rels, ys_loc[i]*xs_loc[j] - xs_loc[j]*ys_loc[i]);
        Add(basic_rels, xs_loc[i]*ys_loc[j] - ys_loc[j]*xs_loc[i]);
      else
        Add(basic_rels, xs_loc[i]*xs_loc[j]*xs_loc[i]
                        - xs_loc[j]*xs_loc[i]*xs_loc[j]);
        Add(basic_rels, ys_loc[i]*ys_loc[j]*ys_loc[i]
                        - ys_loc[j]*ys_loc[i]*ys_loc[j]);
      fi;
    od;
  od;
  for i in [1..n-2] do
    Add(basic_rels, xs_loc[i]*xs_loc[i+1]*ys_loc[i]
                    - ys_loc[i+1]*xs_loc[i]*xs_loc[i+1]);
    Add(basic_rels, ys_loc[i]*ys_loc[i+1]*xs_loc[i]
                    - xs_loc[i+1]*ys_loc[i]*ys_loc[i+1]);
    Add(basic_rels, xs_loc[i]*ys_loc[i+1]*ys_loc[i]
                    - ys_loc[i+1]*ys_loc[i]*xs_loc[i+1]);
    Add(basic_rels, ys_loc[i]*xs_loc[i+1]*xs_loc[i]
                    - xs_loc[i+1]*xs_loc[i]*ys_loc[i+1]);
  od;

  # sliders S[level][strand]
  S_loc := List([1..n], a -> []);
  for j in [1..n-1] do
    S_loc[2][j] := q_loc*ys_loc[j] + (1-q_loc)*one_loc - xs_loc[j];
  od;
  for i in [3..n] do
    for j in [1..n-i+1] do
      yProd := ys_loc[j];
      for k in [j+1..j+i-2] do yProd := yProd * ys_loc[k]; od;
      xProd := xs_loc[j];
      for k in [j+1..j+i-2] do xProd := xProd * xs_loc[k]; od;
      S_loc[i][j] := (q_loc^(i-1)*yProd - xProd) * S_loc[i-1][j];
    od;
  od;

  # complex relations (zippers + untwistings)
  complex_rels := [];
  for j in [1..n-2] do
    Add(complex_rels,
        (q_loc*ys_loc[j+1] + (1-q_loc)*one_loc - xs_loc[j+1]
         - q_loc*ys_loc[j]*ys_loc[j+1]
         + xs_loc[j]*xs_loc[j+1]) * S_loc[2][j]);
  od;
  for i in [3..n-1] do
    for j in [1..n-i] do
      yProd_after := ys_loc[j+1];
      for k in [j+2..j+i-1] do yProd_after := yProd_after * ys_loc[k]; od;
      yProd_full := ys_loc[j] * yProd_after;
      xProd_after := xs_loc[j+1];
      for k in [j+2..j+i-1] do xProd_after := xProd_after * xs_loc[k]; od;
      xProd_full := xs_loc[j] * xProd_after;
      Add(complex_rels,
          (q_loc^(i-1)*yProd_after - xProd_after
           - q_loc^(i-1)*yProd_full + xProd_full) * S_loc[i][j]);
    od;
  od;
  for j in [1..n-1] do
    Add(complex_rels, (xs_loc[j] - t_loc*one_loc) * S_loc[2][j]);
  od;
  for i in [3..n] do
    for j in [1..n-i+1] do
      xProd_reverse := xs_loc[j+i-2];
      yProd_reverse := ys_loc[j+i-2];
      for k in [2..i-1] do
        xProd_reverse := xProd_reverse * xs_loc[j+i-1-k];
        yProd_reverse := yProd_reverse * ys_loc[j+i-1-k];
      od;
      Add(complex_rels,
          (xProd_reverse - q_loc^(i-1)*yProd_reverse) * S_loc[i][j]);
    od;
  od;

  # pre-reduce complex relations modulo invertibility
  inv_GB := GP2NPList(inv_rels);
  complex_np := GP2NPList(complex_rels);
  complex_simp := List(complex_np, p_pol -> StrongNormalFormNP(p_pol, inv_GB));
  complex_simp := Filtered(complex_simp, p_pol -> p_pol <> []);

  I_loc := Concatenation(inv_GB,
                         GP2NPList(basic_rels),
                         complex_simp);
  target_loc := GP2NPList([S_loc[n][1]])[1];

  return rec(I := I_loc, target := target_loc, R := R_loc);
end;


# === Helpers for building/rebuilding NP polynomials ===

# Scalar-multiply an NP polynomial.
MulScalarNP := function(p_pol, c)
  if p_pol = [] then return []; fi;
  return [p_pol[1], List(p_pol[2], x -> c*x)];
end;

# Lift an NP polynomial's coefficients from one field to another by preserving
# the integer representative of each ZmodpZObj. Used to "transport" anchor
# basis-element pol fields to a fresh sample field as targets for replay.
LiftPolToField := function(p_pol, R_new)
  if p_pol = [] then return []; fi;
  return [p_pol[1], List(p_pol[2], c -> Int(c) * One(R_new))];
end;

# Rebuild  sum_k c_vec[k] * left_k * I[i_k] * right_k  as an NP polynomial.
RebuildFromCoeffs := function(trace_X_, I_in, c_vec, R_in)
  local T_, k, sum_pol, term, oneR;
  T_ := Length(trace_X_);
  oneR := One(R_in);
  sum_pol := [];
  for k in [1..T_] do
    if not IsZero(c_vec[k]) then
      term := BimulNP(trace_X_[k][1], I_in[trace_X_[k][2]], trace_X_[k][3]);
      term := MulScalarNP(term, c_vec[k]);
      sum_pol := AddNP(sum_pol, term, oneR, oneR);
    fi;
  od;
  return sum_pol;
end;


# === Replay function ===
# Given a trace structure, new inputs, new target, and the new ground field,
# solve the linear system for new coefficients c'_k satisfying
#     target = sum_k c'_k * left_k * I_new[i_k] * right_k.
ReplayTrace := function(trace_X_, I_new, target_new, R_new)
  local T_, building, all_mons, k, idx, lookup, A, b, sol;

  T_ := Length(trace_X_);

  # Step 1: building blocks B_k = left_k * I_new[i_k] * right_k.
  building := List(trace_X_, tup -> BimulNP(tup[1], I_new[tup[2]], tup[3]));

  # Step 2: collect all distinct monomials appearing in target plus the
  # building blocks; this is the row index of the linear system.
  all_mons := [];
  if target_new <> [] then UniteSet(all_mons, target_new[1]); fi;
  for k in [1..T_] do
    if building[k] <> [] then UniteSet(all_mons, building[k][1]); fi;
  od;
  lookup := function(mon) return PositionSorted(all_mons, mon); end;

  # Step 3: build M x T matrix A and length-M vector b over R_new such that
  # A_{m,k} = (coeff of monomial m in B_k)  and  b_m = (coeff of m in target).
  A := NullMat(Length(all_mons), T_, R_new);
  b := List([1..Length(all_mons)], m -> Zero(R_new));
  if target_new <> [] then
    for k in [1..Length(target_new[1])] do
      idx := lookup(target_new[1][k]);
      b[idx] := target_new[2][k];
    od;
  fi;
  for k in [1..T_] do
    if building[k] <> [] then
      for idx in [1..Length(building[k][1])] do
        A[lookup(building[k][1][idx])][k] := building[k][2][idx];
      od;
    fi;
  od;

  # Step 4: solve A * c = b. SolutionMat(M, v) finds row vector x with x*M = v,
  # so we transpose A and treat b as a row vector.
  sol := SolutionMat(TransposedMat(A), b);
  if sol = fail then
    Print("    WARNING: linear system has no solution at this point.\n");
    return fail;
  fi;
  return sol;
end;


# === Extract trace at the anchor ===
# StrongNormalFormTraceDiff returns rec(pol := X - nf(X), trace := ...). When
# X is in the ideal (nf = 0), pol = X and trace expresses X over the input I.
anchor := BuildInputsAndTarget(p, Int(q), Int(t));
Print("Anchor inputs match the saved I: ", anchor.I = I, "\n");

target_orig_np := anchor.target;
B_pol := List(GBT, r -> r.pol);
nf := StrongNormalFormNP(target_orig_np, B_pol);
Print("Target X_n is in the ideal? ", nf = [], "\n");


# === Build the list of (label, trace, target_at_anchor) triples to replay ===
# In "single" mode there's just one entry: the X_n target.
# In "full"   mode there's one per basis element of GBT.
replay_items := [];
if replay_mode = "single" then
  res := StrongNormalFormTraceDiff(target_orig_np, GBT);
  Add(replay_items, rec(label := "X_n",
                        trace := res.trace,
                        anchor_pol := target_orig_np));
elif replay_mode = "full" then
  for k in [1..Length(GBT)] do
    Add(replay_items,
        rec(label := Concatenation("B[", String(k), "]"),
            trace := GBT[k].trace,
            anchor_pol := GBT[k].pol));
  od;
else
  Error("unknown replay_mode: \"", replay_mode, "\"");
fi;
Print("Replay mode: ", replay_mode,
      "  (", Length(replay_items), " trace(s) per sample)\n");
total_trace_tuples := Sum(replay_items, item -> Length(item.trace));
Print("Total trace tuples to replay per sample: ", total_trace_tuples, "\n\n");


# === Replay at every sample point ===
# Each result has the form
#   rec(sample, ok_count, total_count, items := [ rec(label, ok, c), ... ])
# where ok_count = number of items where rebuild = target.
results := [];
for sample in samples do
  Print("--- Sample p=", sample[1], ", q=", sample[2], ", t=", sample[3],
        " ---\n");
  fresh := BuildInputsAndTarget(sample[1], sample[2], sample[3]);
  items := [];
  ok_count := 0;
  for item in replay_items do
    # Lift the anchor target to the sample field. For X_n the lift is the
    # freshly-built target (computed via the slider construction at sample
    # parameters); for arbitrary basis elements we numerically lift their
    # anchor pol via LiftPolToField, preserving the integer representatives.
    if item.label = "X_n" then
      target_sample := fresh.target;
    else
      target_sample := LiftPolToField(item.anchor_pol, fresh.R);
    fi;
    c_vec := ReplayTrace(item.trace, fresh.I, target_sample, fresh.R);
    if c_vec = fail then
      Add(items, rec(label := item.label, ok := false, c := fail));
    else
      rebuilt := RebuildFromCoeffs(item.trace, fresh.I, c_vec, fresh.R);
      matches := (rebuilt = target_sample);
      if matches then ok_count := ok_count + 1; fi;
      Add(items, rec(label := item.label, ok := matches, c := c_vec));
    fi;
  od;
  Print("    ", ok_count, "/", Length(items), " replays verified\n\n");
  Add(results, rec(sample := sample,
                   ok_count := ok_count,
                   total_count := Length(items),
                   items := items));
od;

Print("=== Summary ===\n");
for r in results do
  Print("  p=", r.sample[1], " q=", r.sample[2], " t=", r.sample[3],
        "  ", r.ok_count, "/", r.total_count, " ok\n");
od;


# === Save GAP-loadable data file ===
# Contents:
#   anchor metadata (n, pStr_anchor, q_int_anchor, t_int_anchor)
#   replay_mode             - "single" or "full"
#   replay_items            - list of rec(label, trace, anchor_pol)
#                             one entry per replayed trace
#   replay_results          - list of rec(sample, items := [rec(label, ok, c)])
#                             one entry per sample
# Loading this file plus rebuilding I via BuildInputsAndTarget is sufficient
# for downstream CRT/interpolation work.
dataFile := Concatenation("logs_and_traces/gap_traces/replay_data-", replay_mode, "-",
                          String(n), "-(", pStr, ")-",
                          String(Int(q)), "-", String(Int(t)), ".gap");
PrintTo(dataFile,
  "# Trace-replay data produced by trace_replay.gap.\n",
  "# Anchor: n strands, prime pStr_anchor, parameters",
       " (q_int_anchor, t_int_anchor).\n",
  "n := ", n, ";\n",
  "pStr_anchor := \"", pStr, "\";\n",
  "q_int_anchor := ", Int(q), ";\n",
  "t_int_anchor := ", Int(t), ";\n",
  "replay_mode := \"", replay_mode, "\";\n",
  "replay_items := [\n");
for i in [1..Length(replay_items)] do
  item := replay_items[i];
  AppendTo(dataFile,
    "  rec( label := \"", item.label,
    "\",\n       trace := ", item.trace,
    ",\n       anchor_pol := ", item.anchor_pol, " )");
  if i < Length(replay_items) then AppendTo(dataFile, ",\n");
                              else AppendTo(dataFile, "\n");
  fi;
od;
AppendTo(dataFile, "];\nreplay_results := [\n");
for i in [1..Length(results)] do
  r := results[i];
  AppendTo(dataFile,
    "  rec( sample := ", r.sample,
    ",\n       ok_count := ", r.ok_count,
    ", total_count := ", r.total_count,
    ",\n       items := [\n");
  for j in [1..Length(r.items)] do
    it := r.items[j];
    AppendTo(dataFile,
      "         rec( label := \"", it.label,
      "\", ok := ", it.ok,
      ", c := ", it.c, " )");
    if j < Length(r.items) then AppendTo(dataFile, ",\n");
                           else AppendTo(dataFile, "\n");
    fi;
  od;
  AppendTo(dataFile, "       ] )");
  if i < Length(results) then AppendTo(dataFile, ",\n");
                         else AppendTo(dataFile, "\n");
  fi;
od;
AppendTo(dataFile, "];\n");
Print("\nWrote loadable data to ", dataFile, "\n");

LogTo();   # close the logfile

# Notes for downstream use:
#   * If "ok" is true at every sample, the linear solver found valid witness
#     coefficients at each point. This does NOT guarantee that the recovered
#     c-vectors lie on the same solution branch (the system is generally
#     underdetermined). To check: try interpolating each c_k(q,t) via Stage-2
#     of the CRT/interpolation pipeline and see whether the interpolant
#     reproduces the values at fresh held-out sample points.
#   * The list `results` holds one record per sample with fields {sample, c,
#     ok}. Subsequent CRT and bivariate-rational reconstruction code should
#     consume this.
