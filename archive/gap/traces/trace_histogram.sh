#!/usr/bin/env bash
# Summarize a saved SGrobnerTrace output produced by zipper.gap.
#
# Usage:  trace_histogram.sh path/to/grobner_trace_zipper-...gap
#
# Streams over the file once and reports:
#   * file size and basis size
#   * trace-tuple usage histogram by input-relation index, with each I[k]
#     labeled by what it actually is (matched to zipper.gap's input order)
#   * per-basis-element table of pol-term count and trace-tuple count
#   * aggregate left/right word-length distribution across all trace tuples
#
# The trace tuple format saved by GBNP is
#     [ <leftWord>, <inputIdx>, <rightWord>, <coefficient> ]
# so a `], <number>, [` pattern uniquely picks out the inputIdx field.

set -euo pipefail


# === Argument parsing ===

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <grobner_trace_zipper-*.gap>" >&2
  exit 2
fi

trace_file="$1"
if [[ ! -f "$trace_file" ]]; then
  echo "error: file not found: $trace_file" >&2
  exit 1
fi


# === File metadata ===

echo "=== File ==="
ls -lh "$trace_file" \
  | awk '{printf "  path: %s\n  size: %s\n", $9, $5}'

# Extract n (number of strands) from the saved file's `n := <int>;` line.
# Determines the input layout and what each I[k] is.
n_strands=$(grep -m1 -E '^n := [0-9]+;' "$trace_file" \
            | grep -oE '[0-9]+' | head -1)
if [[ -z "${n_strands:-}" ]]; then
  echo "  warning: could not extract n; labels fall back to bare indices" >&2
  n_strands=0
fi
echo "  n = $n_strands strands"


# === Build input-relation labels matching zipper.gap's input order ===
#
# The order is:
#   1. Invertibility:                x_i*y_i = 1, y_i*x_i = 1   (i = 1..n-1)
#   2. Far-commutativity and braid:  one (i, j) double loop with j > i
#   3. Length-3 braid identities:    one i = 1..n-2 loop, 4 each
#   4. Per-pair zippers:             one j = 1..n-2 loop
#   5. i-strand zippers (i >= 3):    one (i, j) loop, j in [1..n-i]
#   6. Per-strand unzippings:        one j = 1..n-1 loop
#   7. Bigelow untwistings (i >= 3): one (i, j) loop, j in [1..n-i+1]
#
# Labels follow the slider-centric naming from zipper.gap: zip(i;j) and
# unzip(i;j) where i is the slider's level (= i in S[i][j]) and j is the
# anchor strand.

labels=()

if (( n_strands >= 2 )); then
  # 1. invertibility
  for ((i=1; i<n_strands; i++)); do
    labels+=("x_${i}*y_${i}=1")
    labels+=("y_${i}*x_${i}=1")
  done

  # 2. far-commutativity + braid (one (i, j) loop)
  for ((i=1; i<n_strands-1; i++)); do
    for ((j=i+1; j<n_strands; j++)); do
      if (( j - i >= 2 )); then
        labels+=("[x_${i},x_${j}]")
        labels+=("[y_${i},y_${j}]")
        labels+=("[y_${i},x_${j}]")
        labels+=("[x_${i},y_${j}]")
      else
        labels+=("x-braid(${i},${j})")
        labels+=("y-braid(${i},${j})")
      fi
    done
  done

  # 3. length-3 braid identities (mixed-sign).
  # Labels are the leading-monomial pattern (each is unique).
  for ((i=1; i<n_strands-1; i++)); do
    j=$((i+1))
    labels+=("xxy(${i},${j})")    # x_i*x_j*y_i = y_j*x_i*x_j
    labels+=("yyx(${i},${j})")    # y_i*y_j*x_i = x_j*y_i*y_j
    labels+=("xyy(${i},${j})")    # x_i*y_j*y_i = y_j*y_i*x_j
    labels+=("yxx(${i},${j})")    # y_i*x_j*x_i = x_j*x_i*y_j
  done

  # 4. per-pair zippers, using S[2][j]
  for ((j=1; j<=n_strands-2; j++)); do
    labels+=("zip(2;${j})")
  done

  # 5. i-strand zippers (i >= 3), using S[i][j]
  for ((i=3; i<=n_strands-1; i++)); do
    for ((j=1; j<=n_strands-i; j++)); do
      labels+=("zip(${i};${j})")
    done
  done

  # 6. per-strand unzippings, using S[2][j]
  for ((j=1; j<=n_strands-1; j++)); do
    labels+=("unzip(2;${j})")
  done

  # 7. Bigelow untwistings (i >= 3), using S[i][j]
  for ((i=3; i<=n_strands; i++)); do
    for ((j=1; j<=n_strands-i+1; j++)); do
      labels+=("unzip(${i};${j})")
    done
  done
fi


# === Basis size ===

echo
echo "=== Basis elements ==="
basis_size=$(grep -c 'pol :=' "$trace_file" || true)
echo "  basis size: $basis_size"


# === Trace-tuple usage histogram by input-relation index ===

echo
echo "=== Trace-tuple usage by input-relation index ==="

# Two scratch files: one for the count-by-idx output, one for the labels.
counts_file=$(mktemp)
labels_file=$(mktemp)
trap 'rm -f "$counts_file" "$labels_file"' EXIT

grep -oE '\], [0-9]+, \[' "$trace_file" \
  | grep -oE '[0-9]+' \
  | sort -n | uniq -c \
  > "$counts_file"

total_tuples=$(awk '{s+=$1} END {print s+0}' "$counts_file")
echo "  total trace tuples: $total_tuples"
if (( basis_size > 0 && total_tuples > 0 )); then
  awk -v n="$basis_size" -v t="$total_tuples" \
    'BEGIN {printf "  mean tuples per basis element: %.0f\n", t/n}'
fi

# Write labels (one per line, indexed by line number = I[k]).
: > "$labels_file"
for label in "${labels[@]}"; do
  printf '%s\n' "$label" >> "$labels_file"
done

echo
printf "  %-9s %-9s %-7s %s\n" count share "I[idx]" relation
printf "  %-9s %-9s %-7s %s\n" -------- -------- ------- --------
awk -v t="$total_tuples" -v lf="$labels_file" '
  BEGIN {
    while ((getline line < lf) > 0) labels[++n_labels] = line
    close(lf)
  }
  {
    pct = (t > 0) ? (100.0 * $1 / t) : 0
    label = ($2 in labels) ? labels[$2 + 0] : "?"
    printf "  %9d %8.2f%% I[%-4s] %s\n", $1, pct, $2, label
  }
' "$counts_file"


# === Per-basis-element analysis and word-length distribution ===

echo
echo "=== Per-basis-element analysis ==="
echo "    (pol terms, trace length per element + word-length histogram)"

gawk '
  # State machine over the GBT section of the file.  Each basis element is a
  # rec(pol := <NP poly>, trace := <list of tuples>).  We track which section
  # we are in and accumulate counters; on transitioning out, we record per-
  # element values into pol_size[] and trace_size[].

  /pol := / {
    if (state == "trace") {
      element_idx++
      pol_size[element_idx] = saved_pol_terms
      trace_size[element_idx] = trace_tuples
    }
    state = "pol"
    pol_terms = 0
    trace_tuples = 0
  }

  /trace := / {
    saved_pol_terms = pol_terms
    state = "trace"
    trace_tuples = 0
  }

  {
    if (state == "pol") {
      # Each monomial has exactly one ZmodpZObj coefficient.
      s = $0
      while (match(s, /ZmodpZObj\(/)) {
        pol_terms++
        s = substr(s, RSTART + RLENGTH)
      }
    } else if (state == "trace") {
      # Count trace tuples via the unique idx-field pattern.
      s = $0
      while (match(s, /\], [0-9]+, \[/)) {
        trace_tuples++
        s = substr(s, RSTART + RLENGTH)
      }

      # Extract left/right word lengths per tuple.
      # Pattern: [ [ <leftContents> ], <idx>, [ <rightContents> ], ZmodpZObj
      # The contents are integers separated by commas (or whitespace if empty).
      s = $0
      while (match(s,
          /\[ \[([^][]*)\], [0-9]+, \[([^][]*)\], ZmodpZObj/, captured)) {
        left  = captured[1]
        right = captured[2]
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", left)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", right)
        if (left  == "") { l_len = 0 } \
                    else { tmp = left;  gsub(/[^,]/, "", tmp); \
                           l_len = length(tmp) + 1 }
        if (right == "") { r_len = 0 } \
                    else { tmp = right; gsub(/[^,]/, "", tmp); \
                           r_len = length(tmp) + 1 }
        left_word_dist[l_len]++
        right_word_dist[r_len]++
        s = substr(s, RSTART + RLENGTH)
      }
    }
  }

  END {
    if (state == "trace") {
      element_idx++
      pol_size[element_idx] = saved_pol_terms
      trace_size[element_idx] = trace_tuples
    }

    # Per-element table.
    printf "\n  %-6s %10s  %14s\n", "elem", "pol_terms", "trace_tuples"
    printf "  %-6s %10s  %14s\n",   "----", "---------", "------------"
    for (i = 1; i <= element_idx; i++) {
      printf "  B[%-3d] %10d  %14d\n", i, pol_size[i], trace_size[i]
    }

    # Summary stats over per-element trace and pol sizes.
    if (element_idx > 0) {
      t_min = trace_size[1]; t_max = t_min; t_sum = 0
      p_min = pol_size[1];   p_max = p_min; p_sum = 0
      for (i = 1; i <= element_idx; i++) {
        if (trace_size[i] < t_min) t_min = trace_size[i]
        if (trace_size[i] > t_max) t_max = trace_size[i]
        t_sum += trace_size[i]
        if (pol_size[i] < p_min) p_min = pol_size[i]
        if (pol_size[i] > p_max) p_max = pol_size[i]
        p_sum += pol_size[i]
      }
      printf "\n  trace length:  min=%d  mean=%d  max=%d\n",
             t_min, t_sum/element_idx, t_max
      printf "  pol terms:     min=%d  mean=%d  max=%d\n",
             p_min, p_sum/element_idx, p_max
    }

    # Word-length distribution across all trace tuples.
    print ""
    print "  word len | left tuples |  right tuples"
    print "  -------- | ----------- |  ------------"
    max_word_len = 0
    for (l in left_word_dist)  if (l+0 > max_word_len) max_word_len = l+0
    for (l in right_word_dist) if (l+0 > max_word_len) max_word_len = l+0
    for (l = 0; l <= max_word_len; l++) {
      if ((l in left_word_dist) || (l in right_word_dist)) {
        l_count = (l in left_word_dist)  ? left_word_dist[l]  : 0
        r_count = (l in right_word_dist) ? right_word_dist[l] : 0
        printf "  %8d | %11d |  %12d\n", l, l_count, r_count
      }
    }
  }
' "$trace_file"
