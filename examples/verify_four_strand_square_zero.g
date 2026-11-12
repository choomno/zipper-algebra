# Independent exact GAP/GBNP verification of the four-strand square-zero case.
# Run from the repository root:
#   gap -A -q -r --quitonbreak examples/verify_four_strand_square_zero.g
# Set SquareZeroProofFile with GAP's -c option to select another proof.txt.
# Generated data are parsed as arithmetic and words; they are never executed.
# GBNP's IsGrobnerBasis checks the complete noncommutative ambiguity criterion.

SquareZero := rec();;
# T4 is the 31st defining relation. Every derivation must use only the first
# 30 inputs, so reducing all 31 against the proved basis also certifies that
# adjoining T4 does not change the algebra.
SquareZero.auxiliaryInputCount := 30;;
SquareZero.fail := function(message)
    PrintTo("*errout*", "FAIL: ", message, "\n");
    QUIT_GAP(1);
end;;
if LoadPackage("gbnp") <> true then SquareZero.fail("GBNP is required"); fi;

SquareZero.fields := function(text)
    local answer, current, c;
    answer := []; current := "";
    for c in text do
        if c = '|' then Add(answer, current); current := "";
        else Add(current, c); fi;
    od;
    Add(answer, current);
    return answer;
end;;

SquareZero.integer := function(text)
    if Length(text) = 0 or not ForAll(text, c -> c in "0123456789") then
        SquareZero.fail("expected a nonnegative integer");
    fi;
    return Int(text);
end;;

SquareZero.word := function(text)
    local pieces, word;
    if text = "" then return []; fi;
    if text[Length(text)] <> ',' then
        SquareZero.fail("word must end with a comma");
    fi;
    pieces := SplitString(text, ",");
    word := List(pieces, SquareZero.integer);
    if not ForAll(word, i -> i >= 1 and i <= 6) then
        SquareZero.fail("word generator outside 1..6");
    fi;
    return word;
end;;

# Recursive descent: sums, products, unary signs, powers, parentheses, q,t.
# Every scalar belongs to the same cancelling rational-function field.
SquareZero.scalar := function(text, field)
    local chars, pos, vars, expression, product, unary, power, primary, value;
    chars := Filtered(text, c -> not c in " \t\r\n");
    if Length(chars) = 0 or Length(chars) > 100000 then
        SquareZero.fail("empty or excessively long coefficient");
    fi;
    pos := 1;
    vars := IndeterminatesOfFunctionField(field);
    primary := function()
        local start, answer;
        if pos > Length(chars) then SquareZero.fail("incomplete coefficient"); fi;
        if chars[pos] = '(' then
            pos := pos + 1; answer := expression();
            if pos > Length(chars) or chars[pos] <> ')' then
                SquareZero.fail("unbalanced coefficient parentheses");
            fi;
            pos := pos + 1; return answer;
        elif chars[pos] in "qt" then
            answer := vars[Position("qt", chars[pos])];
            pos := pos + 1; return answer;
        fi;
        start := pos;
        while pos <= Length(chars) and chars[pos] in "0123456789" do
            pos := pos + 1;
        od;
        if start = pos then SquareZero.fail("invalid coefficient token"); fi;
        return One(field) * SquareZero.integer(chars{[start..pos-1]});
    end;
    power := function()
        local answer, exponent, start, sign;
        answer := primary();
        if pos <= Length(chars) and chars[pos] = '^' then
            pos := pos + 1; sign := 1;
            if pos <= Length(chars) and chars[pos] = '-' then
                sign := -1; pos := pos + 1;
            fi;
            start := pos;
            while pos <= Length(chars) and chars[pos] in "0123456789" do
                pos := pos + 1;
            od;
            if start = pos then SquareZero.fail("expected an integer exponent"); fi;
            exponent := SquareZero.integer(chars{[start..pos-1]});
            if exponent > 10000 then SquareZero.fail("exponent exceeds limit"); fi;
            if sign = -1 and IsZero(answer) then SquareZero.fail("division by zero"); fi;
            answer := answer^(sign * exponent);
        fi;
        return answer;
    end;
    unary := function()
        local sign;
        sign := 1;
        while pos <= Length(chars) and chars[pos] in "+-" do
            if chars[pos] = '-' then sign := -sign; fi;
            pos := pos + 1;
        od;
        return sign * power();
    end;
    product := function()
        local answer, op, factor;
        answer := unary();
        while pos <= Length(chars) and chars[pos] in "*/" do
            op := chars[pos]; pos := pos + 1; factor := unary();
            if op = '*' then answer := answer * factor;
            else
                if IsZero(factor) then SquareZero.fail("division by zero"); fi;
                answer := answer / factor;
            fi;
        od;
        return answer;
    end;
    expression := function()
        local answer, op, term;
        answer := product();
        while pos <= Length(chars) and chars[pos] in "+-" do
            op := chars[pos]; pos := pos + 1; term := product();
            if op = '+' then answer := answer + term;
            else answer := answer - term; fi;
        od;
        return answer;
    end;
    value := expression();
    if pos <= Length(chars) then SquareZero.fail("trailing coefficient token"); fi;
    return value;
end;;

SquareZero.read := function(path, field)
    local stream, line, fields, result, counts, ended, row, index, poly, term;
    stream := InputTextFile(path);
    if stream = fail then SquareZero.fail("cannot open proof data"); fi;
    line := ReadLine(stream);
    if line = fail then SquareZero.fail("missing META record"); fi;
    fields := SquareZero.fields(Chomp(line));
    if Length(fields) <> 4 or fields[1] <> "META" then
        SquareZero.fail("missing META record");
    fi;
    counts := List(fields{[2..4]}, SquareZero.integer);
    if counts[1] <> 31 or counts[2] < 1 or counts[2] > counts[3]
        or counts[3] > 10000 then SquareZero.fail("invalid proof counts"); fi;
    result := rec(inputs := List([1..counts[1]], i -> [[], []]),
        nodes := List([1..counts[3]], i -> [[], []]),
        proofs := List([1..counts[3]], i -> []), basisCount := counts[2]);
    ended := false;
    while true do
        line := ReadLine(stream);
        if line = fail then break; fi;
        fields := SquareZero.fields(Chomp(line));
        if ended then SquareZero.fail("records after END"); fi;
        if fields[1] = "END" then
            if Length(fields) <> 4 or List(fields{[2..4]}, SquareZero.integer) <> counts then
                SquareZero.fail("invalid END counts");
            fi;
            ended := true;
        elif fields[1] in ["I", "G"] then
            if Length(fields) <> 6 or fields[5] <> "0" or fields[6] <> "" then
                SquareZero.fail("invalid polynomial record");
            fi;
            index := SquareZero.integer(fields[2]);
            if fields[1] = "I" then row := result.inputs; else row := result.nodes; fi;
            if index < 1 or index > Length(row) then SquareZero.fail("polynomial index outside bounds"); fi;
            poly := row[index];
            Add(poly[1], SquareZero.word(fields[4]));
            Add(poly[2], SquareZero.scalar(fields[3], field));
        elif fields[1] = "D" then
            # Validate raw references even when their scalar is zero.
            if Length(fields) <> 7 then SquareZero.fail("invalid derivation record"); fi;
            index := SquareZero.integer(fields[2]);
            if index < 1 or index > counts[3] then SquareZero.fail("derivation index outside bounds"); fi;
            term := [SquareZero.scalar(fields[3], field), SquareZero.word(fields[4]),
                SquareZero.integer(fields[5]), SquareZero.integer(fields[6]),
                SquareZero.word(fields[7])];
            if term[4] < 1 or (term[3] = 0 and term[4] > SquareZero.auxiliaryInputCount)
                or (term[3] = 1 and term[4] > counts[3]) or not term[3] in [0, 1] then
                SquareZero.fail("derivation reference outside bounds");
            fi;
            Add(result.proofs[index], term);
        else SquareZero.fail("unknown proof record");
        fi;
    od;
    CloseStream(stream);
    if not ended then SquareZero.fail("missing END record"); fi;
    result.inputs := List(result.inputs, CleanNP);
    result.nodes := List(result.nodes, CleanNP);
    if ForAny(Concatenation(result.inputs, result.nodes), p -> p = [[], []])
        or ForAny(result.proofs, p -> Length(p) = 0) then
        SquareZero.fail("missing or zero polynomial or derivation");
    fi;
    return result;
end;;

# Independent construction of the 31 ordered defining relations, X3 and w.
SquareZero.presentation := function(field)
    local z, q, t, unit, scalar, plus, minus, times, word, x, y, rels,
        sliders, i, j, level, ascending, descending, s, p;
    z := IndeterminatesOfFunctionField(field); q := z[1]; t := z[2];
    unit := One(field);
    scalar := c -> CleanNP([[[]], [unit*c]]);
    plus := function(p, r) return AddNP(p, r, unit, unit); end;
    minus := function(p, r) return AddNP(p, r, unit, -unit); end;
    times := function(p, r)
        local answer, a, b;
        answer := [[], []];
        for a in [1..Length(p[1])] do for b in [1..Length(r[1])] do
            Add(answer[1], Concatenation(p[1][a], r[1][b]));
            Add(answer[2], p[2][a]*r[2][b]);
        od; od;
        return CleanNP(answer);
    end;
    word := w -> [ [w], [unit] ];
    x := i -> word([i]); y := i -> word([i+3]);
    rels := [];
    for i in [1..3] do
        Add(rels, minus(times(x(i), y(i)), scalar(1)));
        Add(rels, minus(times(y(i), x(i)), scalar(1)));
    od;
    for i in [1..2] do for j in [i+1..3] do
        if j-i = 1 then
            Add(rels, minus(word([i,j,i]), word([j,i,j])));
            Add(rels, minus(word([i+3,j+3,i+3]), word([j+3,i+3,j+3])));
        else
            Add(rels, minus(word([i,j]), word([j,i])));
            Add(rels, minus(word([i+3,j+3]), word([j+3,i+3])));
            Add(rels, minus(word([i+3,j]), word([j,i+3])));
            Add(rels, minus(word([i,j+3]), word([j+3,i])));
        fi;
    od; od;
    for i in [1..2] do j := i+1;
        Add(rels, minus(word([i,j,i+3]), word([j+3,i,j])));
        Add(rels, minus(word([i+3,j+3,i]), word([j,i+3,j+3])));
        Add(rels, minus(word([i,j+3,i+3]), word([j+3,i+3,j])));
        Add(rels, minus(word([i+3,j,i]), word([j,i,j+3])));
    od;
    sliders := [];
    sliders[2] := List([1..3], j -> minus(plus(times(scalar(q), y(j)), scalar(1-q)), x(j)));
    for level in [3..4] do
        sliders[level] := [];
        for j in [1..5-level] do
            ascending := [j..j+level-2];
            s := minus(times(scalar(q^(level-1)), word(ascending+3)), word(ascending));
            sliders[level][j] := times(s, sliders[level-1][j]);
        od;
    od;
    for j in [1..2] do
        p := minus(minus(sliders[2][j+1], times(scalar(q), word([j+3,j+4]))),
            times(scalar(-1), word([j,j+1])));
        Add(rels, times(p, sliders[2][j]));
    od;
    p := minus(minus(times(scalar(q^2), word([5,6])), word([2,3])),
        minus(times(scalar(q^2), word([4,5,6])), word([1,2,3])));
    Add(rels, times(p, sliders[3][1]));
    for j in [1..3] do
        p := minus(x(j), scalar(t));
        Add(rels, times(times(p,p), sliders[2][j]));
    od;
    for level in [3..4] do for j in [1..5-level] do
        descending := Reversed([j..j+level-2]);
        p := minus(word(descending), times(scalar(q^(level-1)), word(descending+3)));
        Add(rels, times(p, sliders[level][j]));
    od; od;
    return rec(inputs := rels, x3 := sliders[3][1],
        w := times(minus(x(1), scalar(t)), sliders[2][1]), times := times);
end;;

SquareZero.verifyIdentities := function(data, field)
    local pending, proved, progress, index, term, source, residual;
    pending := [1..Length(data.nodes)]; proved := [];
    while Length(pending) > 0 do
        progress := false;
        for index in ShallowCopy(pending) do
            if ForAll(data.proofs[index], term -> term[3] = 0 or term[4] in proved) then
                residual := [[], []];
                for term in data.proofs[index] do
                    if term[3] = 0 then source := data.inputs[term[4]];
                    else source := data.nodes[term[4]]; fi;
                    residual := AddNP(residual, BimulNP(term[2], source, term[5]), One(field), term[1]);
                od;
                if residual <> data.nodes[index] then SquareZero.fail(Concatenation("local identity ", String(index))); fi;
                Add(proved, index); Remove(pending, Position(pending, index)); progress := true;
            fi;
        od;
        if not progress then SquareZero.fail("cyclic derivation references"); fi;
    od;
    return Length(proved);
end;;

SquareZero.verifyQuotient := function(basis, relations, x3, w, dimension, generators, multiply)
    local words, product, term, count;
    if not ForAll(relations, p -> StrongNormalFormNP(p, basis) = [[], []]) then
        SquareZero.fail("original relation has nonzero normal form");
    fi;
    Print("PASS|original_relations_reduce_to_zero|", Length(relations), "\n");
    if not IsGrobnerBasis(basis) then SquareZero.fail("candidate is not a Groebner basis"); fi;
    Print("PASS|complete_groebner_basis|", Length(basis), "\n");
    # Cap exceeds the expected dimension; hitting it would not prove finiteness.
    words := BaseQA(basis, generators, dimension+1);
    if Length(words) <> dimension then SquareZero.fail("unexpected standard-word count"); fi;
    Print("PASS|dimension_over_Q(q,t)|", dimension, "\n");
    if StrongNormalFormNP(x3, basis) <> [[], []] then SquareZero.fail("X3 is nonzero"); fi;
    w := StrongNormalFormNP(w, basis);
    if w = [[], []] then SquareZero.fail("w is zero"); fi;
    count := 0;
    for term in words do
        product := multiply(BimulNP([], w, term[1][1]), w);
        if StrongNormalFormNP(product, basis) <> [[], []] then SquareZero.fail("w*b*w is nonzero"); fi;
        count := count + 1;
    od;
    Print("PASS|X3_zero|w_nonzero|square_zero_ideal_tests|", count, "\n");
    return count;
end;;

SquareZero.check := function(path)
    local field, data, defining, started, proved, basis;
    field := FunctionField(Rationals, ["q", "t"]);
    started := Runtime(); data := SquareZero.read(path, field);
    defining := SquareZero.presentation(field);
    if data.inputs <> defining.inputs then SquareZero.fail("input presentation differs"); fi;
    Print("PASS|original_presentation|31\n");
    proved := SquareZero.verifyIdentities(data, field);
    Print("PASS|acyclic_local_identities|", proved, "|ms|", Runtime()-started, "\n");
    Print("PASS|auxiliary_input_bound|", SquareZero.auxiliaryInputCount, "\n");
    basis := data.nodes{[1..data.basisCount]};
    SquareZero.verifyQuotient(basis, data.inputs, defining.x3, defining.w, 186, 6, defining.times);
    Print("PASS: exact four-strand square-zero verification over Q(q,t) | ms=", Runtime()-started, "\n");
end;;

if not IsBound(SquareZeroLibraryOnly) or SquareZeroLibraryOnly <> true then
    if not IsBound(SquareZeroProofFile) then
        SquareZeroProofFile := "runs/four-strand-square-zero/proof.txt";
    fi;
    SquareZero.check(SquareZeroProofFile);
    QUIT_GAP(0);
fi;
