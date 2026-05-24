"""Count formal connected vacuum Wick contractions of density-density vertices.

Interaction expansion denominator is n! 2**n.  Fix Coulomb perfect matching
(0,1),(2,3),... and enumerate the fermion Wick permutations.  Their cycles
are directed fermion loops.  Only one C4 and n-2 C2 loops are retained.
All retained diagrams at given n have the same fermion-loop sign.
This checks combinatorial magnitudes, not frequency integrals or fixed-density
counterterm diagrams.  Spin traces belong inside the loop factors.
"""
import itertools
import math
from collections import Counter

# Exact contraction counts for the topologies stated in the manuscript.
expected_counts = {
    2: {'bare_exchange': 2},
    3: {'one_chord_decorated': 24},
    4: {'one_chord_decorated': 192, 'both_chords_one_bubble_each': 96},
}

for n in (2, 3, 4):
    N = 2*n
    counts = Counter()
    for p in itertools.permutations(range(N)):
        seen = set()
        cycles = []
        for i in range(N):
            if i not in seen:
                cyc = []
                j = i
                while j not in seen:
                    cyc.append(j)
                    seen.add(j)
                    j = p[j]
                cycles.append(cyc)
        if sorted(map(len, cycles)) != sorted([4]+[2]*(n-2)):
            continue
        core = next(c for c in cycles if len(c) == 4)
        bubbles = [c for c in cycles if len(c) == 2]
        owner = {v:k for k,c in enumerate([core]+bubbles) for v in c}
        adj = {k:[] for k in range(n-1)}
        bad = False
        bare = []
        for i in range(0,N,2):
            a, b = owner[i], owner[i+1]
            adj[a].append(b)
            adj[b].append(a)
            if a == b:
                if a != 0:
                    bad = True
                    break
                bare.append((i,i+1))
        if bad:
            continue
        # Require connected quotient graph of fermion loops.
        seen_nodes = {0}
        stack = [0]
        while stack:
            for j in adj[stack.pop()]:
                if j not in seen_nodes:
                    seen_nodes.add(j)
                    stack.append(j)
        if len(seen_nodes) != n-1:
            continue
        if n == 2:
            # Opposite chords distinguish exchange from the other C4 closures.
            if len(bare) == 2 and all(
                (core.index(a)-core.index(b)) % 4 == 2 for a,b in bare
            ):
                counts['bare_exchange'] += 1
        elif len(bare) == 1 and (
            core.index(bare[0][0])-core.index(bare[0][1])
        ) % 4 == 2:
            counts['one_chord_decorated'] += 1
        elif n == 4 and len(bare) == 0 and all(
            set(adj[k]) == {0} for k in (1,2)
        ):
            # Both bubbles must connect opposite core vertices.
            ok = True
            for b in bubbles:
                endpoints = [v^1 for v in b]
                if (core.index(endpoints[0])-core.index(endpoints[1])) % 4 != 2:
                    ok = False
            if ok:
                counts['both_chords_one_bubble_each'] += 1
    denom = math.factorial(n)*2**n
    assert dict(counts) == expected_counts[n], (n, dict(counts), expected_counts[n])
    print('order',n,'Wick denominator',denom,'counts',dict(counts),
          'coefficients',{k:v/denom for k,v in counts.items()})
