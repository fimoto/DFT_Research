#!/usr/bin/env python3
"""Enumerate Wick contractions for a C4 exchange core and C2 bubble chains.

No third-party packages are needed.  All contractions are generated directly
with cycle type (4,2,...,2), avoiding a factorial scan of irrelevant permutations.
The Coulomb matching is fixed to (0,1),(2,3),... . A directed C4 is generated once
by placing its smallest label first. Each residual perfect matching is one set
of directed C2 loops. Core endpoints are classified as opposite (exchange) or
adjacent (self-energy-type closure). A disconnected bubble ring is rejected.

The interaction expansion denominator is 2**n * n!. Signs, spin traces, and
integrals are common diagram factors and are NOT evaluated. Exact fractions
below are unsigned symmetry coefficients in the connected-vacuum expansion.

Usage:
    python verify_exchange_chain_coefficients.py
    python verify_exchange_chain_coefficients.py --max-order 6 --json report.json
    python verify_exchange_chain_coefficients.py --core-pattern exchange
"""
import argparse
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def pairings(items):
    if not items:
        yield ()
        return
    a, *rest = items
    for j, b in enumerate(rest):
        reduced = rest[:j] + rest[j+1:]
        for tail in pairings(reduced):
            yield ((a, b),) + tail


def classify(core, bubbles, n):
    core_pos = {v: i for i, v in enumerate(core)}
    bubble_mate = {}
    for a, b in bubbles:
        bubble_mate[a] = b
        bubble_mate[b] = a
    used_core = set()
    used_bubble_vertices = set()
    lengths = []
    separations = []
    for start in core:
        if start in used_core:
            continue
        current = start ^ 1  # partner on the fixed Coulomb line
        length = 0
        path_vertices = set()
        while current not in core_pos:
            if current in path_vertices:
                return None
            other = bubble_mate[current]
            path_vertices.update((current, other))
            current = other ^ 1
            length += 1
        if current == start or current in used_core:
            return None
        separations.append((core_pos[current] - core_pos[start]) % 4)
        used_core.update((start, current))
        used_bubble_vertices.update(path_vertices)
        lengths.append(length)
    if len(used_core) != 4 or len(used_bubble_vertices) != 2*n-4:
        return None
    if all(d == 2 for d in separations):
        pattern = 'exchange'
    elif all(d in (1, 3) for d in separations):
        pattern = 'adjacent'
    else:
        raise AssertionError(('inconsistent_core_pairing', core, separations))
    return pattern, tuple(sorted(lengths))


def enumerate_order(n, selected_pattern):
    vertices = tuple(range(2*n))
    counts = Counter()
    generated = 0
    for core_set in itertools.combinations(vertices, 4):
        remaining = [v for v in vertices if v not in core_set]
        for tail in itertools.permutations(core_set[1:]):
            core = (core_set[0],) + tail
            for bubbles in pairings(remaining):
                generated += 1
                key = classify(core, bubbles, n)
                if key is not None and selected_pattern in ('both', key[0]):
                    counts[key] += 1
    denominator = (2**n) * math.factorial(n)
    rows = []
    for (pattern, lengths), count in sorted(counts.items()):
        coefficient = Fraction(count, denominator)
        base = 4 if pattern == 'exchange' else 2
        expected = Fraction(1, base if lengths[0] == lengths[1] else base//2)
        if coefficient != expected:
            raise AssertionError((n, lengths, count, coefficient, expected))
        rows.append({
            'core_pattern': pattern,
            'bubble_chain_lengths': list(lengths),
            'wick_contractions': count,
            'expansion_denominator': denominator,
            'symmetry_coefficient': str(coefficient),
            'automorphism_prediction': str(expected),
        })
    patterns = ('exchange', 'adjacent') if selected_pattern == 'both' else (selected_pattern,)
    expected_keys = {(pattern, (a, n-2-a)) for pattern in patterns
                     for a in range((n-2)//2 + 1)}
    if set(counts) != expected_keys:
        raise AssertionError(('missing_topologies', n, set(counts), expected_keys))
    return {'interaction_order': n, 'generated_cycle_permutations': generated,
            'topologies': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-order', type=int, default=4)
    parser.add_argument('--core-pattern', choices=('exchange', 'adjacent', 'both'), default='both')
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    if not 2 <= args.max_order <= 7:
        parser.error('--max-order must be between 2 and 7 (cost grows rapidly)')
    report = {'status': 'passed', 'convention': 'unsigned connected-vacuum symmetry coefficients',
              'orders': [enumerate_order(n, args.core_pattern) for n in range(2, args.max_order+1)]}
    print('n  core      (chain lengths)  Wick count / (2^n n!)  coefficient')
    for order in report['orders']:
        for row in order['topologies']:
            print(f"{order['interaction_order']}  {row['core_pattern']:9} "
                  f"{str(tuple(row['bubble_chain_lengths'])):15} "
                  f"{row['wick_contractions']:7} / {row['expansion_denominator']:<8} "
                  f"{row['symmetry_coefficient']}")
    if args.json:
        args.json.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
