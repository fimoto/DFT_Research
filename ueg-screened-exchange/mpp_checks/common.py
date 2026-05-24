"""Shared helpers for the mpp_checks scripts.

Each script reproduces numbers printed in the Mellin-pole manuscript
"Revisiting the r_s expansion of the uniform electron gas in d dimensions:
log r_s terms from coincident Mellin poles" (manuscript version v4.21) and its
Supporting Information (SI).  Every reproduced number is stored with the place
where it is printed, the printed value, the value computed here and a verdict.
"""
import json, os, sys, time, platform

MANUSCRIPT = "MPP v4.21 (submitted version): main text Eqs. (1)-(59), SI Eqs. (S1)-(S104)"
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


class Recorder:
    """Collects checks of one script and writes results/<name>.json."""

    def __init__(self, name, title):
        self.name, self.title = name, title
        self.items, self.t0 = [], time.time()

    def _add(self, d):
        self.items.append(d)
        flag = "PASS" if d["pass"] else "FAIL"
        print(f"[{flag}] {d['where']}: {d['what']}\n       computed {d['computed']}  printed {d['printed']}  ({d['criterion']})", flush=True)

    def digits(self, where, what, computed, printed):
        """The printed decimal string must be the correctly rounded or the truncated
        value of `computed` (any sign).  `computed` may be an mpmath number."""
        s = printed.strip()
        neg = s.startswith("-")
        body = s.lstrip("+-")
        nd = len(body.split(".")[1]) if "." in body else 0
        import mpmath as mp
        with mp.workdps(60):
            c = mp.mpf(computed)
            p = mp.mpf(s)
            scale = mp.mpf(10) ** (-nd)
            rounded = abs(c - p) <= scale / 2 * (1 + mp.mpf(10) ** -30)
            trunc = (abs(c) - abs(p)) >= -scale * mp.mpf(10) ** -30 and (abs(c) - abs(p)) < scale and (c < 0) == neg
            ok = bool(rounded or trunc)
        self._add(dict(where=where, what=what, computed=mp.nstr(c, nd + 6), printed=printed,
                       criterion="printed digits = rounded or truncated value", pass_kind=("rounded" if rounded else "truncated" if trunc else "no"), **{"pass": ok}))
        return ok

    def close(self, where, what, computed, target, tol, rel=False, printed=None):
        import mpmath as mp
        with mp.workdps(60):
            c, t = mp.mpf(computed), mp.mpf(target)
            dev = abs(c - t) / (abs(t) if rel else 1)
            ok = bool(dev <= tol)
        self._add(dict(where=where, what=what, computed=mp.nstr(c, 25), printed=printed if printed is not None else mp.nstr(t, 25),
                       deviation=mp.nstr(dev, 3), criterion=f"{'relative' if rel else 'absolute'} deviation <= {mp.nstr(mp.mpf(tol), 3)}", **{"pass": ok}))
        return ok

    def bound(self, where, what, value, limit, printed):
        """value (a deviation found here) must not exceed the printed bound."""
        import mpmath as mp
        ok = bool(mp.mpf(value) <= mp.mpf(limit))
        self._add(dict(where=where, what=what, computed=mp.nstr(mp.mpf(value), 4), printed=printed,
                       criterion=f"deviation found here <= printed bound {mp.nstr(mp.mpf(limit), 3)}", **{"pass": ok}))
        return ok

    def exact(self, where, what, computed, target):
        ok = bool(computed == target)
        self._add(dict(where=where, what=what, computed=str(computed), printed=str(target), criterion="exact equality", **{"pass": ok}))
        return ok

    def note(self, where, what, value):
        self.items.append(dict(where=where, what=what, computed=str(value), printed=None, criterion="informational", **{"pass": True}))

    def write(self, extra=None):
        os.makedirs(RESULTS, exist_ok=True)
        import mpmath, numpy, scipy
        out = dict(script=self.name + ".py", title=self.title, manuscript=MANUSCRIPT,
                   python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, mpmath=mpmath.__version__,
                   seconds=round(time.time() - self.t0, 1), n_checks=len(self.items),
                   n_fail=sum(not i["pass"] for i in self.items), checks=self.items)
        if extra:
            out.update(extra)
        with open(os.path.join(RESULTS, self.name + ".json"), "w") as f:
            json.dump(out, f, indent=1)
        print(f"== {self.name}: {len(self.items)} checks, {out['n_fail']} failed, {out['seconds']} s", flush=True)
        return out["n_fail"]
