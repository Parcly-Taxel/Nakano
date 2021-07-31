import re
from math import gcd, isqrt
from collections import Counter
from functools import reduce
import numpy as np
import lifelib
rng = np.random.default_rng()
rule_re = re.compile(r"rule\s*=\s*([A-Za-z0-9/_-]+)\s*$", re.M)

def factors(n):
    """Yield the factors of n in ascending order."""
    rtn = isqrt(n)
    smalls = list(filter(lambda k: n % k == 0, range(1, rtn + 1)))
    larges = [n // k for k in smalls]
    if rtn * rtn == n:
        smalls.pop()
    yield from smalls
    yield from reversed(larges)

def equalundertf(pat, tf, dt=0):
    """Checks whether the given pattern equals itself when transformed
    and optionally evolved a given number of generations."""
    return pat.centre() == pat[dt](tf).centre()

def getsymmetries(pat, p):
    """Return the periodic pattern's apgsearch symmetry (looking at only one phase),
    temporal symmetry (considering half/quarter-period transformations if applicable)
    and mod."""
    b = equalundertf(pat, "flip_x")
    f = equalundertf(pat, "flip_y")
    r = equalundertf(pat, "swap_xy")
    l = equalundertf(pat, "swap_xy_flip")
    hp = 0 if p % 2 else p // 2
    qp = 0 if p % 4 else p // 4
    if b and r:
        sym = "D8_" + "41"[pat.bounding_box[2] % 2]
        return (sym, sym, p)
    if b and f:
        dx, dy = pat.bounding_box[2:]
        sym = "D4_+" + "421"[dx%2 + dy%2]
        if hp and equalundertf(pat, "rot90", hp):
            return (sym, "D8_" + sym[4], hp)
        return (sym, sym, p)
    if r and l:
        sym = "D4_x" + "41"[pat.bounding_box[2] % 2]
        if hp and equalundertf(pat, "rot90", hp):
            return (sym, "D8_" + sym[4], hp)
        return (sym, sym, p)
    if b or f:
        sym = "D2_+" + "21"[pat.bounding_box[2+f] % 2]
        if hp and equalundertf(pat, "flip_y" if b else "flip_x", hp):
            dx, dy = (pat + pat[hp]).bounding_box[2:]
            return (sym, "D4_+" + "421"[dx%2 + dy%2], hp)
        return (sym, sym, p)
    if r or l:
        sym = "D2_x"
        if hp and equalundertf(pat, "swap_xy_flip" if r else "swap_xy", hp):
            return (sym, "D4_x" + "41"[(pat + pat[hp]).bounding_box[2] % 2], hp)
        return (sym, sym, p)
    if equalundertf(pat, "rot90"):
        sym = "C4_" + "41"[pat.bounding_box[2] % 2]
        if hp and equalundertf(pat, "flip_x", hp):
            return (sym, "D8_" + sym[3], hp)
        return (sym, sym, p)
    if equalundertf(pat, "rot180"):
        dx, dy = pat.bounding_box[2:]
        sym = "C2_" + "421"[dx%2 + dy%2]
        if hp:
            if equalundertf(pat, "rot90", hp):
                return (sym, "C4_" + sym[3], hp)
            if equalundertf(pat, "flip_x", hp):
                return (sym, "D4_+" + sym[3], hp)
            if equalundertf(pat, "swap_xy", hp):
                return (sym, "D4_x" + sym[3], hp)
        return (sym, sym, p)
    sym = "C1"
    if qp and (equalundertf(pat, "rot90", qp) or equalundertf(pat, "rot270", qp)):
        pat4 = pat + pat[qp] + pat[2*qp] + pat[3*qp]
        return (sym, "C4_" + "41"[pat4.bounding_box[2] % 2], qp)
    if hp:
        if equalundertf(pat, "rot180", hp):
            dx, dy = (pat + pat[hp]).bounding_box[2:]
            return (sym, "C2_" + "421"[dx%2 + dy%2], hp)
        if equalundertf(pat, "flip_x", hp):
            return (sym, "D2_+" + "21"[(pat + pat[hp]).bounding_box[2] % 2], hp)
        if equalundertf(pat, "flip_y", hp):
            return (sym, "D2_+" + "21"[(pat + pat[hp]).bounding_box[3] % 2], hp)
        if equalundertf(pat, "swap_xy", hp) or equalundertf(pat, "swap_xy_flip", hp):
            return (sym, "D2_x", hp)
    return (sym, sym, p)

dxys = ((-1,-1), (0,-1), (1,-1), (-1,0), (1,0), (-1,1), (0,1), (1,1))
mooreconfs = np.array((0, 1, 2, 6, 1, 3, 6, 13, 2, 6, 4, 12, 5, 14, 17, 22, 2, 5, 4, 17, 6, 14, 12, 22, 7, 18, 10, 28, 18, 23, 28, 36, 1, 3, 5, 14, 8, 9, 16, 24, 6, 13, 17, 22, 16, 24, 30, 35, 5, 15, 11, 21, 16, 25, 26, 40, 18, 29, 27, 37, 31, 41, 39, 45, 2, 5, 7, 18, 5, 15, 18, 29, 4, 17, 10, 28, 11, 21, 27, 37, 4, 11, 10, 27, 17, 21, 28, 37, 10, 27, 20, 32, 27, 38, 32, 42, 6, 14, 18, 23, 16, 25, 31, 41, 12, 22, 28, 36, 26, 40, 39, 45, 17, 21, 27, 38, 30, 34, 39, 44, 28, 37, 32, 42, 39, 44, 47, 48, 1, 8, 5, 16, 3, 9, 14, 24, 5, 16, 11, 26, 15, 25, 21, 40, 6, 16, 17, 30, 13, 24, 22, 35, 18, 31, 27, 39, 29, 41, 37, 45, 3, 9, 15, 25, 9, 19, 25, 33, 14, 24, 21, 40, 25, 33, 34, 43, 14, 25, 21, 34, 24, 33, 40, 43, 23, 41, 38, 44, 41, 46, 44, 49, 6, 16, 18, 31, 14, 25, 23, 41, 17, 30, 27, 39, 21, 34, 38, 44, 12, 26, 28, 39, 22, 40, 36, 45, 28, 39, 32, 47, 37, 44, 42, 48, 13, 24, 29, 41, 24, 33, 41, 46, 22, 35, 37, 45, 40, 43, 44, 49, 22, 40, 37, 44, 35, 43, 45, 49, 36, 45, 42, 48, 45, 49, 48, 50))
hsegments = ((0, 1), (1, 2), (3, 6), (9, 10), (19, 13), (32, 10), (42, 6), (48, 2), (50, 1))

def henselcode(transtr):
    """Encode the given transition specification (a length-51 string of 0s and 1s)
    in Hensel notation."""
    res = []
    for (n, (start, length)) in enumerate(hsegments):
        sn = str(n)
        seg = transtr[start:start+length]
        if "0" not in seg:
            res.append(sn)
        elif "1" not in seg:
            continue
        else:
            target, prefix = ("1", "") if seg.count("1") < (length + 3) // 2 else ("0", "-")
            res.append(sn)
            res.append(prefix)
            noted = sorted("cekainyqjrtwz"[i] for i in range(length) if seg[i] == target)
            res.append("".join(noted))
    return "".join(res)

def minmaxrule(pat, ngens):
    """Compute the isotropic (and outer-totalistic if applicable) minrule and maxrule
    of the given not-necessarily-periodic pattern over the given number of generations."""
    moore = pat.owner.pattern("3o$3o$3o")(-1,-1)
    prevpat, nextpat = pat[0], pat[1]
    reqtrans = {-1}
    for _ in range(ngens):
        cc = prevpat.convolve(moore).coords()
        A = prevpat[np.concatenate([cc + dxy for dxy in dxys])].reshape(8, -1)
        B = mooreconfs[np.matmul((128, 64, 32, 16, 8, 4, 2, 1), A, dtype=int)]
        reqtrans.update((2*B + prevpat[cc].astype(int) + 1) *
                        (2*nextpat[cc].astype(int) - 1))
        prevpat, nextpat = nextpat, nextpat[1]
    reqB = "".join("1" if 2*t+1 in reqtrans else "0" if -(2*t+1) in reqtrans else "-" for t in range(51))
    reqS = "".join("1" if 2*t+2 in reqtrans else "0" if -(2*t+2) in reqtrans else "-" for t in range(51))
    minB = henselcode(reqB.replace("-", "0"))
    minS = henselcode(reqS.replace("-", "0"))
    isominrule = f"B{minB}/S{minS}"
    maxB = henselcode(reqB.replace("-", "1"))
    maxS = henselcode(reqS.replace("-", "1"))
    isomaxrule = f"B{maxB}/S{maxS}"

    segsB = [reqB[start:start+length] for (start, length) in hsegments]
    segsS = [reqS[start:start+length] for (start, length) in hsegments]
    B_otincompat = any("0" in seg and "1" in seg for seg in segsB)
    S_otincompat = any("0" in seg and "1" in seg for seg in segsS)
    if B_otincompat or S_otincompat:
        return (isominrule, isomaxrule, None, None)
    reqB = "".join(seg.replace("-", "0" if "0" in seg else "1" if "1" in seg else "-") for seg in segsB)
    reqS = "".join(seg.replace("-", "0" if "0" in seg else "1" if "1" in seg else "-") for seg in segsS)
    minB = henselcode(reqB.replace("-", "0"))
    minS = henselcode(reqS.replace("-", "0"))
    otminrule = f"B{minB}/S{minS}"
    maxB = henselcode(reqB.replace("-", "1"))
    maxS = henselcode(reqS.replace("-", "1"))
    otmaxrule = f"B{maxB}/S{maxS}"
    return (isominrule, isomaxrule, otminrule, otmaxrule)

def analyse(pat_string, rule="b3s23"):
    """Compute cell periods and other statistics of the given oscillator/spaceship (as an RLE or apggcode)
    like the old Oscillizer. The pattern may be in any 2-state isotropic rule not containing B0.
    Return a dictionary containing relevant statistics and other information:
    res["p"] -> period
    res["dxy"] -> displacement
    res["apgcode"] -> apgcode
    res["syms"] -> (apgsearch symmetry, temporal symmetry, mod); see getsymmetries()
    res["pops"] -> list of successive populations of the object
    res["popstats"] -> (min. population, max. population, average population, sum of populations across a period)
    res["heats"] -> list of successive heats of the object
    res["heatstats"] -> (min. heat, max. heat, average heat, sum of heats across a period)
    res["bb"] -> bounding box (width, height, cell count)
    res["rules"] -> (isotropic minrule, isotropic maxrule, OT minrule, OT maxrule);
    the last two entries are None if the pattern does not work in any OT rule

    For oscillators only:
    res["cellperiods"] -> dictionary mapping active cells to their periods (on cells in the provided pattern
    have their periods negated)
    res["cellpcounts"] -> dictionary mapping seen cell periods to frequency
    res["cellstats"] -> (rotor cells, stator cells, active cells, strict (full-period) rotor cells)
    res["volstats"] -> (volatility, strict volatility)
    res["tempstats"] -> (temperature, rotor temperature)"""
    if m := rule_re.search(pat_string):
        rule = m[1]
    lt = lifelib.load_rules(rule).lifetree(n_layers=1)
    empty = lt.pattern()
    pat = lt.pattern(pat_string)

    res = {}
    lcmp = res["p"] = pat.period
    dxy = res["dxy"] = pat.displacement
    res["apgcode"] = pat.apgcode
    osc = dxy == (0, 0)
    res["syms"] = getsymmetries(pat, lcmp)
    phases = [pat[i] for i in range(lcmp)]
    pops = [phase.population for phase in phases]
    heats = [(phases[i] ^ phases[i-1]).population for i in range(lcmp)]
    sp = sum(pops)
    sh = sum(heats)
    res["pops"] = pops
    res["popstats"] = (min(pops), max(pops), sp / lcmp, sp)
    res["heats"] = heats
    res["heatstats"] = (min(heats), max(heats), sh / lcmp, sh)

    if osc:
        cellperiods = {}
        pcounts = Counter()
        remcells = sum(phases, start=empty)
        actives = remcells.population
        width, height = remcells.bounding_box[2:]
        res["bb"] = (width, height, width * height)
        for subp in factors(lcmp):
            if subp == 1:
                subpcells = reduce(lambda x, y: x & y, phases)
            if subp == lcmp:
                if remcells:
                    cellperiods.update({(x, y): subp for (x, y) in remcells.coords().tolist()})
                    pcounts[subp] = remcells.population
                break
            else:
                subpcells = remcells - sum((phases[i] ^ phases[i-subp] for i in range(lcmp-subp)), start=empty)
            if subpcells:
                cellperiods.update({(x, y): subp for (x, y) in subpcells.coords().tolist()})
                pcounts[subp] = subpcells.population
                remcells -= subpcells
        for (c, p) in cellperiods.items():
            if pat[c]:
                cellperiods[c] = -p
        res["cellperiods"] = cellperiods
        res["cellpcounts"] = [(p, pcounts[p]) for p in sorted(pcounts, reverse=True)]
        stators = pcounts[1]
        rotors = actives - stators
        strict_rotors = min(rotors, pcounts[lcmp])
        res["cellstats"] = (rotors, stators, actives, strict_rotors)
        res["volstats"] = (rotors / actives, strict_rotors / actives)
        res["tempstats"] = (sh / (lcmp * actives), sh / (lcmp * rotors) if rotors else 0)
    else:
        bb = max((pat[i].bounding_box[2:] for i in range(lcmp)), key=lambda x: x[0] * x[1])
        res["bb"] = (bb[0], bb[1], bb[0] * bb[1])

    res["rules"] = minmaxrule(pat, res["syms"][2])
    return res

def speedstring(dmaj, dmin, period):
    """Return a short string describing the speed with the given parameters.
    Conditions are dmaj >= dmin >= 0, dmaj > 0, period >= 1."""
    if dmaj > dmin > 0:
        numerator = f"({dmaj},{dmin})"
        terminator = ""
    elif dmaj == dmin:
        numerator = f"{dmaj}" if dmaj > 1 else ""
        terminator = "d"
    else:
        numerator = f"{dmaj}" if dmaj > 1 else ""
        terminator = "o"
    denominator = f"/{period}" if period > 1 else ""
    return f"{numerator}c{denominator}{terminator}"

def resultprint(res):
    """Pretty-print the summary statistics of the given results dictionary.
    The formatting for oscillators closely matches the old Oscillizer."""
    period = res["p"]
    apgcode = res["apgcode"]
    if res["dxy"] == (0, 0):
        print(f"p{period} oscillator" if period > 1 else "still life", "–", apgcode)
        sym, tempsym, mod = res["syms"]
        tempstr = "" if sym == tempsym else f" (temporal {tempsym})"
        print(f"Symmetry {sym}{tempstr}, mod {mod}")
        print("Cell periods:")
        for (cp, count) in res["cellpcounts"]:
            print(f"p{cp} – {count}")
        minpop, maxpop, avgpop, _ = res["popstats"]
        print(f"Population {minpop}–{maxpop}, average {avgpop:.2f}")
        rotors, stators, total, strict_rotors = res["cellstats"]
        print(f"{rotors} rotor cells ({strict_rotors} full-period), {stators} stator cells, {total} total")
        vol, svol = res["volstats"]
        print(f"Volatility {vol:.4f} ({svol:.4f} strict)")
        minheat, maxheat, avgheat, _ = res["heatstats"]
        print(f"Heat {minheat}–{maxheat}, average {avgheat:.2f}")
        temp, rtemp = res["tempstats"]
        print(f"Temperature {temp:.4f} ({rtemp:.4f} rotor)")
    else:
        dx, dy = res["dxy"]
        dmaj, dmin = max(abs(dx), abs(dy)), min(abs(dx), abs(dy))
        unsimp_speed = speedstring(dmaj, dmin, period)
        if (k := gcd(dx, dy, period)) > 1:
            simp_speed = speedstring(dmaj//k, dmin//k, period//k)
            print(f"{unsimp_speed} ({simp_speed} simplified) spaceship – {apgcode}")
        else:
            print(f"{unsimp_speed} spaceship – {apgcode}")
        sym, tempsym, mod = res["syms"]
        tempstr = " (glide-symmetric)" if mod < period else ""
        print(f"Symmetry {sym}{tempstr}, mod {mod}")
        minpop, maxpop, avgpop, _ = res["popstats"]
        print(f"Population {minpop}–{maxpop}, average {avgpop:.2f}")
        minheat, maxheat, avgheat, _ = res["heatstats"]
        print(f"Heat {minheat}–{maxheat}, average {avgheat:.2f}")
    bw, bh, bc = res["bb"]
    print(f"Bounding box {bw}×{bh} = {bc}")
    isominrule, isomaxrule, otminrule, otmaxrule = res["rules"]
    isostr = f"{isominrule}" if isominrule == isomaxrule else f"{isominrule} – {isomaxrule}"
    if otminrule == None:
        otstr = "no OT"
    else:
        span = f"{otminrule}" if otminrule == otmaxrule else f"{otminrule} – {otmaxrule}"
        otstr = "OT " + span
    print(f"Rules {isostr} ({otstr})")

start_colours = ("ffffff", "1c92cd", "0ab87b", "e86075", "f8f290",
                 "ba4117", "d91e9b", "aeaeae", "bef0e9", "428c28",
                 "6f1044", "76adf4", "2f5963", "d9b790")

def blackorwhite(col):
    """Given an RGB colour specified as a 3-sequence of integers from 0 to 255,
    determine whether white or black text is more appropriate for it."""
    nd = [x/255 for x in col]
    dd = [x/12.92 if x <= 0.03928 else ((x+0.055)/1.055) ** 2.4 for x in nd]
    return 0 if 0.2126*dd[0] + 0.7152*dd[1] + 0.0722*dd[2] > 0.1791 else 255

def periodmap(cellperiods, outfn="osc.png", scale=16):
    """Save a colour-coded map of the active cells and print the
    corresponding periods. Requires Pillow; scale sets the size of a cell."""
    from PIL import Image
    if scale < 1:
        return
    cells = np.array([[c[0], c[1], p] for (c, p) in cellperiods.items()])
    cells[:,0] -= min(cells[:,0])
    cells[:,1] -= min(cells[:,1])
    width = max(cells[:,0]) + 1
    height = max(cells[:,1]) + 1
    unique_periods = np.unique(abs(cells[:,2]))[::-1]
    colourmap = {}
    for (n, p) in enumerate(unique_periods):
        if p == 1:
            colourmap[p] = ((0, 0, 0), (0, 0, 0))
        elif n < 14:
            col = [int(start_colours[n][i:i+2], 16) for i in (0, 2, 4)]
            colourmap[p] = (tuple(col), blackorwhite(col))
        else:
            col = rng.integers(16, 240, 3)
            colourmap[p] = (tuple(col), blackorwhite(col))
    for (k, v) in colourmap.items():
        print(f"p{k} – {v[0][0]:02x}{v[0][1]:02x}{v[0][2]:02x}")
    field = np.full(((height+2)*scale, (width+2)*scale, 3), 204, dtype=np.uint8)
    large_pixel = np.ones((scale, scale, 3), dtype=np.uint8)
    border = max((scale+1) // 4, 1)
    sps = max(scale - 2*border, 0)
    small_pixel = np.ones((sps, sps, 3), dtype=np.uint8)
    for row in cells:
        col, bw = colourmap[abs(row[2])]
        field[(row[1]+1)*scale:(row[1]+2)*scale,(row[0]+1)*scale:(row[0]+2)*scale,:] = large_pixel * col
        if row[2] < -1:
            field[(row[1]+1)*scale+border:(row[1]+2)*scale-border,
                  (row[0]+1)*scale+border:(row[0]+2)*scale-border,:] = small_pixel * bw
    Image.fromarray(field).save(outfn)

def n(pat_string, rule="b3s23", outfn="osc.png", scale=16):
    """Interactive function to analyse a pattern and display the results."""
    res = analyse(pat_string, rule)
    resultprint(res)
    if res["dxy"] == (0, 0):
        periodmap(res["cellperiods"], outfn, scale)
