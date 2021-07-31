![15 default subperiod colours used by Nakano](/logo.svg)

The function `n()` in `nakano.py` is a complete replacement and enhancement of Jason Summers's [Oscillizer](http://entropymine.com/jason/life/oscillizer), which went offline sometime in April 2021 (but that program's source code was still available). Like [Skopje](https://github.com/Parcly-Taxel/Skopje), the whole program started out as a standalone script in [Shinjuku](https://gitlab.com/parclytaxel/Shinjuku) before being split into its own repository and is named after a real-world location – [Nakano](https://en.wikipedia.org/wiki/Nakano,_Tokyo) is a special ward of Tokyo bordering Shinjuku (and the repository split occurred during the Tokyo Olympics). In addition to requiring [lifelib](https://gitlab.com/apgoucher/lifelib) and [NumPy](https://numpy.org), [Pillow](https://python-pillow.org) is needed to generate the colour-coded maps of cell subperiods that were a hallmark of Oscillizer.

```
>>> from nakano import n
>>> n("""x = 41, y = 29, rule = B3/S23
... 12boo$12boobo15bo$16bo13bobo$8bo4bo11bo$7bobo4boboo5b3o4bo7boo$16boo4b
... o7boo6bobo$9bo12boo3b5o3boobobo$8boo17booboo3bobobo$8b5o17bo6bo$boo5b
... ooboo22bobbo$obo6bo28bo$oboboo28bo3bo$bobobo12b3o8bobobbo3bo$3bo14bobo
... 8bobbo5bo$bbobbo3bobo6b3o8bobo3bobbo$bbo5bobbo25bo$bbo3bobbobo23bobobo
... $bbo3bo28boobobo$bbo28bo6bobo$bbobbo22booboo5boo$3bo6bo17b5o$bobobo3b
... ooboo17boo$oboboo3b5o3boo12bo$obo6boo7bo4boo$boo7bo4b3o5boobo4bobo$15b
... o11bo4bo$8bobo13bo$9bo15boboo$27boo!""")
p96 oscillator – xp96_0ggy142ky0oggy1gk2gb3wg8gz34jc2mx7f266x23y011xoo8usy1ggzwvwg608k0sy4sksy2g0hxm2cj43zcis346wgo0ggyc302106gwvzy3g3n1110g0ooy0ckgg0664vex643siczy41xcd042y41y0242
Symmetry C1 (temporal C2_1), mod 48
Cell periods:
p96 – 114
p48 – 1
p12 – 56
p8 – 84
p4 – 56
p3 – 88
p2 – 12
p1 – 72
Population 177–240, average 210.77
411 rotor cells (114 full-period), 72 stator cells, 483 total
Volatility 0.8509 (0.2360 strict)
Heat 115–196, average 162.33
Temperature 0.3361 (0.3950 rotor)
Bounding box 41×33 = 1353
Rules B3/S23-k – B34q5c6ci/S234cy6i7e (OT B3/S23)

p96 – ffffff
p48 – 1c92cd
p12 – 0ab87b
p8 – e86075
p4 – f8f290
p3 – ba4117
p2 – d91e9b
p1 – 000000
```

Running the above command also produces this image:

![p96 oscillator used as an Oscillizer example](/osc-example.png)

----

Nakano provides all the features of Oscillizer and has the following improvements:

* Instead of an RLE an apgcode can be provided.
* Periods, pattern sizes and number of distinct subperiods are essentially unlimited (Oscillizer had a limit of p300, 400×400 and 10 subperiods).
* Spaceships are also acceptable input (i.e. Nakano can handle any periodic object). In that case a subperiod map will not be produced, but relevant statistics are still computed.
* Any isotropic rule not containing B0 can be used. The rule to use can either be specified in the RLE or through the `rule=` parameter of `n()`.
* The image size can be changed using the `scale=` parameter in `n()`; image output can be disabled by setting `scale=0`. Where that image gets saved is also configurable using the `outfn=` option.
* Statistics of interest to contemporary cellular automata enthusiasts, like mod and symmetry type, are computed alongside the traditional ones.

`n()` is a convenience function for the command line, interfacing with the function `analyse(pat_string, rule="b3s23")` that should be used when many patterns need to be analysed together. The latter function returns a dictionary with the following string keys:

Key | Value
--- | -----
`p` | The overall oscillator period. There may be no cells with this period – their presence or absence is indicated by strict volatility (see below).
`dxy` | A pair of integers representing displacement in the x and y directions per period.
`apgcode` | The pattern's apgcode.
`syms` | A triple (base symmetry, temporal symmetry, mod). The first two elements are individually one of the 16 possible "official" [symmetries](https://conwaylife.com/wiki/Symmetry) that may be specified in an apgsearch (excluding `D8_+2` which Conway's Life does not respect). Base symmetry only looks at individual phases while temporal symmetry considers all of them together; mod is the least number of generations for the pattern to reappear, possibly under a transformation. There are 43 possibilities for the first two elements together (see table below).
`pops` | List of successive populations of the object, of length `p`. The presence of this field allows more complicated statistics to be performed on the population (standard deviation, etc.)
`popstats` | A quadruple (minimum population, maximum population, average population, sum of populations across a period). The last element does not convey any further _substantial_ information – the average population is a standard floating-point number, so the last element allows more accurate computation of the average population should such accuracy be desired.
`heats` | List of successive heats (births + deaths) of the object, of length `p`.
`heatstats` | A quadruple (minimum heat, maximum heat, average heat, sum of heats across a period). See the descriptions of `pops` and `popstats` for justification of the seemingly extraneous data.
`bb` | Measurements of the bounding box: (width, height, cell count [width × height]).
`rules` | A quadruple (isotropic minrule, isotropic maxrule, outer-totalistic (OT) minrule, OT maxrule). The last two entries are `None` if the pattern does not work in any OT rule.

The 43 possible combinations of base and temporal symmetries were enumerated in [Dean Hickerson's oscillator stamp collection](http://radicaleye.com/DRH/stamps.html):

Base symmetry (mod = period) | Half-period temporal symmetries | Quarter-period
---------------------------- | --------------------------------| --------------
`C1` | `C2_1`, `C2_2`, `C2_4`, `D2_+1`, `D2_+2`, `D2_x` | `C4_1`, `C4_4`
`C2_1` | `C4_1`, `D4_+1`, `D4_x1`
`C2_2` | `D4_+2`
`C2_4` | `C4_4`, `D4_+4`, `D4_x4`
`C4_1` | `D8_1`
`C4_4` | `D8_4`
`D2_+1` | `D4_+1`, `D4_+2`
`D2_+2` | `D4_+2`, `D4_+4`
`D2_x` | `D4_x1`, `D4_x4`
`D4_+1` | `D8_1`
`D4_+2` |
`D4_+4` | `D8_4`
`D4_x1` | `D8_1`
`D4_x4` | `D8_4`
`D8_1` |
`D8_4` |

The following fields are only computed for oscillators (and still lifes, which are p1 oscillators):

Key | Value
--- | -----
`cellperiods` | A dictionary mapping active cells (those positions that are ever alive) to their subperiods. Cells which are alive in the provided pattern have their subperiods _negated_; this subtlety is used to mark alive cells for the image.
`cellpcounts` | A dictionary mapping occurring cell subperiods to their corresponding frequencies.
`cellstats` | A quadruple (rotor cells, stator cells, active cells, strict (full-period) rotor cells). The stator consists of those cells which are always alive, with the remaining active cells forming the rotor.
`volstats` | A pair of floats (volatility, strict volatility) representing the proportion of (strict) rotor cells to active cells. The latter element is not provided by LifeViewer, which was one motivation for writing Nakano.
`tempstats` | A pair of floats (temperature, rotor temperature) representing the proportion of heat to active (rotor) cells. These two elements are also not provided by LifeViewer.

----

For auditing purposes the script that generated the `mooreconfs` array, essentially a lookup table mapping neighbourhood configurations to Hensel codes, is included as `mooreconfgen.py`. The function `minmaxrule()`, accepting a lifelib Pattern and a number of generations, is not limited to periodic objects and may be used for any pattern whatsoever – the only condition is that the provided rule and number of generations be of interest.
