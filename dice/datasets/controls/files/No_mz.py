#
# mz to No for NoDB
#

import math


def l10(x):
    return math.log10(x)


def e10(x):
    return 10.0 ** x


min = 300.00
max = 2000.00

lmin = l10(min)
lmax = l10(max)

ppm = 1.0e-6

hw = l10(1.0 + 20.0 * ppm)

N = int((lmax - lmin) / hw + 0.5)

# N can be used in range(1,N) calls since the actual range includes N-1 but not N.
# Does not mean that there are N bins -- there are only N-1 bins since NoDB always starts at 1!

domain = (min, max)
# range is meant to be inclusive, so not exactly Pythonic.
range = (1, N - 1)


def mz(numero):
    return [e10(lmin + hw * (numero - 1)), e10(lmin + hw * (numero + 1))]


def numero(mz):
    out = numeros(mz)
    if out:
        return out[0]
    else:
        return 0


def numeros(mz):
    lmz = l10(mz)
    if mz < min or mz > max:
        return []
    y = int(round((lmz - lmin) / hw))
    if lmz - (lmin + y * hw) >= 0.0:
        d = 1
    else:
        d = -1
    if y == 0:
        return [1]
    if y == N:
        return [N - 1]
    return [y, y + d]
