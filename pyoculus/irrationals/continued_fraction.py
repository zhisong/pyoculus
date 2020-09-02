## @file continued_fraction.py
#  @brief Contains functions related to the continued fraction expansion
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#
import numpy as np


def expandcf(irrational, n, thres_ai=1000):
    """! Expand an (positive) irrational using the continued fraction expansion
    @param irrational -- the positive irrational number to expand. an absolute will be taken if negative
    @param n -- the number of terms
    @param thres_ai=1000 -- terminate the expansion if some ai>thres_ai

    @param an integer sequence contains the continued fraction expansion of irrational up to the nth term
    """

    ai = np.zeros(n, dtype=np.int)

    residue = np.abs(irrational)

    for ii in range(n):

        ai[ii] = np.floor(residue).astype(np.int)
        nterms = ii + 1

        if ai[ii] > thres_ai:
            # ai is beyond threshold, the expansion has terminated
            ai[ii] = 0
            nterms = ii
            break

        residue = residue - ai[ii]

        if np.abs(residue) < 1 / float(thres_ai):
            # residue is too small, the expansion has terminated
            break

        residue = 1.0 / residue

    ai = ai[0:nterms]

    return ai


def fromcf(ai):
    """! Obtain the fraction pp/qq of the continued fraction ai
    @param ai an integer array contains ai for the continued fraction expansion
    @returns (pp, qq), the fraction (pp/qq)
    """

    numer = ai[-1]
    denom = 1

    for jj in range(len(ai) - 1):
        tmpi = denom
        denom = numer
        numer = tmpi
        numer = numer + ai[-jj - 2] * denom

    pp = int(numer)
    qq = int(denom)

    return (pp, qq)
