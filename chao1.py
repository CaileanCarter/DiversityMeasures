# From pycogent: https://github.com/pycogent/pycogent/tree/f720cc3753429d130f9e9bc0756b8878c3d50ef2

def osd(counts):
    """Returns observed, singles and doubles from counts.

    Handy for diversity calculations."""
    return (counts!=0).sum(), (counts==1).sum(), (counts==2).sum()

def chao1_bias_corrected(observed, singles, doubles):
    """Calculates bias-corrected chao1 given counts: Eq. 2 in EstimateS manual.

    Formula: chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1)) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the bias-corrected formulat from Chao 1987, Eq. 2 in the
    EstimateS manual.
    """
    return observed + singles*(singles-1) / (2.0*(doubles+1))

def chao1(counts):
    """Calculates chao1 according to table in EstimateS manual.

    Specifically, uses bias-corrected version unless bias_corrected is set
    to False _and_ there are both singletons and doubletons."""
    o, s, d = osd(counts)
    return chao1_bias_corrected(o, s, d)

