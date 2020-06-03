"""Estimated Glomerular Filtration Rate (eGFR) calculator

Currently only the CKD-EPI equation is implemented.

"""
KAPPA_FEMALE = 0.7
KAPPA_MALE = 0.9
ALPHA_FEMALE = -0.329
ALPHA_MALE = -0.411
FEMALE = 1.018
BLACK = 1.159

def ckdepi(scr, age, is_female, is_black):
    """Calculate estimated GFR using the CKD-EPI equation

    Parameters
    ----------
    scr : float
        serum creatinine, mg/dL
    age : int
        age, years
    is_female : bool
    is_black : bool

    Returns
    -------
    float
        estimated GFR, mL/min per 1.73 m2

    Notes
    -----
    Implemented as described in [1]_ (table 2, footnotes).

    References
    ----------
    .. [1] Levey, AS, Stevens, LA, Schmid, CH, Zhang, YL, Castro, AF, Feldman,
       HI, Kusek, JW, Eggers, P, Van Lente, F, Greene, T, Coresh, J (2009).
       A new equation to estimate glomerular filtration rate. Ann. Intern.
       Med., 150, 9:604-12.

    """
    if is_female:
        kappa = KAPPA_FEMALE
        alpha = ALPHA_FEMALE
        f = FEMALE
    else:
        kappa = KAPPA_MALE
        alpha = ALPHA_MALE
        f = 1

    if is_black:
        b = BLACK
    else:
        b = 1

    x = scr/kappa
    return 141*min(x, 1)**alpha*max(x, 1)**(-1.209)*0.993**age*f*b
