"""Generate the two pySecDec packages used as an external check of the PRB
coefficient integrals. Run once; then `sh make_all.sh`; then `python3 run.py`.

Both integrals are written in the polynomial-power form that sector
decomposition accepts, with an analytic regulator `delta` on the endpoint
monomial. Nothing from src/, production/ or the notebooks is used.
"""
from pySecDec import make_package

# 1. Cross-sector finite part, PRB Eq. (A3cross-profile):
#    Fp int_0^inf dx x^-4 J_C(x^-4) with J_C(g) = 2 int_0^1 l dl (1+l g)^-1/2,
#    i.e. 2 int_0^1 l dl int_0^inf dx x^(-2+delta) (x^4+l)^(-1/2) at delta -> 0.
#    The x=0 power divergence 4/(3 eps) becomes an explicit pole structure that
#    must cancel (delta^-1 coefficient 0); the delta^0 term is the finite part.
#    Map x = t/(1-t):  2 l t^(-2+delta) (1-t)^(2-delta) (t^4 + l (1-t)^4)^(-1/2).
make_package(
    name='fp_cross',
    integration_variables=['t','l'],
    regulators=['delta'],
    requested_orders=[0],
    polynomials_to_decompose=['t**(-2+delta)','(1-t)**(2-delta)','(t**4 + l*(1-t)**4)**(-1/2)'],
    remainder_expression='2*l',
    decomposition_method='iterative_no_primary',
)

# 2. Euler integral of the Lindhard piece, PRB Appendix A:
#    M_E = (16/15) int_0^1 dw w^(-1/2) (1-w)^(3/4) (1+w)^(-13/4) (2 - w + 4w^2 + w^3).
make_package(
    name='me_euler',
    integration_variables=['w'],
    regulators=['delta'],
    requested_orders=[0],
    polynomials_to_decompose=['w**(-1/2+delta)','(1-w)**(3/4)','(1+w)**(-13/4)'],
    remainder_expression='2 - w + 4*w**2 + w**3',
    prefactor='16/15',
    decomposition_method='iterative_no_primary',
)
