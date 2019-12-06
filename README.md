# Chemical-Equilibrium with Julia
Gas phase chemical equilibrium calculation.
Function MicPV calculates chemical equilibrium of the gas mixture provides the temperature T and pressure p or volume v are specified.
Data input:
m-number of chemical elements;
k - number of substances;
b-amounts of chemical elements;
Nji-atomic matrix, number of atoms of element j in substance i;
F- Gibbs energy G divided by (RT), R - gas constant, F=-G/RT, G=H-TS, H - enthalpy, S - entropy.

Computed values:
equilibrium composition, mole;
volume, if pressure is specified,
pressure, if volume is specified.
