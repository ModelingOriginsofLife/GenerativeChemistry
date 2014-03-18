GenerativeChemistry
===================

This is a well-mixed reactor model of a 'generative' chemistry on strings - that is, the reaction network is generated from stoichiometry and underlying rules (e.g. bond energies, local chemical structure). This type of chemistry effectively self-extends its own reaction network at need.

[Click
here!](http://modelingoriginsoflife.github.io/GenerativeChemistry/autocata/)

This particular version implements a chemistry with a weakened
'non-locality' and normalized chemical environments. This tends to
(for whatever reason) produce more polymerization reactions.

There's an odd effect of this however:

B=C is stable in solution, as is A=D, but a single molecule of A=D
will consume the entirety of a bath of B=C to produce a family of B-C
variant polymers. This is due to the B=C molecule condensing on either
side with existing B-C polymers as well as exchanging across the
double bond. This means that the B-C polymers both grow through
condensation and shrink through splitting, giving rise to some sort of
population of molecules that can 'eat' B-C.

This is more complex a phenomenon than something like fire (another
run-away reaction), but lacks the perfect preservation of identity of
something like a molecule that is (individually) autocatalytic, so you
have some sort of state between two 'types' of amplification.