GenerativeChemistry
===================

This is a well-mixed reactor model of a 'generative' chemistry on strings - that is, the reaction network is generated from stoichiometry and underlying rules (e.g. bond energies, local chemical structure). This type of chemistry effectively self-extends its own reaction network at need.

[Click here!](http://modelingoriginsoflife.github.io/GenerativeChemistry/)

Interesting reactions to explore:

A-B-C-B-A-A-D: Catalyzes A-A-D-D-C-B -> A-D + D-C-B-A

A=B-C-C: Makes catalyst to produce A-D-D-C-B from A-A-D-D-C-B

A=C: Condensation agent with D-C-B-A to form A-C and D-C-B-A-D-C-B-A-...

A=B: Condensation agent with D-C-B-A to form A-B and D-C-B-A-D-C-B-A

B=A-A=B: Condensation agent with D-C-B-A -> A-B-D-C-B-A + D-C-B-A-B-A-D-C-B-A + ...
