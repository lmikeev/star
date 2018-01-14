
/*
 * Network of two genes with overlapping 
 * promotor region and mutual repression
 * (see Barzel, B. and Biham, O., 2008. Calculation of 
 * switching times in the genetic toggle switch and 
 * other bistable systems. Physical Review E, 78(4), p.041919.)
 */

const
  k1 = 0.20,
  k2 = 0.50,
  d1 = 0.005,
  d2 = 0.005,
  b1 = 0.005,
  b2 = 0.002,
  u1 = 0.02,
  u2 = 0.02

var
  P1, P2 : species,
  dna, dnaP1, dnaP2 : boolean

chemical_reactions
  // production of P1
  dna -> dna + P1  @ mass_action(k1)

  // production of P2
  dna -> dna + P2  @ mass_action(k2)

  // degradation of P1
  P1 -> 0  @ mass_action(d1)

  // degradation of P2
  P2 -> 0  @ mass_action(d2)

  // binding/unbinding of P1
  dna + P1 <-> dnaP1  @ mass_action(b1), mass_action(u1)

  // binding/unbinding of P2
  dna + P2 <-> dnaP2  @ mass_action(b2), mass_action(u2)

  // production of P1
  dnaP1 -> dnaP1 + P1  @ mass_action(k1)

  // production of P2
  dnaP2 -> dnaP2 + P2  @ mass_action(k2)
end

// initial conditions
init
  dna = true
end
