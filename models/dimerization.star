
/*
 * Simple dimerization reaction 
 * (see Wilkinson, D.J., 2011. Stochastic modelling
 * for systems biology. CRC press.)
 */

const
  c1 = 1.66e-3,
  c2 = 0.2

var
  P, P2 : species

chemical_reactions
  2 P -> P2  @ c1*P*(P-1)/2
  P2 -> 2 P  @ c2*P2
end

init
  P = 301, P2 = 0
end
