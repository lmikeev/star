
/*
 * Three-stage gene expression model 
 * (see Shahrezaei, V. and Swain, P.S., 2008. Analytical 
 * distributions for stochastic gene expression. Proceedings 
 * of the National Academy of Sciences, 105(45), pp.17256-17261.)
 */

const tau_on = 0.05
const tau_off = 0.05
const k_r = 10.0
const gamma_r = 1.0
const k_p = 4.0
const gamma_p = 1.0
const tau_on_p = 0.015

var R, P : species
var D_off, D_on : boolean

chemical_reactions
  D_off <-> D_on  @ mass_action(tau_on), mass_action(tau_off)
  D_on -> D_on + R  @ mass_action(k_r)
  R -> 0  @ mass_action(gamma_r)
  R -> R + P  @ mass_action(k_p)
  P -> 0  @ mass_action(gamma_p)
  P + D_off -> P + D_on  @ mass_action(tau_on_p)
end

init
  D_off = 1
end
