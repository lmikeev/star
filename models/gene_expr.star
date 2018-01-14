
/*
 * Simple gene expression model
 * (see Golding, I., Paulsson, J., Zawilski, S.M. and 
 * Cox, E.C., 2005. Real-time kinetics of gene activity in 
 * individual bacteria. Cell, 123(6), pp.1025-1036.)
 */

var
  R : species,
  D_on, D_off : boolean

// kinetic constants
const
  c1 = 0.0270,
  c2 = 0.1667,
  c3 = 0.40

chemical_reactions
  D_off <-> D_on  @ mass_action(c1), mass_action(c2)
  D_on -> D_on + R  @ mass_action(c3)
end

// initial conditions
init
  D_off = true
end
