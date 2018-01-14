
/*
 * Oscillatory p53 system
 * (see Geva‐Zatorsky, N., Rosenfeld, N., Itzkovitz, S., 
 * Milo, R., Sigal, A., Dekel, E., Yarnitzky, T., Liron, Y., 
 * Polak, P., Lahav, G. and Alon, U., 2006.  Oscillations and 
 * variability in the p53 system. Molecular systems biology, 2(1).)
 */

const
  k1 = 90.0,
  k2 = 0.002,
  k3 = 1.7,
  k4 = 1.1,
  k5 = 0.93,
  k6 = 0.96,
  k7 = 0.01

var
  p53, preMdm2, Mdm2

chemical_reactions
  0 <-> p53  @ mass_action(k1), mass_action(k2)
  p53 -> 0  @ k3 * Mdm2 * p53 / (p53 + k7)
  0 -> preMdm2  @ k4 * p53
  preMdm2 -> Mdm2  @ mass_action(k5)
  Mdm2 -> 0  @ mass_action(k6)
end

init
  p53 = 70, preMdm2 = 30, Mdm2 = 60
end
