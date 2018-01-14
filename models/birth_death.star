
const c1 = 1.0
const c2 = 0.10

var A

chemical_reactions
  0 <-> A  @ mass_action(c1), mass_action(c2)
end

init
  A = 1000
end
