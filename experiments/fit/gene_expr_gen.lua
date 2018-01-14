
-- generate sample data

transient({
  model = 'models/gene_expr.star',
  tspan = linspace(0, 300, 30),

  method = 'ssa',
  nrepeat = 1,
  ntrajectories = 5,
  dump_trajectories = true,

  observables = {'D_off', 'D_on', 'R'},
  obs_err = {
    R = 1.0
  }
})
