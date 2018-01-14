
transient({
  model = 'models/birth_death.star',
  method = 'rk45',
  tspan = linspace(0, 50),
  abs_tol = 1e-15,
  rel_tol = 1e-3,
  dump_moments = {
    nmoments = 3,
    cmp_path = 'experiments/transient/birth_death_exact'
  },
  dump_distr = {
    cmp_path = 'experiments/transient/birth_death_exact'
  },
  stats = true,
  plot = true,
})
