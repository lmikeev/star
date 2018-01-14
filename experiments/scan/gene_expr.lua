
scan({
  model = 'models/gene_expr.star',
  objf = 'R=10 and D_on=1',
  tspan = {0, 300},
  nrepeat = 1000,
  param = {
    c1 = uniform(0.0010,0.10),
    c2 = uniform(0.010,1.0),
  },
  plot_objf_2d = {param1 = 'c1', param2 = 'c2'},
  stats = true
})
