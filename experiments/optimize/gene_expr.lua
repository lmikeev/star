
optimize({
  model = 'models/gene_expr.star',
  objf = 'R=10 and D_on=1',
  tspan = {0, 300},
  param = {
    c1 = {
      min = 0.0010,
      max = 0.10,
      init = 0.050
    },
    c2 = {
      min = 0.10,
      max = 1.0,
      init = 0.50
    }
  },
  maximize = true,
  opt_method = 'nlopt-ld-mma',
  plot_objf_2d = {param1='c1', param2='c2'},
  stats = true
})
