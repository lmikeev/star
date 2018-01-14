
transient({
  model = 'models/three_stage.star',
  tspan = linspace(0, 10),
  plot = {
    {'D_on', 'D_off'},
    {'R', 'P'},
  },
  plot_distr = {
    {'R', 'P'},
    {vars = {'R', 'P'}, cond = 'D_on=1'},
    {vars = {'R', 'P'}, cond = 'D_off=1'},
  },
  plot_distr_2d = {
    {var1 = 'R', var2 = 'P', cond = 'D_on=1'},
    {var1 = 'R', var2 = 'P', cond = 'D_off=1'},
  }
})
