
set_model('models/dimerization.star')
set_tspan(linspace(0, 10, 20))
set_hybrid()
set('nmoments', 3)
set('plot', true)
set('stats', true)

transient({
  hybrid_vars_force_stoch = {'P'},
  plot_distr = {vars = 'P'},
})

transient({
  hybrid_vars_force_stoch = {'P2'},
  plot_distr = {'P2'},
})
