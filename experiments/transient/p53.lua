
set_model('models/p53.star')
set_tspan(linspace(0, 40, 100))

set_det()

set('plot', true)
set('stats', true)

transient({nmoments = 2})
transient({nmoments = 4})
transient({nmoments = 6})

transient({method = 'ssa', ntrajectories = 1000})
