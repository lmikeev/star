
set_tspan(linspace(0, 50))

set('method', 'ssa')
set('ntrajectories', 1000)

set('plot', true)
set('stats', true)

transient({
  model_sbml = 'models/sbml/dsmts-001-01.xml'
})

transient({
  model_sbml = 'models/sbml/dsmts-002-01.xml'
})

transient({
  model_sbml = 'models/sbml/dsmts-003-01.xml'
})
