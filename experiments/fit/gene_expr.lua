
-- estimate the kinetic constants
-- of the gene_expression model
-- based on time series data

fit({
  model = 'models/gene_expr.star',
  time_series_data = 'experiments/fit/gene_expr_data',

  opt_method = 'nlopt-ld-mma',

  -- specify range for the unknown parameters
  param = {
    c1 = {
      min = 0.0010,
      max = 0.10
    },
    c2 = {
      min = 0.010,
      max = 1.0
    },
    c3 = {
      min = 0.010,
      max = 1.0
    }
  },

  init = {
    R = {
      min = 0,
      max = 5
    }
  },

  obs_err = {
    R = {
      min = 0.00010,
      max = 2.0
    }
  },

  plot_objf = {{'c1', 'c2'}, {'c3'}},

  plot_objf_2d = {param1 = 'c1', param2 = 'c2'},

  stats = true
})
