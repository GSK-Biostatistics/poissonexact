# ep_network's verbose = TRUE prints per-param/stratum progress

    Code
      ep_network(cases ~ group, data = toy2, strata = "strat", exact_vars = "group",
      verbose = TRUE)
    Message
      
      SAS rule (2): EXACTMAX(no var) or no OFFSET/EXACTMAX =>
                      S_i = intercept sufficient statistic (within stratum).
      
      === EXACT for parameter column: group1 (network algorithm) ===
        stratum 1: n=2, M=1
        stratum 2: n=2, M=1
        after finishing stratum 1/2 (1): live states = 2
        after finishing stratum 2/2 (2): live states = 3
    Output
      
      Exact conditional Poisson (categorical + STRATA; network algorithm)
      -------------------------------------------------------------------
       Parameter Estimate StdError Alpha   LowerCL  UpperCL PValue SuffStat
          group1        0 1.414214  0.05 -4.363058 4.363058      1        1

# ep_network's print method identifies the network algorithm

    Code
      print(fit)
    Output
      
      Exact conditional Poisson (categorical + STRATA; network algorithm)
      -------------------------------------------------------------------
       Parameter Estimate StdError Alpha   LowerCL  UpperCL PValue SuffStat
          group1        0 1.414214  0.05 -4.363058 4.363058      1        1

# .network_one_parameter prints its parameter header when verbose

    Code
      .network_one_parameter(model, jcol, slevels_ord, condition_on = "model", alpha = 0.05,
        ci_type = "two-sided", verbose = TRUE)
    Message
      
      === EXACT for parameter column: group1 (network algorithm) ===
        stratum 1: n=2, M=1
        stratum 2: n=2, M=1
        after finishing stratum 1/2 (1): live states = 2
        after finishing stratum 2/2 (2): live states = 3
    Output
      $parm
      [1] "group1"
      
      $outdist_entry
        SuffStat prob LogWeight
      1        0 0.25 0.0000000
      2        1 0.50 0.6931472
      3        2 0.25 0.0000000
      
      $result_row
        Parameter Estimate StdError Alpha   LowerCL  UpperCL PValue SuffStat
      1    group1        0 1.414214  0.05 -4.363058 4.363058      1        1
      

