# print.ep_fit prints the stored algorithm label and a matching rule

    Code
      print(fake)
    Output
      
      Exact conditional Poisson regression (multivariate shift)
      ---------------------------------------------------------
       Parameter Estimate StdError Alpha LowerCL UpperCL PValue SuffStat
           group      0.5      0.1  0.05     0.3     0.7   0.01       12

# print.ep_fit falls back to a default label when none is set

    Code
      print(fake)
    Output
      
      Exact conditional Poisson regression
      ------------------------------------
       Parameter Estimate StdError Alpha LowerCL UpperCL PValue SuffStat
               x        0        0  0.05       0       0      1        0

# print.ep_fit returns its argument invisibly

    Code
      ret <- withVisible(print(fake))
    Output
      
      label
      -----
       Parameter Estimate StdError Alpha LowerCL UpperCL PValue SuffStat
               x        0        0  0.05       0       0      1        0

