# viztest returns a well-formed viztest object

    
    Correspondents of PW Tests with CI Tests
      level psame pdiff      easy  method
    1  0.99     1   0.5 -2.050901  Lowest
    2  0.99     1   0.5 -2.050901  Middle
    3  0.99     1   0.5 -2.050901 Highest
    4  0.99     1   0.5 -2.050901 Easiest
    
    Missed Tests (n=2 of 10)
      bigger smaller pw_test ci_olap
    6   zero    cyl6   Insig      No
    7   zero    cyl8   Insig      No

---

    # A tibble: 5 x 6
        vij  s_zb min_zb max_zb    zb    ci
      <int> <dbl>  <dbl>  <dbl> <dbl> <dbl>
    1     1 0.506  0.945   1.92  1.44 0.850
    2     2 0.381  0.945   1.83  1.45 0.853
    3     3 0.325  1.07    1.79  1.49 0.863
    4     4 0.412  0.907   1.83  1.51 0.870
    5     5 0.474  0.907   1.92  1.61 0.893

