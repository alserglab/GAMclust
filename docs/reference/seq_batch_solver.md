# Function for batch solving SGMWCS problems

Creates a wrapper for sequentially solving multiple SGMWCS optimization
instances. The returned function applies a selected MWCS solver to each
network independently and collects the results into a single list.

## Usage

``` r
seq_batch_solver(mwcs_solver)
```

## Arguments

- mwcs_solver:

  SGMWCS solver from mwcsr package
