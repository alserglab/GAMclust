# Calculate cosine similarity between rows of `a` and `b` (`b` could be a vector)

Calculate cosine similarity between rows of `a` and `b` (`b` could be a
vector)

## Usage

``` r
cosine(a, b = NULL)
```

## Arguments

- a:

  matrix

- b:

  matrix or vector, if null, then `b` is set to be equal to `a`

## Value

matrix of cosine distances between rows of `a` and rows of `b`
