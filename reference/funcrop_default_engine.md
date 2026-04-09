# Get or Set the Default Estimation Engine

Query or set the default engine used by funcrop fitting functions. When
called without arguments, returns the current default. When called with
an engine name, sets it as the session default.

## Usage

``` r
funcrop_default_engine(engine = NULL)
```

## Arguments

- engine:

  Character string specifying the engine to use as default: `"asreml"`
  or `"bayesreml"`. If `NULL` (the default), returns the current default
  engine without changing it.

## Value

Character string – the current (or newly set) default engine name.
Returned invisibly when setting.

## Details

The default engine is stored in `options("funcrop.engine")`. If the
option is not set, auto-detection is used: ASReml-R if installed,
otherwise bayesreml. If neither is installed, an error is raised.

Setting an engine that is not installed raises an error immediately, so
downstream fitting functions can rely on the default being valid.

## Examples

``` r
# Query current default
if (FALSE) { # \dontrun{
funcrop_default_engine()

# Set bayesreml as default
funcrop_default_engine("bayesreml")
} # }
```
