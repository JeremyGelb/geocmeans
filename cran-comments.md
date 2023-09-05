# Version 0.1.1

## Test environments

* local R installation, R 4.0.2
* GitHub Actions - (ubuntu-18.04): release, devel
* GitHub Actions - (ubuntu-20.04): release, devel
* GitHub Actions - (windows): release, devel
* GitHub Actions - (macOS): release

Note : not ran on macOS devel, the package sf required by spdep can not be found.

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTES.

```
checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Jeremy Gelb <jeremy.gelb@ucs.inrs.ca>'
```

# Round 1 (after automatic checks and human validation) the 19/04/2021

**problem**:

 Found the following (possibly) invalid file URIs:
     URI: CONTRIBUTING.md
       From: README.md
     URI: CONDUCT.md
       From: README.md
     URI: LICENSE.txt
       From: README.md

**correction** : 
The three local links have been replaced by fully specified URLs.

# Round 2 (after automatic checks and human validation) the 20/04/2021

**problem** : 

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

**correction**

The three main references are added in the DESCRIPTION file.

**Problem** : 

Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or
'F' as vector names.)

**correction**

All the 'T' and 'F' are replaced by TRUE and FALSE accordingly.

**problem** : 

Please add \\value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      barPlots.Rd: \\value
      eval_parameters.Rd: \\value
      evaluateMatrices.Rd: \\value
      kppCenters.Rd: \\value
      sanity_check.Rd: \\value
      spiderPlots.Rd: \\value
      violinPlots.Rd: \\value
      
**correction**

For each function above, a value tag is added (with @return in Roxygen)

**problem** : 

\\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \\dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace
\\dontrun{} with \\donttest{}.

**correction** : 

\\dontrun{} was used to avoid running the function selectParameters because it is quite long (>5 sec). As suggested, it is now replaced by \\donttest{}


# Version 0.2.0

## Test environments

* local R installation, R 4.0.5
* GitHub Actions - (windows): release
* GitHub Actions - (macOS): release
* rhub - (ubuntu): release

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTES.

# Round 1, After automatic check and contact with CRAN team (2021/08/22)

**problem** : 

```
Found the following assignments to the global environment:
File ‘geocmeans/R/shinyapp.R’:
  assign("shiny_data", shiny_data, .GlobalEnv)
```

suggestion from CRAN team:

"Please use, e.g., an environment in your package and not the GlobalEnv
to assign to."

**solution**: 

An environment called `geocmeans_env` is created in package the geocmeans. It could be used in the future for other purposes.


# Version 0.2.1

## Test environments

* local R installation (windows), R Under development 4.3
* GitHub Actions - (windows): release
* rhub - (ubuntu): release
* rhub - (macos): release

## R CMD check results
0 ERRORs | 0 WARNINGs | 1 NOTES.


# Version 0.3.3

## Test environments

* local R installation (windows), R Under development 4.3
* GitHub Actions - (windows): release
* rhub - (ubuntu): release
* rhub - (macos): release

── R CMD check results ────────────────────────────────────────────── geocmeans 0.3.3 ────
Duration: 6m 48.8s

❯ checking installed package size ... NOTE
    installed size is  7.8Mb
    sub-directories of 1Mb or more:
      doc       1.8Mb
      extdata   3.0Mb
      libs      1.4Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## round 1 after manual checking from CRAN

```
 Package CITATION file contains call(s) to old-style personList() or
 as.personList().  Please use c() on person objects instead.
 Package CITATION file contains call(s) to old-style citEntry() or
 citHeader()/citFooter().  Please use bibentry() instead, possibly with
 arguments 'header' and 'footer'.
```

**solution**: 

The CITATION file was rewritten as suggested.


# Version 0.3.4

## Test environments

* local R installation (windows), R Under development (unstable) (2023-08-21 r84998 ucrt) -- "Unsuffered Consequences"
* GitHub Actions - (windows): release
* rhub - (ubuntu): release
* rhub - (macos): release

── R CMD check results ───────────────────────────────────────────────────────────────────────────────────────────────── geocmeans 0.3.4 ────
Duration: 6m 47.4s

❯ checking installed package size ... NOTE
    installed size is  7.7Mb
    sub-directories of 1Mb or more:
      doc       1.7Mb
      extdata   3.0Mb
      libs      1.4Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Round 1, After automatic check

**problem** : 

```
* checking S3 generic/method consistency ... NOTE
Mismatches for methods registered for non-generic:
is:
  function(object, class2)
is.FCMres:
  function(x)
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
```

**solution**:

The is.FCMres method has now two parameters: object and class2. The documentation has also be changed.

**problem** : 

```
* checking examples ... [65s/27s] NOTE
Examples with CPU time > 2.5 times elapsed time
               user system elapsed ratio
spConsistency 1.925  0.056    0.58 3.416
```

**solution**:

The example has been modified to take less time (25 replications instead of 50).
