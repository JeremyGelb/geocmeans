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
