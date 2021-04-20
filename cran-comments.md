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
