#!/bin/bash

#mkdir pdffiles
#find . -maxdepth 1 -name '*.pdf' -print0 | xargs -0 -I {} "pdffiles/"
find . -maxdepth 1 -name '*.pdf' -type f | xargs -n 1 -I '{}' mv {}  pdffiles/
find . -maxdepth 1 -name '*.log' -type f | xargs -n 1 -I '{}' mv {}  logfiles/
find . -maxdepth 1 -name '*gene*' -type f | xargs -n 1 -I '{}' mv {}  genefiles/

