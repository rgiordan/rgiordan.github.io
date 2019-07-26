#!/bin/bash
#pandoc --bibliography publications.bib --csl mycsl.csl -o citations.html -s publications_bibtex.md
pandoc --bibliography publications.bib -o citations.html publications_bibtex.md
#pandoc --bibliography publications.bib -t markdown -o citations.md -s publications_bibtex.md
