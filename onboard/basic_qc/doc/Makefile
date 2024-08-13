PREFIX=code_project
BIBPREFIX=refs_si

$(PREFIX).pdf: $(PREFIX).tex
	mkdir -p out
	pdflatex --output-directory=out $(PREFIX).tex
	preview out/$(PREFIX).pdf

pdf: $(PREFIX).tex
	mkdir -p out
	pdflatex --output-directory=out $(PREFIX).tex
	pdflatex --output-directory=out $(PREFIX).tex
	preview out/$(PREFIX).pdf

all: $(PREFIX).tex
	mkdir -p out
	pdflatex --output-directory=out $(PREFIX).tex
	cp $(BIBPREFIX).bib out; cd out; bibtex $(PREFIX).aux; rm $(BIBPREFIX).bib; cd ..
	pdflatex --output-directory=out $(PREFIX).tex
	pdflatex --output-directory=out $(PREFIX).tex
	preview out/$(PREFIX).pdf

all_sha: $(PREFIX).tex
	pdflatex $(PREFIX).tex
	biber $(PREFIX)
	pdflatex $(PREFIX).tex
	pdflatex $(PREFIX).tex


view: out/$(PREFIX).pdf
	preview out/$(PREFIX).pdf

clean:
	rm -rf out
