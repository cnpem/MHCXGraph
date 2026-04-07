SPHINXBUILD = sphinx-build
SOURCEDIR   = source
BUILDDIR    = build

.PHONY: help html clean github

SPHINXBUILD = sphinx-build
SOURCEDIR   = source
BUILDDIR    = build

.PHONY: help html clean clean-root github

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

clean:
	rm -rf "$(BUILDDIR)"

clean-root:
	for f in * .*; do \
		case "$$f" in \
			. .. source Makefile .git .gitignore) ;; \
			*) rm -rf "$$f" ;; \
		esac \
	done

github: clean html
	$(MAKE) clean-root
	cp -a "$(BUILDDIR)/html/." .
	touch .nojekyll
