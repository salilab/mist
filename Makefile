PREFIX=/usr/local
BINDIR=${PREFIX}/bin

.PHONY: install test

FILES=${BINDIR}/MiST.py

install: ${FILES}

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@

test:
	cd test && nosetests
