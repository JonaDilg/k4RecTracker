MAKEFLAGS += --no-print-directory

.PHONY: make install clean

make:
	@ mkdir -p build install ; \
	cd build ; \
	cmake .. -DCMAKE_INSTALL_PREFIX=../install ; \
	make install -j8 ; \
	cd .. ; \
	export LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:$$LD_LIBRARY_PATH ; \
	export PYTHONPATH=${CURDIR}/install/python:$$PYTHONPATH ; \
	printf "#!/bin/bash\nif [ -n \"\$$KEY4HEP_STACK\" ];\nthen\n  echo '----> Info: Key4hep stack already set up. Skipping...'\nelse\n source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh\nfi\nexport LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:\$$LD_LIBRARY_PATH\nexport PYTHONPATH=${CURDIR}/install/python:\$$PYTHONPATH\n" > ${CURDIR}/setup.sh ; \
	chmod +x ${CURDIR}/setup.sh

install:
	@ mkdir -p build install ; \
	cd build ; \
	make install -j8 ; \
	cd .. ; \
	export LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:$$LD_LIBRARY_PATH ; \
	export PYTHONPATH=${CURDIR}/install/python:$$PYTHONPATH ; \
	printf "#!/bin/bash\nif [ -n \"\$$KEY4HEP_STACK\" ];\nthen\n  echo '----> Info: Key4hep stack already set up. Skipping...'\nelse\n source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh\nfi\nexport LD_LIBRARY_PATH=${CURDIR}/install/lib:${CURDIR}/install/lib64:\$$LD_LIBRARY_PATH\nexport PYTHONPATH=${CURDIR}/install/python:\$$PYTHONPATH\n" > ${CURDIR}/setup.sh ; \
	chmod +x ${CURDIR}/setup.sh

clean:
	@ (rm -r build install && rm setup.sh) || true
