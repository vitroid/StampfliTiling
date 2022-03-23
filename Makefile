.DELETE_ON_ERROR:

test: 0.cage 1.cage
	echo Done.
%.cage:
	genice2 stampfli[$*] -f cage[python] | awk '/^1/&&(NF==4){c[$$1]++}END{for(s in c){print s,c[s]}}'

#test4
#12 136
#14 16
#15 16
#16 56
pep8:
	autopep8 -r -a -a -i ./


GENICE=genice2
PACKAGE=genice2_stampfli


# %: temp_% replacer.py $(BASE)/lattices/cif.py $(BASE)/lattices/zeolite.py $(BASE)/__init__.py
# 	pip install genice2_dev
# 	python replacer.py < $< > $@


# prepare: # might require root privilege.
# 	pip install genice cif2ice validators


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install pillow
	pip install --index-url https://test.pypi.org/simple/ $(PACKAGE)



install:
	./setup.py install
uninstall:
	-pip uninstall -y $(PACKAGE)
build: README.md $(wildcard genice2_cif/lattices/*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload dist/*
check:
	./setup.py check
clean:
	-rm *~ *gro *cif
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
