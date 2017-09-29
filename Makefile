.PHONY: build time_python time_cython

build:
	python3 setup.py build_ext --inplace

time_python:
	time python3 sim.py

time_cython:
	time python3 sim_cython.py
