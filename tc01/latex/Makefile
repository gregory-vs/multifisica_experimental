PROJECT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MAIN ?= main.tex

.PHONY: all pdf watch clean

all: pdf

pdf:
	cd "$(PROJECT_DIR)" && latexmk "$(MAIN)"

watch:
	cd "$(PROJECT_DIR)" && latexmk -pvc "$(MAIN)"

clean:
	cd "$(PROJECT_DIR)" && latexmk -c "$(MAIN)"
