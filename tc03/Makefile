PROJECT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

.PHONY: all pdf watch clean

all: pdf

pdf:
	cd "$(PROJECT_DIR)" && latexmk main.tex

watch:
	cd "$(PROJECT_DIR)" && latexmk -pvc main.tex

clean:
	cd "$(PROJECT_DIR)" && latexmk -c main.tex
