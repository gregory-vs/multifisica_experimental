# Template LaTeX

Projeto base para reutilizar na criacao de novos documentos, seguindo a mesma
estrutura usada em `repos/latex/latex`.

## Estrutura

```text
template-latex/
├── main.tex
├── Makefile
├── latexmkrc
├── referencias.bib
├── capa/
│   └── capa.tex
├── relatorio/
│   └── conteudo.tex
├── imagens/
│   └── logos/
├── gerais/
│   ├── pacotes.tex
│   └── comandos.tex
└── build/
```

## Como reutilizar

1. Copie a pasta `template-latex/` para o nome do novo projeto.
2. Edite `gerais/comandos.tex` com titulo, disciplina, autores e data.
3. Ajuste `capa/capa.tex` se quiser trocar o layout da capa.
4. Escreva o texto principal em `relatorio/conteudo.tex`.
5. Coloque imagens em `imagens/` e logos em `imagens/logos/`.
6. Atualize `referencias.bib` se houver bibliografia.

## Compilacao

```bash
make -C tc01/latex
```

Execute esse comando na raiz do repositorio. Dentro de `tc01/latex`, basta
executar `make`. Tambem e possivel usar `make -f tc01/latex/Makefile` na raiz.
O Makefile sempre compila a partir da pasta do relatorio.

O PDF sera gerado em `tc01/latex/build/main.pdf` (o `main.pdf` fora de `build`
nao e atualizado).

Para compilar diretamente da raiz, use:

```bash
latexmk -cd -pdf -shell-escape -outdir=build tc01/latex/main.tex
```

A opcao `-cd` evita o erro `File gerais/pacotes.tex not found` ao compilar
fora da pasta do documento. No VS Code, use a receita `latexmk (shell-escape)`
do LaTeX Workshop. Os arquivos incluidos indicam `main.tex` como documento raiz.

Dependencias: `make`, `latexmk`, uma instalacao TeX Live com os pacotes usados
no preambulo e `pygmentize` (Pygments, exigido pelo pacote `minted` instalado).

## Compilacao continua

```bash
make -C tc01/latex watch
```

Esses comandos usam `latexmk` com `-shell-escape`, igual ao projeto base.

## Limpeza

```bash
make -C tc01/latex clean
```

A limpeza remove auxiliares e preserva o PDF gerado.
