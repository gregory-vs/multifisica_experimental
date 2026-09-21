# Relatório TC03

O relatório de capacitor está em `main.tex`, na raiz de `tc03`.

Para compilar, execute nesta pasta:

```bash
make
```

O PDF atualizado é `build/main.pdf`. São necessários `make`, `latexmk` e
TeX Live com os pacotes usados no documento.

Use `make watch` para recompilar ao salvar e `make clean` para remover
arquivos auxiliares preservando o PDF. No VS Code, abra `tc03` e use a
receita `latexmk` do LaTeX Workshop.

O estilo de `listings` inclui o tratamento dos acentos em português para
compilar os comentários e textos do código Python com pdfLaTeX.
