# gradepath plotting internals

Import anchor for the gradepath figure layer. The plotting code uses
many ggplot2 verbs, so the whole namespace is imported; data-frame
columns are referenced through the `rlang` `.data` pronoun inside
`aes()` to avoid R CMD check "no visible binding for global variable"
notes.
