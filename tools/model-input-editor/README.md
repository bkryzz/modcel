# Model Input Editor

Small Tkinter utilities for preparing structured MODCEL input files.

- `input_editor.py` edits simple `VAR = value` input templates.
- `cell_editor.py` edits cell phenotypes and interaction matrices.

`input_editor.py` can call `cell_editor.py` as the next step in the workflow.
`cell_editor.py` can optionally call a local simulation launcher named
`play.py`. That launcher is not included in this public repository.

## Run

```bash
python3 input_editor.py
```

The included `default` file is a neutral toy template for testing the editor.
It is not a scientific MODCEL parameter set.

