# ORCA Python Library
This library uses `maturin` to manage the Rust library. 

To install first create a virtual environment 
`python -m venv orca-env`

To activate
Windows:
`orca-env\Scripts\activate`

macOS \& Unix:
`source tutorial-env/bin/activate`

Install maturin 
`pip install maturin`

Build the library
`maturin develop`

# Running Tests
From the top level `orca` directory run 

```
pip install pytest

pytest python/tests
```