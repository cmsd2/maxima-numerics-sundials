# numerics-sundials

[![Docs](https://img.shields.io/badge/docs-online-blue)](https://cmsd2.github.io/maxima-numerics-sundials/)

SUNDIALS CVODE integration for [maxima-numerics](https://github.com/cmsd2/maxima-numerics), providing production-grade ODE solving with event detection for Maxima.

## Install

Install locally during development:

```
mxpm install --path . --editable
```

Or copy-install:

```
mxpm install --path .
```

## Usage

```maxima
load("numerics-sundials");
numerics-sundials_hello();
```

## Documentation

Build documentation artifacts (`.info` and help index):

```
mxpm doc build
```

Live preview with mdBook:

```
mxpm doc serve
```

## License

MIT
