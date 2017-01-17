# Eddy current simulations prototype

## Dependencies

- [Arb v2.9.0](http://arblib.org/)

- [GSL v2.1 or later](https://www.gnu.org/software/gsl/). On Ubuntu install with `apt-get install libgsl-dev`

## Compiling

```
g++ -Wall -c hypergeometric_pfq.cpp half-space.cpp && g++ -o half-space.out half-space.o hypergeometric_pfq.o -larb -lflint -lgsl -lgslcblas -lm
```
