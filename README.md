# Eddy current simulations prototype

## Dependencies

- [Arb v2.9.0](http://arblib.org/)

- [GSL v2.1 or later](https://www.gnu.org/software/gsl/). On Ubuntu install with `apt-get install libgsl-dev`

## Compiling

```
g++ -Wall -c lib/hypergeometric_pfq.cpp lib/struve.cpp lib/theo.cpp main.cpp && g++ -o eddy-sim.out hypergeometric_pfq.o struve.o theo.o main.o -larb -lflint -lgsl -lgslcblas -lm
```
