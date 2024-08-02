# Pooling


## In progress

```python
assign_wells_multidim(n_compounds:int, n_dims:int)
```
Here, should we give an indication to the user on how to choose n_dims ? Maybe we could have already a default value which depends on n_compounds ?

```python
assign_wells_random(n_compounds:int, n_compounds_per_well:int, n_wells:int, guesses:int, differentiate:int)
```
Same question here for the additional parameters. Also differentiate is given a default value (1) in other functions, so I'd do the same here ?

My idea with this is that if we get default values for the method-specific parameters, then we can easily run a comparison function that gives the (default) results of all methods at once just based on n_compounds and differentiate as a general easy overview.
Placeholder for this: method_comparison

## Installation


## Usage


## License

[MIT](https://choosealicense.com/licenses/mit/)
