# triangular
A Delauney triangulation library for rust

Calculates the Delauney Triangulation of a set of points

# Examples
By default it expects points in the range of 0.0 - 1.0
```rust
use triangular::Triangulation;

// build a set of points
let points = vec![(0.3, 0.5),
                  (0.1, 0.6),
                  (0.4, 0.8)];

let triangles = Triangulation::new(&points).build::<usize>();
assert_eq!(&*triangles, &[2, 1, 0]);
```
