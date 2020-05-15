//! Calculates the Delauney Triangulation of a set of points
//! 
//! # Examples
//! By default it expects points in the range of 0.0 - 1.0
//! ```
//! use triangular::Triangulation;
//! 
//! // build a set of points
//! let points = vec![(0.3, 0.5),
//!                   (0.1, 0.6),
//!                   (0.4, 0.8)];
//! 
//! let triangles = Triangulation::new(&points).build::<usize>();
//! assert_eq!(&*triangles, &[2, 1, 0]);
//! ```

extern crate cgmath;

use cgmath::prelude::*;
use cgmath::BaseFloat;
use cgmath::Vector2;

use smallvec::{smallvec, SmallVec};

#[derive(Debug, Copy, Clone)]
enum Edge {
    BTOA = 0,
    CTOB = 1,
    ATOC = 2,
}

#[derive(Copy, Clone)]
enum TIndex<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>>,
{
    Point(usize, &'a [V]),
    Bound(usize, &'a [Vector2<T>]),
}

impl<T, V> TIndex<'_, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    fn value(&self) -> Vector2<T> {
        match self {
            Self::Point(i, source) => source[*i].into(),
            Self::Bound(i, source) => source[*i],
        }
    }
}

impl<T, V> TIndex<'_, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    fn get_index(&self) -> usize {
        match self {
            Self::Point(i, _) => *i,
            _ => panic!(),
        }
    }
}

impl<T, V> std::fmt::Debug for TIndex<'_, T, V>
where
    T: BaseFloat + std::fmt::Debug,
    V: Into<Vector2<T>> + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", *self)?;
        Ok(())
    }
}

#[derive(Debug)]
struct TriangleStore<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    store: Vec<LinkedTriangle<'a, T, V>>,
}

impl<'a, 'b, T, V> TriangleStore<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    fn new(size: usize) -> Self {
        Self {
            store: Vec::with_capacity(size),
        }
    }

    fn make(
        &mut self,
        a: TIndex<'a, T, V>,
        b: TIndex<'a, T, V>,
        c: TIndex<'a, T, V>,
    ) -> TriangleRef {
        let t = LinkedTriangle::new(a, b, c);
        self.store.push(t);
        TriangleRef::new(self.store.len() - 1)
    }

    fn get(&self, t: &TriangleRef) -> &LinkedTriangle<'a, T, V> {
        &self.store[t.t]
    }

    fn get_mut(&mut self, t: &TriangleRef) -> &mut LinkedTriangle<'a, T, V> {
        &mut self.store[t.t]
    }

    fn clear(&mut self, t: &TriangleRef) {
        self.get_mut(t).clear();
    }
}

#[derive(Debug)]
struct LinkedTriangle<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>>,
{
    a: TIndex<'a, T, V>,
    b: TIndex<'a, T, V>,
    c: TIndex<'a, T, V>,
    deleted: bool,
    placed: SmallVec<[TriangleRef; 8]>,
    edges: [Option<(TriangleRef, Edge)>; 3],
}

impl<'a, T, V> LinkedTriangle<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    fn new(
        a: TIndex<'a, T, V>,
        b: TIndex<'a, T, V>,
        c: TIndex<'a, T, V>,
    ) -> LinkedTriangle<'a, T, V> {
        LinkedTriangle {
            a,
            b,
            c,
            deleted: false,
            placed: SmallVec::new(),
            edges: [None; 3],
        }
    }

    fn get(&self, edge: Edge) -> &Option<(TriangleRef, Edge)> {
        &self.edges[edge as usize]
    }

    fn get_mut(&mut self, edge: Edge) -> &mut Option<(TriangleRef, Edge)> {
        &mut self.edges[edge as usize]
    }

    fn clear(&mut self) {
        self.deleted = true;
        self.edges = [None; 3];
    }

    fn indices<I>(&self) -> SmallVec<[I; 3]>
    where
        I: From<usize>,
    {
        smallvec![
            self.a.get_index().into(),
            self.b.get_index().into(),
            self.c.get_index().into()
        ]
    }

    #[inline(always)]
    fn max_lsq(&self) -> T {
        self.a.value().distance2(self.b.value()).max(
            self.b
                .value()
                .distance2(self.c.value())
                .max(self.c.value().distance2(self.a.value())),
        )
    }

    #[inline(always)]
    fn no_bounds(&self) -> bool {
        !matches!(self.a, TIndex::Bound(_, _))
            && !matches!(self.b, TIndex::Bound(_, _))
            && !matches!(self.c, TIndex::Bound(_, _))
    }

    #[inline(always)]
    fn in_circumcircle(&self, d: TIndex<T, V>) -> bool {
        math::in_circumcircle(self.a.value(), self.b.value(), self.c.value(), d.value())
    }

    #[inline(always)]
    fn contains(&self, d: TIndex<T, V>) -> bool {
        math::contains(self.a.value(), self.b.value(), self.c.value(), d.value())
    }
}

#[derive(Debug, Copy, Clone)]
struct TriangleRef {
    t: usize,
}

impl TriangleRef {
    fn new(t: usize) -> Self {
        Self { t }
    }

    fn triangle_place<'a, T, V>(
        &self,
        other: &TriangleRef,
        conn: Edge,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        let t = *store.get(self).get(conn);
        if let Some((t, edge)) = t {
            store.get_mut(&other).get_mut(Edge::CTOB).replace((t, edge));
            store
                .get_mut(&t)
                .get_mut(edge)
                .replace((*other, Edge::CTOB));
            other.edge_fix(&t, self, Edge::CTOB, edge, store);
        } else {
            store.get_mut(self).placed.push(*other);
        };
    }

    #[inline(always)]
    fn triangle_edge_build<'a, T, V>(
        &self,
        other: &Self,
        edge: Edge,
        conn: Edge,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        store.get_mut(self).get_mut(edge).replace((*other, conn));
    }

    fn triangle_triple_build<'a, T, V>(
        t1: &Self,
        t2: &Self,
        t3: &Self,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        t1.triangle_edge_build(t3, Edge::BTOA, Edge::ATOC, store);
        t1.triangle_edge_build(t2, Edge::ATOC, Edge::BTOA, store);
        t2.triangle_edge_build(t1, Edge::BTOA, Edge::ATOC, store);
        t2.triangle_edge_build(t3, Edge::ATOC, Edge::BTOA, store);
        t3.triangle_edge_build(t2, Edge::BTOA, Edge::ATOC, store);
        t3.triangle_edge_build(t1, Edge::ATOC, Edge::BTOA, store);
    }

    fn place<'a, T, V>(&self, store: &mut TriangleStore<'a, T, V>, p: TIndex<'a, T, V>)
    where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        let (t1, t2, t3) = {
            let LinkedTriangle { a, b, c, .. } = *store.get(self);
            (
                store.make(p, b, c),
                store.make(p, c, a),
                store.make(p, a, b),
            )
        };

        Self::triangle_triple_build(&t1, &t2, &t3, store);

        self.triangle_place(&t1, Edge::CTOB, store);
        self.triangle_place(&t2, Edge::ATOC, store);
        self.triangle_place(&t3, Edge::BTOA, store);

        store.clear(self);
    }

    fn build_edge_fix<'a, T, V>(
        &self,
        other: &TriangleRef,
        edge: Edge,
        conn: Edge,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        let t = *store.get(self).get(edge);
        if let Some((t, edge)) = t {
            store.get_mut(other).get_mut(conn).replace((t, edge));
            store.get_mut(&t).get_mut(edge).replace((*other, conn));
        }
    }

    fn next_edge_fix<'a, T, V>(
        &self,
        parent: &TriangleRef,
        edge: Edge,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        if let Some((t, conn)) = *store.get(self).get(edge) {
            self.edge_fix(&t, parent, edge, conn, store);
        } else if !store.get(parent).deleted {
            store.get_mut(parent).placed.push(*self);
        }
    }

    fn edge_fix<'a, T, V>(
        &self,
        other: &Self,
        parent: &Self,
        edge: Edge,
        conn: Edge,
        store: &mut TriangleStore<'a, T, V>,
    ) where
        T: BaseFloat,
        V: Into<Vector2<T>> + Copy,
    {
        let (a, b, c, bc, ca) = {
            let t = store.get(self);
            match edge {
                Edge::BTOA => (t.a, t.b, t.c, Edge::CTOB, Edge::ATOC),
                Edge::CTOB => (t.b, t.c, t.a, Edge::ATOC, Edge::BTOA),
                Edge::ATOC => (t.c, t.a, t.b, Edge::BTOA, Edge::CTOB),
            }
        };

        let (d, o1, o2) = {
            let t = store.get(other);
            match conn {
                Edge::BTOA => (t.c, Edge::CTOB, Edge::ATOC),
                Edge::CTOB => (t.a, Edge::ATOC, Edge::BTOA),
                Edge::ATOC => (t.b, Edge::BTOA, Edge::CTOB),
            }
        };

        if store.get(self).deleted || store.get(other).deleted {
            return;
        }

        if !store.get(self).in_circumcircle(d)
            && !matches!(d, TIndex::Bound(_, _))
            && !matches!(c, TIndex::Bound(_, _))
        {
            let (t1, t2) = { (store.make(a, d, c), store.make(d, b, c)) };

            t1.triangle_edge_build(&t2, Edge::CTOB, Edge::ATOC, store);
            t2.triangle_edge_build(&t1, Edge::ATOC, Edge::CTOB, store);

            store.get_mut(parent).placed.extend_from_slice(&[t1, t2]);

            self.build_edge_fix(&t1, ca, Edge::ATOC, store);
            self.build_edge_fix(&t2, bc, Edge::CTOB, store);
            other.build_edge_fix(&t1, o1, Edge::BTOA, store);
            other.build_edge_fix(&t2, o2, Edge::BTOA, store);

            t1.next_edge_fix(other, Edge::BTOA, store);
            t1.next_edge_fix(self, Edge::ATOC, store);
            t2.next_edge_fix(other, Edge::CTOB, store);
            t2.next_edge_fix(self, Edge::BTOA, store);

            store.clear(self);
            store.clear(other);
        } else {
            store.get_mut(parent).placed.push(*self);
        }
    }
}

pub struct Triangulation<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>>,
{
    points: &'a [V],
    bounds: Option<(V, V)>,
    max_len: Option<T>,
}

impl<'a, T, V> Triangulation<'a, T, V>
where
    T: BaseFloat,
    V: Into<Vector2<T>> + Copy,
{
    /// Returns a new Triangulations object that holds the configuration and a reference to the points
    /// ```
    /// use triangular::Triangulation;
    ///
    /// let points = vec![(0.6, 1.1),
    ///                   (0.2, 1.5),
    ///                   (0.8, 1.6)];
    ///
    /// let triangles = Triangulation::new(&points);
    /// ```
    pub fn new(points: &'a [V]) -> Self {
        Self {
            points,
            bounds: None,
            max_len: None,
        }
    }

    /// Return all Triangles calculated within the specified bounds
    /// ```
    /// use triangular::Triangulation;
    ///
    /// // build a set of points
    /// let points = vec![(0.6, 1.1),
    ///                   (0.2, 1.5),
    ///                   (0.8, 1.6)];
    ///
    /// let triangles = Triangulation::new(&points)
    ///                     .with_bounds((0.0, 0.0), (2.0, 2.0))
    ///                     .build::<usize>();
    /// assert_eq!(&*triangles, &[2, 1, 0]);
    /// ```
    pub fn with_bounds(mut self, c1: V, c2: V) -> Self {
        self.bounds = Some((c1, c2));
        self
    }

    
    /// Return all Triangles Calculated with edges < `max_len`
    /// ```
    /// use triangular::Triangulation;
    ///
    /// // build a set of points
    /// let points = vec![(0.6, 1.1),
    ///                   (0.2, 1.5),
    ///                   (0.8, 1.6)];
    ///
    /// let triangles = Triangulation::new(&points)
    ///                     .with_max_len(0.1)
    ///                     .build::<usize>();
    /// assert_eq!(&*triangles, &[]);
    /// ```
    pub fn with_max_len(mut self, max_len: T) -> Self {
        self.max_len = Some(max_len);
        self
    }

    /// Builds the triangulations and returns all indices of the generated triangles
    /// ```
    /// use triangular::Triangulation;
    ///
    /// // build a set of points
    /// let points = vec![(0.6, 1.1),
    ///                   (0.2, 1.5),
    ///                   (0.8, 1.6)];
    ///
    /// let triangles = Triangulation::new(&points)
    ///                     .build::<usize>();
    /// ```
    pub fn build<I>(self) -> Box<[I]>
    where
        I: From<usize>,
    {
        let bounds: Vec<Vector2<T>> = if let Some(bounds) = self.bounds {
            vec![
                (bounds.0.into().x - bounds.1.into().x, bounds.0.into().y).into(),
                (bounds.1.into().x, bounds.0.into().y).into(),
                (
                    bounds.1.into().x,
                    T::from(2.0).unwrap() * (bounds.1.into().y - bounds.0.into().y),
                )
                    .into(),
            ]
        } else {
            vec![
                (-T::one(), T::zero()).into(),
                (T::one(), T::zero()).into(),
                (T::one(), T::one() + T::one()).into(),
            ]
        };

        let a = TIndex::Bound(0, &bounds);
        let b = TIndex::Bound(1, &bounds);
        let c = TIndex::Bound(2, &bounds);

        let mut store = TriangleStore::new(10 * self.points.len());
        store.make(a, b, c);

        let size = (self.points.len() as f64).log2() as usize * 3;
        'placing: for i in 0..self.points.len() {
            let v = TIndex::Point(i, self.points);
            let mut options = Vec::with_capacity(size);
            options.push(TriangleRef::new(0));
            while !options.is_empty() {
                let i = options.pop().unwrap();
                let t = store.get(&i);
                if t.contains(v) {
                    if !t.deleted {
                        i.place(&mut store, v);
                        continue 'placing;
                    } else {
                        options.extend_from_slice(&t.placed);
                    }
                }
            }
        }

        let output: Vec<I> = store
            .store
            .iter()
            .filter(|x| {
                !x.deleted
                    && x.no_bounds()
                    && if let Some(len) = self.max_len {
                        x.max_lsq() < len * len
                    } else {
                        true
                    }
            })
            .flat_map(|x| x.indices())
            .collect();
        output.into_boxed_slice()
    }
}

mod math {
    use cgmath::BaseFloat;
    use cgmath::Vector2;

    #[inline(always)]
    pub(crate) fn in_circumcircle<T: BaseFloat>(
        a: Vector2<T>,
        b: Vector2<T>,
        c: Vector2<T>,
        d: Vector2<T>,
    ) -> bool {
        let a = a - d;
        let b = b - d;
        let c = c - d;
        (a.x.powi(2) + a.y.powi(2)) * (b.x * c.y - c.x * b.y)
            - (b.x.powi(2) + b.y.powi(2)) * (a.x * c.y - c.x * a.y)
            + (c.x.powi(2) + c.y.powi(2)) * (a.x * b.y - b.x * a.y)
            <= T::zero()
    }

    #[inline(always)]
    pub(crate) fn d_area<T: BaseFloat>(a: Vector2<T>, b: Vector2<T>, c: Vector2<T>) -> T {
        (-b.y * c.x + a.y * (-b.x + c.x) + a.x * (b.y - c.y) + b.x * c.y).abs()
    }

    #[inline(always)]
    pub(crate) fn contains<T: BaseFloat>(
        a: Vector2<T>,
        b: Vector2<T>,
        c: Vector2<T>,
        d: Vector2<T>,
    ) -> bool {
        let area = d_area(a, b, c).recip();
        let s = area * (a.y * c.x - a.x * c.y + (c.y - a.y) * d.x + (a.x - c.x) * d.y);
        let t = area * (a.x * b.y - a.y * b.x + (a.y - b.y) * d.x + (b.x - a.x) * d.y);
        s > T::zero() && t > T::zero() && s + t < T::one()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn triangulate() {
        let mut rng = StdRng::from_seed([0; 32]);

        let points: Vec<(f64, f64)> = (0..2500)
            .into_iter()
            .map(|_| (rng.gen_range(0.0, 1.0), rng.gen_range(0.0, 1.0)))
            .collect();
        let trianglution = Triangulation::new(&points);
        let t = trianglution.with_bounds((0.0, 0.0), (1.0, 1.0)).build::<usize>();

        assert_eq!(t.len(), 14508);
    }
}
