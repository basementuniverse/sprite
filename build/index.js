(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else {
		var a = factory();
		for(var i in a) (typeof exports === 'object' ? exports : root)[i] = a[i];
	}
})(self, () => {
return /******/ (() => { // webpackBootstrap
/******/ 	var __webpack_modules__ = ({

/***/ "./node_modules/@basementuniverse/utils/utils.js":
/*!*******************************************************!*\
  !*** ./node_modules/@basementuniverse/utils/utils.js ***!
  \*******************************************************/
/***/ ((module) => {

/**
 * @overview A library of useful functions
 * @author Gordon Larrigan
 */

/**
 * Check if two numbers are approximately equal
 * @param {number} a Number a
 * @param {number} b Number b
 * @param {number} [p=Number.EPSILON] The precision value
 * @return {boolean} True if numbers a and b are approximately equal
 */
const floatEquals = (a, b, p = Number.EPSILON) => Math.abs(a - b) < p;

/**
 * Clamp a number between min and max
 * @param {number} a The number to clamp
 * @param {number} [min=0] The minimum value
 * @param {number} [max=1] The maximum value
 * @return {number} A clamped number
 */
const clamp = (a, min = 0, max = 1) => a < min ? min : (a > max ? max : a);

/**
 * Get the fractional part of a number
 * @param {number} a The number from which to get the fractional part
 * @return {number} The fractional part of the number
 */
const frac = a => a >= 0 ? a - Math.floor(a) : a - Math.ceil(a);

/**
 * Round n to d decimal places
 * @param {number} n The number to round
 * @param {number} [d=0] The number of decimal places to round to
 * @return {number} A rounded number
 */
const round = (n, d = 0) => {
  const p = Math.pow(10, d);
  return Math.round(n * p + Number.EPSILON) / p;
}

/**
 * Do a linear interpolation between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value, should be in the interval [0, 1]
 * @return {number} An interpolated value in the interval [a, b]
 */
const lerp = (a, b, i) => a + (b - a) * i;

/**
 * Get the position of i between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolated value in the interval [a, b]
 * @return {number} The position of i between a and b
 */
const unlerp = (a, b, i) => (i - a) / (b - a);

/**
 * Do a bilinear interpolation
 * @param {number} c00 Top-left value
 * @param {number} c10 Top-right value
 * @param {number} c01 Bottom-left value
 * @param {number} c11 Bottom-right value
 * @param {number} ix Interpolation value along x
 * @param {number} iy Interpolation value along y
 * @return {number} A bilinear interpolated value
 */
const blerp = (c00, c10, c01, c11, ix, iy) => lerp(lerp(c00, c10, ix), lerp(c01, c11, ix), iy);

/**
 * Re-map a number i from range a1...a2 to b1...b2
 * @param {number} i The number to re-map
 * @param {number} a1
 * @param {number} a2
 * @param {number} b1
 * @param {number} b2
 * @return {number}
 */
const remap = (i, a1, a2, b1, b2) => b1 + (i - a1) * (b2 - b1) / (a2 - a1);

/**
 * Do a smooth interpolation between a and b
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value
 * @return {number} An interpolated value in the interval [a, b]
 */
const smoothstep = (a, b, i) => lerp(a, b, 3 * Math.pow(i, 2) - 2 * Math.pow(i, 3));

/**
 * Get an angle in radians
 * @param {number} degrees The angle in degrees
 * @return {number} The angle in radians
 */
const radians = degrees => (Math.PI / 180) * degrees;

/**
 * Get an angle in degrees
 * @param {number} radians The angle in radians
 * @return {number} The angle in degrees
 */
const degrees = radians => (180 / Math.PI) * radians;

/**
 * Get a random float in the interval [min, max)
 * @param {number} min Inclusive min
 * @param {number} max Exclusive max
 * @return {number} A random float in the interval [min, max)
 */
const randomBetween = (min, max) => Math.random() * (max - min) + min;

/**
 * Get a random integer in the interval [min, max]
 * @param {number} min Inclusive min
 * @param {number} max Inclusive max
 * @return {number} A random integer in the interval [min, max]
 */
const randomIntBetween = (min, max) => Math.floor(Math.random() * (max - min + 1)) + min;

/**
 * Get a normally-distributed random number
 * @param {number} [mu=0.5] The mean value
 * @param {number} [sigma=0.5] The standard deviation
 * @param {number} [samples=2] The number of samples
 * @return {number} A normally-distributed random number
 */
const cltRandom = (mu = 0.5, sigma = 0.5, samples = 2) => {
  let total = 0;
  for (let i = samples; i--;) {
    total += Math.random();
  }
  return mu + (total - samples / 2) / (samples / 2) * sigma;
};

/**
 * Get a normally-distributed random integer in the interval [min, max]
 * @param {number} min Inclusive min
 * @param {number} max Inclusive max
 * @return {number} A normally-distributed random integer
 */
const cltRandomInt = (min, max) => Math.floor(min + cltRandom(0.5, 0.5, 2) * (max + 1 - min));

/**
 * Return a weighted random integer
 * @param {Array<number>} w An array of weights
 * @return {number} An index from w
 */
const weightedRandom = w => {
  let total = w.reduce((a, i) => a + i, 0), n = 0;
  const r = Math.random() * total;
  while (total > r) {
    total -= w[n++];
  }
  return n - 1;
};

/**
 * An interpolation function
 * @callback InterpolationFunction
 * @param {number} a The minimum number
 * @param {number} b The maximum number
 * @param {number} i The interpolation value, should be in the interval [0, 1]
 * @return {number} The interpolated value in the interval [a, b]
 */

/**
 * Return an interpolated value from an array
 * @param {Array<number>} a An array of values interpolate
 * @param {number} i A number in the interval [0, 1]
 * @param {InterpolationFunction} [f=Math.lerp] The interpolation function to use
 * @return {number} An interpolated value in the interval [min(a), max(a)]
 */
const lerpArray = (a, i, f = lerp) => {
  const s = i * (a.length - 1);
  const p = clamp(Math.trunc(s), 0, a.length - 1);
  return f(a[p] || 0, a[p + 1] || 0, frac(s));
};

/**
 * Get the dot product of two vectors
 * @param {Array<number>} a Vector a
 * @param {Array<number>} b Vector b
 * @return {number} a ∙ b
 */
const dot = (a, b) => a.reduce((n, v, i) => n + v * b[i], 0);

/**
 * Get the factorial of a number
 * @param {number} a
 * @return {number} a!
 */
const factorial = a => {
  let result = 1;
  for (let i = 2; i <= a; i++) {
    result *= i;
  }
  return result;
};

/**
 * Get the number of permutations of r elements from a set of n elements
 * @param {number} n
 * @param {number} r
 * @return {number} nPr
 */
const npr = (n, r) => factorial(n) / factorial(n - r);

/**
 * Get the number of combinations of r elements from a set of n elements
 * @param {number} n
 * @param {number} r
 * @return {number} nCr
 */
const ncr = (n, r) => factorial(n) / (factorial(r) * factorial(n - r));

/**
 * Generate all combinations of r elements from an array
 *
 * @example
 * ```js
 * combinations([1, 2, 3], 2);
 * ```
 *
 * Output:
 * ```json
 * [
 *   [1, 2],
 *   [1, 3],
 *   [2, 3]
 * ]
 * ```
 * @param {Array<*>} a
 * @param {number} r The number of elements to choose in each combination
 * @return {Array<Array<*>>} An array of combination arrays
 */
const combinations = (a, r) => {
  if (r === 1) {
    return a.map(item => [item]);
  }

  return a.reduce(
    (acc, item, i) => [
      ...acc,
      ...combinations(a.slice(i + 1), r - 1).map(c => [item, ...c]),
    ],
    []
  );
};

/**
 * Get a cartesian product of arrays
 *
 * @example
 * ```js
 * cartesian([1, 2, 3], ['a', 'b']);
 * ```
 *
 * Output:
 * ```json
 * [
 *   [1, "a"],
 *   [1, "b"],
 *   [2, "a"],
 *   [2, "b"],
 *   [3, "a"],
 *   [3, "b"]
 * ]
 * ```
 */
const cartesian = (...arr) =>
  arr.reduce(
    (a, b) => a.flatMap(c => b.map(d => [...c, d])),
    [[]]
  );

/**
 * A function for generating array values
 * @callback TimesFunction
 * @param {number} i The array index
 * @return {*} The array value
 */

/**
 * Return a new array with length n by calling function f(i) on each element
 * @param {TimesFunction} f
 * @param {number} n The size of the array
 * @return {Array<*>}
 */
const times = (f, n) => Array(n).fill(0).map((_, i) => f(i));

/**
 * Return an array containing numbers 0->(n - 1)
 * @param {number} n The size of the array
 * @return {Array<number>} An array of integers 0->(n - 1)
 */
const range = n => times(i => i, n);

/**
 * Zip 2 arrays together, i.e. ([1, 2, 3], [a, b, c]) => [[1, a], [2, b], [3, c]]
 * @param {Array<*>} a
 * @param {Array<*>} b
 * @return {Array<Array<*>>}
 */
const zip = (a, b) => a.map((k, i) => [k, b[i]]);

/**
 * Return array[i] with positive and negative wrapping
 * @param {Array<*>} a
 * @param {number} i The positively/negatively wrapped array index
 * @return {*} An element from the array
 */
const at = (a, i) => a[i < 0 ? a.length - (Math.abs(i + 1) % a.length) - 1 : i % a.length];

/**
 * Return the last element of an array without removing it
 * @param {Array<*>} a
 * @return {*} The last element from the array
 */
const peek = (a) => {
  if (!a.length) {
    return undefined;
  }

  return a[a.length - 1];
};

/**
 * Chop an array into chunks of size n
 * @param {Array<*>} a
 * @param {number} n The chunk size
 * @return {Array<Array<*>>} An array of array chunks
 */
const chunk = (a, n) => times(i => a.slice(i * n, i * n + n), Math.ceil(a.length / n));

/**
 * Randomly shuffle a shallow copy of an array
 * @param {Array<*>} a
 * @return {Array<*>} The shuffled array
 */
const shuffle = a => a.slice().sort(() => Math.random() - 0.5);

/**
 * Flatten an object
 * @param {object} o
 * @param {string} concatenator The string to use for concatenating keys
 * @return {object} A flattened object
 */
const flat = (o, concatenator = '.') => {
  return Object.keys(o).reduce((acc, key) => {
    if (o[key] instanceof Date) {
      return {
        ...acc,
        [key]: o[key].toISOString(),
      };
    }

    if (typeof o[key] !== 'object' || !o[key]) {
      return {
        ...acc,
        [key]: o[key],
      };
    }
    const flattened = flat(o[key], concatenator);

    return {
      ...acc,
      ...Object.keys(flattened).reduce(
        (childAcc, childKey) => ({
          ...childAcc,
          [`${key}${concatenator}${childKey}`]: flattened[childKey],
        }),
        {}
      ),
    };
  }, {});
};

/**
 * Unflatten an object
 * @param {object} o
 * @param {string} concatenator The string to check for in concatenated keys
 * @return {object} An un-flattened object
 */
const unflat = (o, concatenator = '.') => {
  let result = {}, temp, substrings, property, i;

  for (property in o) {
    substrings = property.split(concatenator);
    temp = result;
    for (i = 0; i < substrings.length - 1; i++) {
      if (!(substrings[i] in temp)) {
        if (isFinite(substrings[i + 1])) {
          temp[substrings[i]] = [];
        } else {
          temp[substrings[i]] = {};
        }
      }
      temp = temp[substrings[i]];
    }
    temp[substrings[substrings.length - 1]] = o[property];
  }

  return result;
};

/**
 * A split predicate
 * @callback SplitPredicate
 * @param {any} value The current value
 * @return {boolean} True if the array should split at this index
 */

/**
 * Split an array into sub-arrays based on a predicate
 * @param {Array<*>} array
 * @param {SplitPredicate} predicate
 * @return {Array<Array<*>>} An array of arrays
 */
const split = (array, predicate) => {
  const result = [];
  let current = [];
  for (const value of array) {
    if (predicate(value)) {
      if (current.length) {
        result.push(current);
      }
      current = [value];
    } else {
      current.push(value);
    }
  }
  result.push(current);

  return result;
};

/**
 * Pluck keys from an object
 * @param {object} o
 * @param {...string} keys The keys to pluck from the object
 * @return {object} An object containing the plucked keys
 */
const pluck = (o, ...keys) => {
  return keys.reduce(
    (result, key) => Object.assign(result, { [key]: o[key] }),
    {}
  );
};

/**
 * Exclude keys from an object
 * @param {object} o
 * @param {...string} keys The keys to exclude from the object
 * @return {object} An object containing all keys except excluded keys
 */
const exclude = (o, ...keys) => {
  return Object.fromEntries(
    Object.entries(o).filter(([key]) => !keys.includes(key))
  );
};

if (true) {
  module.exports = {
    floatEquals,
    clamp,
    frac,
    round,
    lerp,
    unlerp,
    blerp,
    remap,
    smoothstep,
    radians,
    degrees,
    randomBetween,
    randomIntBetween,
    cltRandom,
    cltRandomInt,
    weightedRandom,
    lerpArray,
    dot,
    factorial,
    npr,
    ncr,
    combinations,
    cartesian,
    times,
    range,
    zip,
    at,
    peek,
    chunk,
    shuffle,
    flat,
    unflat,
    split,
    pluck,
    exclude,
  };
}


/***/ }),

/***/ "./node_modules/@basementuniverse/vec/vec.js":
/*!***************************************************!*\
  !*** ./node_modules/@basementuniverse/vec/vec.js ***!
  \***************************************************/
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {

const { times, chunk, dot } = __webpack_require__(/*! @basementuniverse/utils */ "./node_modules/@basementuniverse/utils/utils.js");

/**
 * @overview A small vector and matrix library
 * @author Gordon Larrigan
 */

/**
 * A 2d vector
 * @typedef {Object} vec
 * @property {number} x The x component of the vector
 * @property {number} y The y component of the vector
 */

/**
 * Create a new vector
 * @param {number|vec} [x] The x component of the vector, or a vector to copy
 * @param {number} [y] The y component of the vector
 * @return {vec} A new vector
 * @example <caption>Various ways to initialise a vector</caption>
 * let a = vec(3, 2);  // (3, 2)
 * let b = vec(4);     // (4, 4)
 * let c = vec(a);     // (3, 2)
 * let d = vec();      // (0, 0)
 */
const vec = (x, y) => (!x && !y ?
  { x: 0, y: 0 } : (typeof x === 'object' ?
    { x: x.x || 0, y: x.y || 0 } : (y === null || y === undefined ?
      { x: x, y: x } : { x: x, y: y })
  )
);

/**
 * Get the components of a vector as an array
 * @param {vec} a The vector to get components from
 * @return {Array<number>} The vector components as an array
 */
vec.components = a => [a.x, a.y];

/**
 * Return a unit vector (1, 0)
 * @return {vec} A unit vector (1, 0)
 */
vec.ux = () => vec(1, 0);

/**
 * Return a unit vector (0, 1)
 * @return {vec} A unit vector (0, 1)
 */
vec.uy = () => vec(0, 1);

/**
 * Add vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {vec} a + b
 */
vec.add = (a, b) => ({ x: a.x + b.x, y: a.y + b.y });

/**
 * Scale a vector
 * @param {vec} a Vector a
 * @param {number} b Scalar b
 * @return {vec} a * b
 */
vec.mul = (a, b) => ({ x: a.x * b, y: a.y * b });

/**
 * Subtract vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {vec} a - b
 */
vec.sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });

/**
 * Get the length of a vector
 * @param {vec} a Vector a
 * @return {number} |a|
 */
vec.len = a => Math.sqrt(a.x * a.x + a.y * a.y);

/**
 * Get the length of a vector using taxicab geometry
 * @param {vec} a Vector a
 * @return {number} |a|
 */
vec.manhattan = a => Math.abs(a.x) + Math.abs(a.y);

/**
 * Normalise a vector
 * @param {vec} a The vector to normalise
 * @return {vec} ^a
 */
vec.nor = a => {
  let len = vec.len(a);
  return len ? { x: a.x / len, y: a.y / len } : vec();
};

/**
 * Get a dot product of vectors
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {number} a ∙ b
 */
vec.dot = (a, b) => a.x * b.x + a.y * b.y;

/**
 * Rotate a vector by r radians
 * @param {vec} a The vector to rotate
 * @param {number} r The angle to rotate by, measured in radians
 * @return {vec} A rotated vector
 */
vec.rot = (a, r) => {
  let s = Math.sin(r),
    c = Math.cos(r);
  return { x: c * a.x - s * a.y, y: s * a.x + c * a.y };
}

/**
 * Check if two vectors are equal
 * @param {vec} a Vector a
 * @param {vec} b Vector b
 * @return {boolean} True if vectors a and b are equal, false otherwise
 */
vec.eq = (a, b) => a.x === b.x && a.y === b.y;

/**
 * Get the angle of a vector
 * @param {vec} a Vector a
 * @return {number} The angle of vector a in radians
 */
vec.rad = a => Math.atan2(a.y, a.x);

/**
 * Copy a vector
 * @param {vec} a The vector to copy
 * @return {vec} A copy of vector a
 */
vec.cpy = a => vec(a);

/**
 * A function to call on each component of a vector
 * @callback vectorMapCallback
 * @param {number} value The component value
 * @param {'x' | 'y'} label The component label (x or y)
 * @return {number} The mapped component
 */

/**
 * Call a function on each component of a vector and build a new vector from the results
 * @param {vec} a Vector a
 * @param {vectorMapCallback} f The function to call on each component of the vector
 * @return {vec} Vector a mapped through f
 */
vec.map = (a, f) => ({ x: f(a.x, 'x'), y: f(a.y, 'y') });

/**
 * Convert a vector into a string
 * @param {vec} a The vector to convert
 * @param {string} [s=', '] The separator string
 * @return {string} A string representation of the vector
 */
vec.str = (a, s = ', ') => `${a.x}${s}${a.y}`;

/**
 * A matrix
 * @typedef {Object} mat
 * @property {number} m The number of rows in the matrix
 * @property {number} n The number of columns in the matrix
 * @property {Array<number>} entries The matrix values
 */

/**
 * Create a new matrix
 * @param {number} [m=4] The number of rows
 * @param {number} [n=4] The number of columns
 * @param {Array<number>} [entries=[]] Matrix values in reading order
 * @return {mat} A new matrix
 */
const mat = (m = 4, n = 4, entries = []) => ({
  m, n,
  entries: entries.concat(Array(m * n).fill(0)).slice(0, m * n)
});

/**
 * Get an identity matrix of size n
 * @param {number} n The size of the matrix
 * @return {mat} An identity matrix
 */
mat.identity = n => mat(n, n, Array(n * n).fill(0).map((v, i) => +(Math.floor(i / n) === i % n)));

/**
 * Get an entry from a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @return {number} The value at position (i, j) in matrix a
 */
mat.get = (a, i, j) => a.entries[(j - 1) + (i - 1) * a.n];

/**
 * Set an entry of a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @param {number} v The value to set in matrix a
 */
mat.set = (a, i, j, v) => { a.entries[(j - 1) + (i - 1) * a.n] = v; };

/**
 * Get a row from a matrix as an array
 * @param {mat} a Matrix a
 * @param {number} m The row offset
 * @return {Array<number>} Row m from matrix a
 */
mat.row = (a, m) => {
  const s = (m - 1) * a.n;
  return a.entries.slice(s, s + a.n);
};

/**
 * Get a column from a matrix as an array
 * @param {mat} a Matrix a
 * @param {number} n The column offset
 * @return {Array<number>} Column n from matrix a
 */
mat.col = (a, n) => times(i => mat.get(a, (i + 1), n), a.m);

/**
 * Add matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat} a + b
 */
mat.add = (a, b) => a.m === b.m && a.n === b.n && mat.map(a, (v, i) => v + b.entries[i]);

/**
 * Subtract matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat} a - b
 */
mat.sub = (a, b) => a.m === b.m && a.n === b.n && mat.map(a, (v, i) => v - b.entries[i]);

/**
 * Multiply matrices
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {mat|boolean} ab or false if the matrices cannot be multiplied
 */
mat.mul = (a, b) => {
  if (a.n !== b.m) { return false; }
  const result = mat(a.m, b.n);
  for (let i = 1; i <= a.m; i++) {
    for (let j = 1; j <= b.n; j++) {
      mat.set(result, i, j, dot(mat.row(a, i), mat.col(b, j)));
    }
  }
  return result;
};

/**
 * Scale a matrix
 * @param {mat} a Matrix a
 * @param {number} b Scalar b
 * @return {mat} a * b
 */
mat.scale = (a, b) => mat.map(a, v => v * b);

/**
 * Transpose a matrix
 * @param {mat} a The matrix to transpose
 * @return {mat} A transposed matrix
 */
mat.trans = a => mat(a.n, a.m, times(i => mat.col(a, (i + 1)), a.n).flat());

/**
 * Get the minor of a matrix
 * @param {mat} a Matrix a
 * @param {number} i The row offset
 * @param {number} j The column offset
 * @return {mat|boolean} The (i, j) minor of matrix a or false if the matrix is not square
 */
mat.minor = (a, i, j) => {
  if (a.m !== a.n) { return false; }
  const entries = [];
  for (let ii = 1; ii <= a.m; ii++) {
    if (ii === i) { continue; }
    for (let jj = 1; jj <= a.n; jj++) {
      if (jj === j) { continue; }
      entries.push(mat.get(a, ii, jj));
    }
  }
  return mat(a.m - 1, a.n - 1, entries);
};

/**
 * Get the determinant of a matrix
 * @param {mat} a Matrix a
 * @return {number|boolean} |a| or false if the matrix is not square
 */
mat.det = a => {
  if (a.m !== a.n) { return false; }
  if (a.m === 1) {
    return a.entries[0];
  }
  if (a.m === 2) {
    return a.entries[0] * a.entries[3] - a.entries[1] * a.entries[2];
  }
  let total = 0, sign = 1;
  for (let j = 1; j <= a.n; j++) {
    total += sign * a.entries[j - 1] * mat.det(mat.minor(a, 1, j));
    sign *= -1;
  }
  return total;
};

/**
 * Normalise a matrix
 * @param {mat} a The matrix to normalise
 * @return {mat|boolean} ^a or false if the matrix is not square
 */
mat.nor = a => {
  if (a.m !== a.n) { return false; }
  const d = mat.det(a);
  return mat.map(a, i => i * d);
};

/**
 * Get the adjugate of a matrix
 * @param {mat} a The matrix from which to get the adjugate
 * @return {mat} The adjugate of a
 */
mat.adj = a => {
  const minors = mat(a.m, a.n);
  for (let i = 1; i <= a.m; i++) {
    for (let j = 1; j <= a.n; j++) {
      mat.set(minors, i, j, mat.det(mat.minor(a, i, j)));
    }
  }
  const cofactors = mat.map(minors, (v, i) => v * (i % 2 ? -1 : 1));
  return mat.trans(cofactors);
};

/**
 * Get the inverse of a matrix
 * @param {mat} a The matrix to invert
 * @return {mat|boolean} a^-1 or false if the matrix has no inverse
 */
mat.inv = a => {
  if (a.m !== a.n) { return false; }
  const d = mat.det(a);
  if (d === 0) { return false; }
  return mat.scale(mat.adj(a), 1 / d);
};

/**
 * Check if two matrices are equal
 * @param {mat} a Matrix a
 * @param {mat} b Matrix b
 * @return {boolean} True if matrices a and b are identical, false otherwise
 */
mat.eq = (a, b) => a.m === b.m && a.n === b.n && mat.str(a) === mat.str(b);

/**
 * Copy a matrix
 * @param {mat} a The matrix to copy
 * @return {mat} A copy of matrix a
 */
mat.cpy = a => mat(a.m, a.n, [...a.entries]);

/**
 * A function to call on each entry of a matrix
 * @callback matrixMapCallback
 * @param {number} value The entry value
 * @param {number} index The entry index
 * @param {Array<number>} entries The array of matrix entries
 * @return {number} The mapped entry
 */

/**
 * Call a function on each entry of a matrix and build a new matrix from the results
 * @param {mat} a Matrix a
 * @param {matrixMapCallback} f The function to call on each entry of the matrix
 * @return {mat} Matrix a mapped through f
 */
mat.map = (a, f) => mat(a.m, a.n, a.entries.map(f));

/**
 * Convert a matrix into a string
 * @param {mat} a The matrix to convert
 * @param {string} [ms=', '] The separator string for columns
 * @param {string} [ns='\n'] The separator string for rows
 * @return {string} A string representation of the matrix
 */
mat.str = (a, ms = ', ', ns = '\n') => chunk(a.entries, a.n).map(r => r.join(ms)).join(ns);

if (true) {
  module.exports = { vec, mat };
}


/***/ }),

/***/ "./index.ts":
/*!******************!*\
  !*** ./index.ts ***!
  \******************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.Sprite = exports.SpriteAnimationRepeatMode = void 0;
const vec_1 = __webpack_require__(/*! @basementuniverse/vec */ "./node_modules/@basementuniverse/vec/vec.js");
var SpriteAnimationRepeatMode;
(function (SpriteAnimationRepeatMode) {
    /**
     * Loop this animation indefinitely
     */
    SpriteAnimationRepeatMode[SpriteAnimationRepeatMode["Repeat"] = 0] = "Repeat";
    /**
     * Play once and then stop on the last frame
     */
    SpriteAnimationRepeatMode[SpriteAnimationRepeatMode["PlayOnceAndStop"] = 1] = "PlayOnceAndStop";
    /**
     * Play once and then reset back to the first frame
     */
    SpriteAnimationRepeatMode[SpriteAnimationRepeatMode["PlayOnceAndReset"] = 2] = "PlayOnceAndReset";
})(SpriteAnimationRepeatMode = exports.SpriteAnimationRepeatMode || (exports.SpriteAnimationRepeatMode = {}));
class Sprite {
    constructor(options) {
        var _a, _b;
        this.position = (0, vec_1.vec)();
        this.size = (0, vec_1.vec)();
        this.origin = (0, vec_1.vec)();
        this.scale = 1;
        this.rotation = 0;
        this.currentAnimationOptions = null;
        this.currentAnimationState = null;
        this.currentImage = null;
        this.currentAttachmentPoints = null;
        const actualOptions = Object.assign({}, Sprite.DEFAULT_OPTIONS, options !== null && options !== void 0 ? options : {});
        if (!actualOptions.debug || actualOptions.debug === true) {
            actualOptions.debug = {
                showSpriteTransforms: !!actualOptions.debug,
                showSpriteBoundingBox: !!actualOptions.debug,
                showAttachmentPoints: !!actualOptions.debug,
            };
        }
        this.options = actualOptions;
        if (this.options.position) {
            this.position = vec_1.vec.cpy(this.options.position);
        }
        if (this.options.size) {
            this.size = vec_1.vec.cpy(this.options.size);
        }
        else {
            // Default to the size of the base image if one exists
            if (this.options.image) {
                this.size = (0, vec_1.vec)(this.options.image.width, this.options.image.height);
            }
            else {
                // Fall back to the size of the image in the first frame of the first
                // available direction of the default animation if one exists
                const defaultAnimationDirections = Object.values(this.options.animations[this.options.defaultAnimation])[0];
                if (defaultAnimationDirections &&
                    ((_b = (_a = defaultAnimationDirections.images) === null || _a === void 0 ? void 0 : _a.length) !== null && _b !== void 0 ? _b : 0) > 0) {
                    this.size = (0, vec_1.vec)(defaultAnimationDirections.images[0].width, defaultAnimationDirections.images[0].height);
                }
            }
            // Otherwise leave the size as (0, 0)
        }
        if (this.options.origin) {
            this.origin = vec_1.vec.cpy(this.options.origin);
        }
        else {
            // Default to the center of the sprite based on size
            this.origin = vec_1.vec.mul(this.size, 0.5);
        }
        if (this.options.scale) {
            this.scale = this.options.scale;
        }
        if (this.options.rotation) {
            this.rotation = this.options.rotation;
        }
        // Check and initialise direction
        this._direction = this.options.defaultDirection;
        if (this.options.directions.length === 0 ||
            !this.options.directions.includes(this._direction)) {
            throw new Error(`Invalid direction "${this._direction}"`);
        }
        // Check and initialise animation
        this._animation = this.options.defaultAnimation;
        const animations = Object.keys(this.options.animations);
        if (animations.length === 0 ||
            !animations.includes(this._animation)) {
            throw new Error(`Invalid animation "${this._animation}"`);
        }
        // Make sure attachment point keyframes are defined in ascending
        // frame order in all animations
        for (const animation of Object.keys(this.options.animations)) {
            for (const direction of Object.keys(this.options.animations[animation])) {
                if (this.options.animations[animation][direction].attachmentPointKeyframes) {
                    for (const attachmentPoint of Object.keys(this
                        .options
                        .animations[animation][direction]
                        .attachmentPointKeyframes)) {
                        this
                            .options
                            .animations[animation][direction]
                            .attachmentPointKeyframes[attachmentPoint]
                            .sort((a, b) => a.frame - b.frame);
                    }
                }
            }
        }
    }
    get direction() {
        return this._direction;
    }
    set direction(value) {
        if (this.options.directions.includes(value)) {
            this._direction = value;
        }
    }
    get animation() {
        return this._animation;
    }
    set animation(value) {
        if (Object.keys(this.options.animations).includes(value)) {
            this._animation = value;
        }
    }
    playAnimation() {
        if (this.currentAnimationState) {
            this.currentAnimationState.playing = true;
        }
    }
    pauseAnimation() {
        if (this.currentAnimationState) {
            this.currentAnimationState.playing = false;
        }
    }
    resetAnimation() {
        if (this.currentAnimationState) {
            this.currentAnimationState.currentFrame = 0;
            this.currentAnimationState.currentFrameTime = 0;
        }
    }
    getAttachmentPoint(name) {
        var _a, _b;
        return (_b = (_a = this.currentAttachmentPoints) === null || _a === void 0 ? void 0 : _a[name]) !== null && _b !== void 0 ? _b : null;
    }
    update(dt) {
        this.currentAnimationOptions = this.updateAnimationOptions();
        this.currentAnimationState = this.updateAnimationState(dt);
        this.currentImage = this.updateImage();
        this.currentAttachmentPoints = this.updateAttachmentPoints();
    }
    updateAnimationOptions() {
        if (!(this._animation in this.options.animations)) {
            throw new Error(`Invalid animation "${this._animation}"`);
        }
        const directions = Object.keys(this.options.animations[this._animation]);
        if (directions.length === 0) {
            throw new Error(`No directions available for animation "${this._animation}"`);
        }
        if (this._direction in this.options.animations[this._animation]) {
            return this.options.animations[this._animation][this._direction];
        }
        if ('*' in this.options.animations[this._animation]) {
            return this.options.animations[this._animation]['*'];
        }
        return this.options.animations[this._animation][directions[0]];
    }
    updateAnimationState(dt) {
        if (!this.currentAnimationOptions ||
            !this.currentAnimationState) {
            return {
                playing: true,
                currentFrame: 0,
                currentFrameTime: 0,
            };
        }
        if (this.currentAnimationState.playing) {
            const frameTime = 1 / this.currentAnimationOptions.frameRate;
            this.currentAnimationState.currentFrameTime += dt;
            if (this.currentAnimationState.currentFrameTime > frameTime) {
                const frameCount = this.currentAnimationOptions.frameCount;
                this.currentAnimationState.currentFrame++;
                this.currentAnimationState.currentFrameTime = 0;
                if (this.currentAnimationState.currentFrame > frameCount) {
                    switch (this.currentAnimationOptions.mode) {
                        case SpriteAnimationRepeatMode.PlayOnceAndReset:
                            this.currentAnimationState.playing = false;
                            this.currentAnimationState.currentFrame = 0;
                            break;
                        case SpriteAnimationRepeatMode.PlayOnceAndStop:
                            this.currentAnimationState.playing = false;
                            this.currentAnimationState.currentFrame = frameCount - 1;
                            break;
                        case SpriteAnimationRepeatMode.Repeat:
                            this.currentAnimationState.currentFrame = 0;
                            break;
                    }
                }
            }
        }
        return this.currentAnimationState;
    }
    updateImage() {
        var _a, _b, _c;
        if (!this.currentAnimationOptions ||
            !this.currentAnimationState) {
            return null;
        }
        if (!this.currentAnimationOptions.images ||
            this.currentAnimationOptions.images.length === 0) {
            return (_a = this.options.image) !== null && _a !== void 0 ? _a : null;
        }
        return (_c = (_b = this.currentAnimationOptions.images[this.currentAnimationState.currentFrame]) !== null && _b !== void 0 ? _b : this.options.image) !== null && _c !== void 0 ? _c : null;
    }
    updateAttachmentPoints() {
        if (!this.options.attachmentPoints ||
            this.options.attachmentPoints.length === 0) {
            return null;
        }
        if (!this.currentAttachmentPoints) {
            this.currentAttachmentPoints = Object.fromEntries(this.options.attachmentPoints.map(attachmentPoint => [
                attachmentPoint.name,
                attachmentPoint.offset
            ]));
        }
        if (this.currentAnimationOptions &&
            this.currentAnimationOptions.attachmentPointKeyframes &&
            this.currentAnimationState) {
            for (const name of Object.keys(this.currentAttachmentPoints)) {
                if (name in this.currentAnimationOptions.attachmentPointKeyframes &&
                    this.currentAnimationOptions.attachmentPointKeyframes[name].length > 0) {
                    const previousKeyframe = this.findPreviousKeyframe(this.currentAnimationOptions.attachmentPointKeyframes[name], this.currentAnimationState.currentFrame);
                    this.currentAttachmentPoints[name] = previousKeyframe.offset;
                }
            }
        }
        return this.currentAttachmentPoints;
    }
    findPreviousKeyframe(keyframes, currentFrame) {
        const found = [...keyframes].reverse().find(keyframe => keyframe.frame <= currentFrame);
        if (!found) {
            return keyframes[keyframes.length - 1];
        }
        return found;
    }
    draw(context) {
        context.save();
        context.translate(this.position.x, this.position.y);
        context.scale(this.scale, this.scale);
        context.rotate(this.rotation);
        if (this.options.preRender) {
            this.options.preRender(context, this);
        }
        if (this.currentImage) {
            context.drawImage(this.currentImage, -this.origin.x, -this.origin.y, this.currentImage.width, this.currentImage.height);
        }
        if (this.options.postRender) {
            this.options.postRender(context, this);
        }
        if (this.options.debug.showSpriteBoundingBox) {
            context.strokeStyle = Sprite.DEBUG_BOUNDING_BOX_COLOUR;
            context.lineWidth = Sprite.DEBUG_BOUNDING_BOX_LINE_WIDTH;
            context.strokeRect(-this.origin.x, -this.origin.y, this.size.x, this.size.y);
        }
        if (this.options.debug.showSpriteTransforms) {
            this.drawTransformsMarker(context, (0, vec_1.vec)(), Sprite.DEBUG_TRANSFORMS_COLOUR_X, Sprite.DEBUG_TRANSFORMS_COLOUR_Y, Sprite.DEBUG_TRANSFORMS_LINE_WIDTH, Sprite.DEBUG_TRANSFORMS_SIZE);
        }
        if (this.options.debug.showAttachmentPoints &&
            this.currentAttachmentPoints) {
            for (const attachmentPoint of Object.values(this.currentAttachmentPoints)) {
                this.drawCross(context, attachmentPoint, Sprite.DEBUG_ATTACHMENT_POINT_COLOUR, Sprite.DEBUG_ATTACHMENT_POINT_LINE_WIDTH, Sprite.DEBUG_ATTACHMENT_POINT_SIZE);
            }
        }
        context.restore();
    }
    drawTransformsMarker(context, position, xColour, yColour, lineWidth, size) {
        context.save();
        context.lineWidth = lineWidth;
        context.strokeStyle = xColour;
        context.beginPath();
        context.moveTo(position.x, position.y);
        context.lineTo(position.x + size, position.y);
        context.stroke();
        context.strokeStyle = yColour;
        context.beginPath();
        context.moveTo(position.x, position.y);
        context.lineTo(position.x, position.y + size);
        context.stroke();
        context.restore();
    }
    drawCross(context, position, colour, lineWidth, size) {
        context.save();
        context.lineWidth = lineWidth;
        const halfSize = Math.ceil(size / 2);
        context.strokeStyle = colour;
        context.beginPath();
        context.moveTo(position.x - halfSize, position.y - halfSize);
        context.lineTo(position.x + halfSize, position.y + halfSize);
        context.moveTo(position.x - halfSize, position.y + halfSize);
        context.lineTo(position.x + halfSize, position.y - halfSize);
        context.stroke();
        context.restore();
    }
}
exports.Sprite = Sprite;
Sprite.DEFAULT_OPTIONS = {
    directions: ['default'],
    defaultDirection: 'default',
    animations: {
        default: {
            '*': {
                name: 'default',
                frameCount: 1,
                frameRate: 1,
                mode: SpriteAnimationRepeatMode.PlayOnceAndStop,
            },
        },
    },
    defaultAnimation: 'default',
};
Sprite.DEBUG_BOUNDING_BOX_COLOUR = 'green';
Sprite.DEBUG_BOUNDING_BOX_LINE_WIDTH = 2;
Sprite.DEBUG_TRANSFORMS_COLOUR_X = 'red';
Sprite.DEBUG_TRANSFORMS_COLOUR_Y = 'orange';
Sprite.DEBUG_TRANSFORMS_LINE_WIDTH = 1;
Sprite.DEBUG_TRANSFORMS_SIZE = 10;
Sprite.DEBUG_ATTACHMENT_POINT_COLOUR = 'blue';
Sprite.DEBUG_ATTACHMENT_POINT_LINE_WIDTH = 2;
Sprite.DEBUG_ATTACHMENT_POINT_SIZE = 5;
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiaW5kZXguanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlcyI6WyIuLi9pbmRleC50cyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiOzs7QUFBQSwrQ0FBNEM7QUFnSTVDLElBQVkseUJBZVg7QUFmRCxXQUFZLHlCQUF5QjtJQUNuQzs7T0FFRztJQUNILDZFQUFVLENBQUE7SUFFVjs7T0FFRztJQUNILCtGQUFlLENBQUE7SUFFZjs7T0FFRztJQUNILGlHQUFnQixDQUFBO0FBQ2xCLENBQUMsRUFmVyx5QkFBeUIsR0FBekIsaUNBQXlCLEtBQXpCLGlDQUF5QixRQWVwQztBQTBGRCxNQUFhLE1BQU07SUFvRGpCLFlBQW1CLE9BQWdDOztRQWY1QyxhQUFRLEdBQVEsSUFBQSxTQUFHLEdBQUUsQ0FBQztRQUN0QixTQUFJLEdBQVEsSUFBQSxTQUFHLEdBQUUsQ0FBQztRQUVsQixXQUFNLEdBQVEsSUFBQSxTQUFHLEdBQUUsQ0FBQztRQUNwQixVQUFLLEdBQVcsQ0FBQyxDQUFDO1FBQ2xCLGFBQVEsR0FBVyxDQUFDLENBQUM7UUFLcEIsNEJBQXVCLEdBQWtDLElBQUksQ0FBQztRQUM5RCwwQkFBcUIsR0FBZ0MsSUFBSSxDQUFDO1FBQzFELGlCQUFZLEdBQWdELElBQUksQ0FBQztRQUNqRSw0QkFBdUIsR0FBb0MsSUFBSSxDQUFDO1FBR3RFLE1BQU0sYUFBYSxHQUFHLE1BQU0sQ0FBQyxNQUFNLENBQ2pDLEVBQUUsRUFDRixNQUFNLENBQUMsZUFBZSxFQUN0QixPQUFPLGFBQVAsT0FBTyxjQUFQLE9BQU8sR0FBSSxFQUFFLENBQ2QsQ0FBQztRQUVGLElBQUksQ0FBQyxhQUFhLENBQUMsS0FBSyxJQUFJLGFBQWEsQ0FBQyxLQUFLLEtBQUssSUFBSSxFQUFFO1lBQ3hELGFBQWEsQ0FBQyxLQUFLLEdBQUc7Z0JBQ3BCLG9CQUFvQixFQUFFLENBQUMsQ0FBQyxhQUFhLENBQUMsS0FBSztnQkFDM0MscUJBQXFCLEVBQUUsQ0FBQyxDQUFDLGFBQWEsQ0FBQyxLQUFLO2dCQUM1QyxvQkFBb0IsRUFBRSxDQUFDLENBQUMsYUFBYSxDQUFDLEtBQUs7YUFDNUMsQ0FBQztTQUNIO1FBRUQsSUFBSSxDQUFDLE9BQU8sR0FBRyxhQUFvQyxDQUFDO1FBRXBELElBQUksSUFBSSxDQUFDLE9BQU8sQ0FBQyxRQUFRLEVBQUU7WUFDekIsSUFBSSxDQUFDLFFBQVEsR0FBRyxTQUFHLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsUUFBUSxDQUFDLENBQUM7U0FDaEQ7UUFFRCxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsSUFBSSxFQUFFO1lBQ3JCLElBQUksQ0FBQyxJQUFJLEdBQUcsU0FBRyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLElBQUksQ0FBQyxDQUFDO1NBQ3hDO2FBQU07WUFDTCxzREFBc0Q7WUFDdEQsSUFBSSxJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssRUFBRTtnQkFDdEIsSUFBSSxDQUFDLElBQUksR0FBRyxJQUFBLFNBQUcsRUFDYixJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssQ0FBQyxLQUFLLEVBQ3hCLElBQUksQ0FBQyxPQUFPLENBQUMsS0FBSyxDQUFDLE1BQU0sQ0FDMUIsQ0FBQzthQUNIO2lCQUFNO2dCQUNMLHFFQUFxRTtnQkFDckUsNkRBQTZEO2dCQUM3RCxNQUFNLDBCQUEwQixHQUFHLE1BQU0sQ0FBQyxNQUFNLENBQzlDLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsZ0JBQWdCLENBQUMsQ0FDdkQsQ0FBQyxDQUFDLENBQUMsQ0FBQztnQkFDTCxJQUNFLDBCQUEwQjtvQkFDMUIsQ0FBQyxNQUFBLE1BQUEsMEJBQTBCLENBQUMsTUFBTSwwQ0FBRSxNQUFNLG1DQUFJLENBQUMsQ0FBQyxHQUFHLENBQUMsRUFDcEQ7b0JBQ0EsSUFBSSxDQUFDLElBQUksR0FBRyxJQUFBLFNBQUcsRUFDYiwwQkFBMEIsQ0FBQyxNQUFPLENBQUMsQ0FBQyxDQUFDLENBQUMsS0FBSyxFQUMzQywwQkFBMEIsQ0FBQyxNQUFPLENBQUMsQ0FBQyxDQUFDLENBQUMsTUFBTSxDQUM3QyxDQUFDO2lCQUNIO2FBQ0Y7WUFFRCxxQ0FBcUM7U0FDdEM7UUFFRCxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsTUFBTSxFQUFFO1lBQ3ZCLElBQUksQ0FBQyxNQUFNLEdBQUcsU0FBRyxDQUFDLEdBQUcsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLE1BQU0sQ0FBQyxDQUFDO1NBQzVDO2FBQU07WUFDTCxvREFBb0Q7WUFDcEQsSUFBSSxDQUFDLE1BQU0sR0FBRyxTQUFHLENBQUMsR0FBRyxDQUFDLElBQUksQ0FBQyxJQUFJLEVBQUUsR0FBRyxDQUFDLENBQUM7U0FDdkM7UUFFRCxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsS0FBSyxFQUFFO1lBQ3RCLElBQUksQ0FBQyxLQUFLLEdBQUcsSUFBSSxDQUFDLE9BQU8sQ0FBQyxLQUFLLENBQUM7U0FDakM7UUFFRCxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsUUFBUSxFQUFFO1lBQ3pCLElBQUksQ0FBQyxRQUFRLEdBQUcsSUFBSSxDQUFDLE9BQU8sQ0FBQyxRQUFRLENBQUM7U0FDdkM7UUFFRCxpQ0FBaUM7UUFDakMsSUFBSSxDQUFDLFVBQVUsR0FBRyxJQUFJLENBQUMsT0FBTyxDQUFDLGdCQUFnQixDQUFDO1FBQ2hELElBQ0UsSUFBSSxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsTUFBTSxLQUFLLENBQUM7WUFDcEMsQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxRQUFRLENBQUMsSUFBSSxDQUFDLFVBQVUsQ0FBQyxFQUNsRDtZQUNBLE1BQU0sSUFBSSxLQUFLLENBQUMsc0JBQXNCLElBQUksQ0FBQyxVQUFVLEdBQUcsQ0FBQyxDQUFDO1NBQzNEO1FBRUQsaUNBQWlDO1FBQ2pDLElBQUksQ0FBQyxVQUFVLEdBQUcsSUFBSSxDQUFDLE9BQU8sQ0FBQyxnQkFBZ0IsQ0FBQztRQUNoRCxNQUFNLFVBQVUsR0FBRyxNQUFNLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLENBQUM7UUFDeEQsSUFDRSxVQUFVLENBQUMsTUFBTSxLQUFLLENBQUM7WUFDdkIsQ0FBQyxVQUFVLENBQUMsUUFBUSxDQUFDLElBQUksQ0FBQyxVQUFVLENBQUMsRUFDckM7WUFDQSxNQUFNLElBQUksS0FBSyxDQUFDLHNCQUFzQixJQUFJLENBQUMsVUFBVSxHQUFHLENBQUMsQ0FBQztTQUMzRDtRQUVELGdFQUFnRTtRQUNoRSxnQ0FBZ0M7UUFDaEMsS0FBSyxNQUFNLFNBQVMsSUFBSSxNQUFNLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLEVBQUU7WUFDNUQsS0FBSyxNQUFNLFNBQVMsSUFBSSxNQUFNLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLFNBQVMsQ0FBQyxDQUFDLEVBQUU7Z0JBQ3ZFLElBQ0UsSUFBSSxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsU0FBUyxDQUFDLENBQUMsU0FBUyxDQUFDLENBQUMsd0JBQXdCLEVBQ3RFO29CQUNBLEtBQUssTUFBTSxlQUFlLElBQUksTUFBTSxDQUFDLElBQUksQ0FDdkMsSUFBSTt5QkFDRCxPQUFPO3lCQUNQLFVBQVUsQ0FBQyxTQUFTLENBQUMsQ0FBQyxTQUFTLENBQUM7eUJBQ2hDLHdCQUF5QixDQUM3QixFQUFFO3dCQUNELElBQUk7NkJBQ0QsT0FBTzs2QkFDUCxVQUFVLENBQUMsU0FBUyxDQUFDLENBQUMsU0FBUyxDQUFDOzZCQUNoQyx3QkFBeUIsQ0FBQyxlQUFlLENBQUM7NkJBQzFDLElBQUksQ0FDSCxDQUFDLENBQUMsRUFBRSxDQUFDLEVBQUUsRUFBRSxDQUFDLENBQUMsQ0FBQyxLQUFLLEdBQUcsQ0FBQyxDQUFDLEtBQUssQ0FDNUIsQ0FBQztxQkFDTDtpQkFDRjthQUNGO1NBQ0Y7SUFDSCxDQUFDO0lBRUQsSUFBVyxTQUFTO1FBQ2xCLE9BQU8sSUFBSSxDQUFDLFVBQVUsQ0FBQztJQUN6QixDQUFDO0lBRUQsSUFBVyxTQUFTLENBQUMsS0FBYTtRQUNoQyxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLFFBQVEsQ0FBQyxLQUFLLENBQUMsRUFBRTtZQUMzQyxJQUFJLENBQUMsVUFBVSxHQUFHLEtBQUssQ0FBQztTQUN6QjtJQUNILENBQUM7SUFFRCxJQUFXLFNBQVM7UUFDbEIsT0FBTyxJQUFJLENBQUMsVUFBVSxDQUFDO0lBQ3pCLENBQUM7SUFFRCxJQUFXLFNBQVMsQ0FBQyxLQUFhO1FBQ2hDLElBQUksTUFBTSxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxDQUFDLFFBQVEsQ0FBQyxLQUFLLENBQUMsRUFBRTtZQUN4RCxJQUFJLENBQUMsVUFBVSxHQUFHLEtBQUssQ0FBQztTQUN6QjtJQUNILENBQUM7SUFFTSxhQUFhO1FBQ2xCLElBQUksSUFBSSxDQUFDLHFCQUFxQixFQUFFO1lBQzlCLElBQUksQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLEdBQUcsSUFBSSxDQUFDO1NBQzNDO0lBQ0gsQ0FBQztJQUVNLGNBQWM7UUFDbkIsSUFBSSxJQUFJLENBQUMscUJBQXFCLEVBQUU7WUFDOUIsSUFBSSxDQUFDLHFCQUFxQixDQUFDLE9BQU8sR0FBRyxLQUFLLENBQUM7U0FDNUM7SUFDSCxDQUFDO0lBRU0sY0FBYztRQUNuQixJQUFJLElBQUksQ0FBQyxxQkFBcUIsRUFBRTtZQUM5QixJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxHQUFHLENBQUMsQ0FBQztZQUM1QyxJQUFJLENBQUMscUJBQXFCLENBQUMsZ0JBQWdCLEdBQUcsQ0FBQyxDQUFDO1NBQ2pEO0lBQ0gsQ0FBQztJQUVNLGtCQUFrQixDQUFDLElBQVk7O1FBQ3BDLE9BQU8sTUFBQSxNQUFBLElBQUksQ0FBQyx1QkFBdUIsMENBQUcsSUFBSSxDQUFDLG1DQUFJLElBQUksQ0FBQztJQUN0RCxDQUFDO0lBRU0sTUFBTSxDQUFDLEVBQVU7UUFDdEIsSUFBSSxDQUFDLHVCQUF1QixHQUFHLElBQUksQ0FBQyxzQkFBc0IsRUFBRSxDQUFDO1FBQzdELElBQUksQ0FBQyxxQkFBcUIsR0FBRyxJQUFJLENBQUMsb0JBQW9CLENBQUMsRUFBRSxDQUFDLENBQUM7UUFDM0QsSUFBSSxDQUFDLFlBQVksR0FBRyxJQUFJLENBQUMsV0FBVyxFQUFFLENBQUM7UUFDdkMsSUFBSSxDQUFDLHVCQUF1QixHQUFHLElBQUksQ0FBQyxzQkFBc0IsRUFBRSxDQUFDO0lBQy9ELENBQUM7SUFFTyxzQkFBc0I7UUFDNUIsSUFBSSxDQUFDLENBQUMsSUFBSSxDQUFDLFVBQVUsSUFBSSxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxFQUFFO1lBQ2pELE1BQU0sSUFBSSxLQUFLLENBQUMsc0JBQXNCLElBQUksQ0FBQyxVQUFVLEdBQUcsQ0FBQyxDQUFDO1NBQzNEO1FBRUQsTUFBTSxVQUFVLEdBQUcsTUFBTSxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxJQUFJLENBQUMsVUFBVSxDQUFDLENBQUMsQ0FBQztRQUN6RSxJQUFJLFVBQVUsQ0FBQyxNQUFNLEtBQUssQ0FBQyxFQUFFO1lBQzNCLE1BQU0sSUFBSSxLQUFLLENBQ2IsMENBQTBDLElBQUksQ0FBQyxVQUFVLEdBQUcsQ0FDN0QsQ0FBQztTQUNIO1FBRUQsSUFBSSxJQUFJLENBQUMsVUFBVSxJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLElBQUksQ0FBQyxVQUFVLENBQUMsRUFBRTtZQUMvRCxPQUFPLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLElBQUksQ0FBQyxVQUFVLENBQUMsQ0FBQyxJQUFJLENBQUMsVUFBVSxDQUFDLENBQUM7U0FDbEU7UUFFRCxJQUFJLEdBQUcsSUFBSSxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxJQUFJLENBQUMsVUFBVSxDQUFDLEVBQUU7WUFDbkQsT0FBTyxJQUFJLENBQUMsT0FBTyxDQUFDLFVBQVUsQ0FBQyxJQUFJLENBQUMsVUFBVSxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUM7U0FDdEQ7UUFFRCxPQUFPLElBQUksQ0FBQyxPQUFPLENBQUMsVUFBVSxDQUFDLElBQUksQ0FBQyxVQUFVLENBQUMsQ0FBQyxVQUFVLENBQUMsQ0FBQyxDQUFDLENBQUMsQ0FBQztJQUNqRSxDQUFDO0lBRU8sb0JBQW9CLENBQUMsRUFBVTtRQUNyQyxJQUNFLENBQUMsSUFBSSxDQUFDLHVCQUF1QjtZQUM3QixDQUFDLElBQUksQ0FBQyxxQkFBcUIsRUFDM0I7WUFDQSxPQUFPO2dCQUNMLE9BQU8sRUFBRSxJQUFJO2dCQUNiLFlBQVksRUFBRSxDQUFDO2dCQUNmLGdCQUFnQixFQUFFLENBQUM7YUFDcEIsQ0FBQztTQUNIO1FBRUQsSUFBSSxJQUFJLENBQUMscUJBQXFCLENBQUMsT0FBTyxFQUFFO1lBQ3RDLE1BQU0sU0FBUyxHQUFHLENBQUMsR0FBRyxJQUFJLENBQUMsdUJBQXVCLENBQUMsU0FBUyxDQUFDO1lBQzdELElBQUksQ0FBQyxxQkFBcUIsQ0FBQyxnQkFBZ0IsSUFBSSxFQUFFLENBQUM7WUFFbEQsSUFBSSxJQUFJLENBQUMscUJBQXFCLENBQUMsZ0JBQWdCLEdBQUcsU0FBUyxFQUFFO2dCQUMzRCxNQUFNLFVBQVUsR0FBRyxJQUFJLENBQUMsdUJBQXVCLENBQUMsVUFBVSxDQUFDO2dCQUMzRCxJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxFQUFFLENBQUM7Z0JBQzFDLElBQUksQ0FBQyxxQkFBcUIsQ0FBQyxnQkFBZ0IsR0FBRyxDQUFDLENBQUM7Z0JBRWhELElBQUksSUFBSSxDQUFDLHFCQUFxQixDQUFDLFlBQVksR0FBRyxVQUFVLEVBQUU7b0JBQ3hELFFBQVEsSUFBSSxDQUFDLHVCQUF1QixDQUFDLElBQUksRUFBRTt3QkFDekMsS0FBSyx5QkFBeUIsQ0FBQyxnQkFBZ0I7NEJBQzdDLElBQUksQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLEdBQUcsS0FBSyxDQUFDOzRCQUMzQyxJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxHQUFHLENBQUMsQ0FBQzs0QkFDNUMsTUFBTTt3QkFFUixLQUFLLHlCQUF5QixDQUFDLGVBQWU7NEJBQzVDLElBQUksQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLEdBQUcsS0FBSyxDQUFDOzRCQUMzQyxJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxHQUFHLFVBQVUsR0FBRyxDQUFDLENBQUM7NEJBQ3pELE1BQU07d0JBRVIsS0FBSyx5QkFBeUIsQ0FBQyxNQUFNOzRCQUNuQyxJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxHQUFHLENBQUMsQ0FBQzs0QkFDNUMsTUFBTTtxQkFDVDtpQkFDRjthQUNGO1NBQ0Y7UUFFRCxPQUFPLElBQUksQ0FBQyxxQkFBcUIsQ0FBQztJQUNwQyxDQUFDO0lBRU8sV0FBVzs7UUFDakIsSUFDRSxDQUFDLElBQUksQ0FBQyx1QkFBdUI7WUFDN0IsQ0FBQyxJQUFJLENBQUMscUJBQXFCLEVBQzNCO1lBQ0EsT0FBTyxJQUFJLENBQUM7U0FDYjtRQUVELElBQ0UsQ0FBQyxJQUFJLENBQUMsdUJBQXVCLENBQUMsTUFBTTtZQUNwQyxJQUFJLENBQUMsdUJBQXVCLENBQUMsTUFBTSxDQUFDLE1BQU0sS0FBSyxDQUFDLEVBQ2hEO1lBQ0EsT0FBTyxNQUFBLElBQUksQ0FBQyxPQUFPLENBQUMsS0FBSyxtQ0FBSSxJQUFJLENBQUM7U0FDbkM7UUFFRCxPQUFPLE1BQUEsTUFBQSxJQUFJLENBQUMsdUJBQXVCLENBQUMsTUFBTSxDQUN4QyxJQUFJLENBQUMscUJBQXFCLENBQUMsWUFBWSxDQUN4QyxtQ0FBSSxJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssbUNBQUksSUFBSSxDQUFDO0lBQ2xDLENBQUM7SUFFTyxzQkFBc0I7UUFDNUIsSUFDRSxDQUFDLElBQUksQ0FBQyxPQUFPLENBQUMsZ0JBQWdCO1lBQzlCLElBQUksQ0FBQyxPQUFPLENBQUMsZ0JBQWdCLENBQUMsTUFBTSxLQUFLLENBQUMsRUFDMUM7WUFDQSxPQUFPLElBQUksQ0FBQztTQUNiO1FBRUQsSUFBSSxDQUFDLElBQUksQ0FBQyx1QkFBdUIsRUFBRTtZQUNqQyxJQUFJLENBQUMsdUJBQXVCLEdBQUcsTUFBTSxDQUFDLFdBQVcsQ0FDL0MsSUFBSSxDQUFDLE9BQU8sQ0FBQyxnQkFBZ0IsQ0FBQyxHQUFHLENBQUMsZUFBZSxDQUFDLEVBQUUsQ0FBQztnQkFDbkQsZUFBZSxDQUFDLElBQUk7Z0JBQ3BCLGVBQWUsQ0FBQyxNQUFNO2FBQ3ZCLENBQUMsQ0FDSCxDQUFDO1NBQ0g7UUFFRCxJQUNFLElBQUksQ0FBQyx1QkFBdUI7WUFDNUIsSUFBSSxDQUFDLHVCQUF1QixDQUFDLHdCQUF3QjtZQUNyRCxJQUFJLENBQUMscUJBQXFCLEVBQzFCO1lBQ0EsS0FBSyxNQUFNLElBQUksSUFBSSxNQUFNLENBQUMsSUFBSSxDQUFDLElBQUksQ0FBQyx1QkFBdUIsQ0FBQyxFQUFFO2dCQUM1RCxJQUNFLElBQUksSUFBSSxJQUFJLENBQUMsdUJBQXVCLENBQUMsd0JBQXdCO29CQUM3RCxJQUFJLENBQUMsdUJBQXVCLENBQUMsd0JBQXdCLENBQUMsSUFBSSxDQUFDLENBQUMsTUFBTSxHQUFHLENBQUMsRUFDdEU7b0JBQ0EsTUFBTSxnQkFBZ0IsR0FBRyxJQUFJLENBQUMsb0JBQW9CLENBQ2hELElBQUksQ0FBQyx1QkFBdUIsQ0FBQyx3QkFBd0IsQ0FBQyxJQUFJLENBQUMsRUFDM0QsSUFBSSxDQUFDLHFCQUFxQixDQUFDLFlBQVksQ0FDeEMsQ0FBQztvQkFDRixJQUFJLENBQUMsdUJBQXVCLENBQUMsSUFBSSxDQUFDLEdBQUcsZ0JBQWdCLENBQUMsTUFBTSxDQUFDO2lCQUM5RDthQUNGO1NBQ0Y7UUFFRCxPQUFPLElBQUksQ0FBQyx1QkFBdUIsQ0FBQztJQUN0QyxDQUFDO0lBRU8sb0JBQW9CLENBQzFCLFNBQTBDLEVBQzFDLFlBQW9CO1FBRXBCLE1BQU0sS0FBSyxHQUFHLENBQUMsR0FBRyxTQUFTLENBQUMsQ0FBQyxPQUFPLEVBQUUsQ0FBQyxJQUFJLENBQ3pDLFFBQVEsQ0FBQyxFQUFFLENBQUMsUUFBUSxDQUFDLEtBQUssSUFBSSxZQUFZLENBQzNDLENBQUM7UUFFRixJQUFJLENBQUMsS0FBSyxFQUFFO1lBQ1YsT0FBTyxTQUFTLENBQUMsU0FBUyxDQUFDLE1BQU0sR0FBRyxDQUFDLENBQUMsQ0FBQztTQUN4QztRQUVELE9BQU8sS0FBSyxDQUFDO0lBQ2YsQ0FBQztJQUVNLElBQUksQ0FBQyxPQUFpQztRQUMzQyxPQUFPLENBQUMsSUFBSSxFQUFFLENBQUM7UUFDZixPQUFPLENBQUMsU0FBUyxDQUNmLElBQUksQ0FBQyxRQUFRLENBQUMsQ0FBQyxFQUNmLElBQUksQ0FBQyxRQUFRLENBQUMsQ0FBQyxDQUNoQixDQUFDO1FBQ0YsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsS0FBSyxFQUFFLElBQUksQ0FBQyxLQUFLLENBQUMsQ0FBQztRQUN0QyxPQUFPLENBQUMsTUFBTSxDQUFDLElBQUksQ0FBQyxRQUFRLENBQUMsQ0FBQztRQUU5QixJQUFJLElBQUksQ0FBQyxPQUFPLENBQUMsU0FBUyxFQUFFO1lBQzFCLElBQUksQ0FBQyxPQUFPLENBQUMsU0FBUyxDQUFDLE9BQU8sRUFBRSxJQUFJLENBQUMsQ0FBQztTQUN2QztRQUVELElBQUksSUFBSSxDQUFDLFlBQVksRUFBRTtZQUNyQixPQUFPLENBQUMsU0FBUyxDQUNmLElBQUksQ0FBQyxZQUFZLEVBQ2pCLENBQUMsSUFBSSxDQUFDLE1BQU0sQ0FBQyxDQUFDLEVBQ2QsQ0FBQyxJQUFJLENBQUMsTUFBTSxDQUFDLENBQUMsRUFDZCxJQUFJLENBQUMsWUFBWSxDQUFDLEtBQUssRUFDdkIsSUFBSSxDQUFDLFlBQVksQ0FBQyxNQUFNLENBQ3pCLENBQUM7U0FDSDtRQUVELElBQUksSUFBSSxDQUFDLE9BQU8sQ0FBQyxVQUFVLEVBQUU7WUFDM0IsSUFBSSxDQUFDLE9BQU8sQ0FBQyxVQUFVLENBQUMsT0FBTyxFQUFFLElBQUksQ0FBQyxDQUFDO1NBQ3hDO1FBRUQsSUFBSSxJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssQ0FBQyxxQkFBcUIsRUFBRTtZQUM1QyxPQUFPLENBQUMsV0FBVyxHQUFHLE1BQU0sQ0FBQyx5QkFBeUIsQ0FBQztZQUN2RCxPQUFPLENBQUMsU0FBUyxHQUFHLE1BQU0sQ0FBQyw2QkFBNkIsQ0FBQztZQUN6RCxPQUFPLENBQUMsVUFBVSxDQUNoQixDQUFDLElBQUksQ0FBQyxNQUFNLENBQUMsQ0FBQyxFQUNkLENBQUMsSUFBSSxDQUFDLE1BQU0sQ0FBQyxDQUFDLEVBQ2QsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDLEVBQ1gsSUFBSSxDQUFDLElBQUksQ0FBQyxDQUFDLENBQ1osQ0FBQztTQUNIO1FBRUQsSUFBSSxJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssQ0FBQyxvQkFBb0IsRUFBRTtZQUMzQyxJQUFJLENBQUMsb0JBQW9CLENBQ3ZCLE9BQU8sRUFDUCxJQUFBLFNBQUcsR0FBRSxFQUNMLE1BQU0sQ0FBQyx5QkFBeUIsRUFDaEMsTUFBTSxDQUFDLHlCQUF5QixFQUNoQyxNQUFNLENBQUMsMkJBQTJCLEVBQ2xDLE1BQU0sQ0FBQyxxQkFBcUIsQ0FDN0IsQ0FBQztTQUNIO1FBRUQsSUFDRSxJQUFJLENBQUMsT0FBTyxDQUFDLEtBQUssQ0FBQyxvQkFBb0I7WUFDdkMsSUFBSSxDQUFDLHVCQUF1QixFQUM1QjtZQUNBLEtBQUssTUFBTSxlQUFlLElBQUksTUFBTSxDQUFDLE1BQU0sQ0FBQyxJQUFJLENBQUMsdUJBQXVCLENBQUMsRUFBRTtnQkFDekUsSUFBSSxDQUFDLFNBQVMsQ0FDWixPQUFPLEVBQ1AsZUFBZSxFQUNmLE1BQU0sQ0FBQyw2QkFBNkIsRUFDcEMsTUFBTSxDQUFDLGlDQUFpQyxFQUN4QyxNQUFNLENBQUMsMkJBQTJCLENBQ25DLENBQUM7YUFDSDtTQUNGO1FBRUQsT0FBTyxDQUFDLE9BQU8sRUFBRSxDQUFDO0lBQ3BCLENBQUM7SUFFTyxvQkFBb0IsQ0FDMUIsT0FBaUMsRUFDakMsUUFBYSxFQUNiLE9BQWUsRUFDZixPQUFlLEVBQ2YsU0FBaUIsRUFDakIsSUFBWTtRQUVaLE9BQU8sQ0FBQyxJQUFJLEVBQUUsQ0FBQztRQUVmLE9BQU8sQ0FBQyxTQUFTLEdBQUcsU0FBUyxDQUFDO1FBRTlCLE9BQU8sQ0FBQyxXQUFXLEdBQUcsT0FBTyxDQUFDO1FBQzlCLE9BQU8sQ0FBQyxTQUFTLEVBQUUsQ0FBQztRQUNwQixPQUFPLENBQUMsTUFBTSxDQUFDLFFBQVEsQ0FBQyxDQUFDLEVBQUUsUUFBUSxDQUFDLENBQUMsQ0FBQyxDQUFDO1FBQ3ZDLE9BQU8sQ0FBQyxNQUFNLENBQUMsUUFBUSxDQUFDLENBQUMsR0FBRyxJQUFJLEVBQUUsUUFBUSxDQUFDLENBQUMsQ0FBQyxDQUFDO1FBQzlDLE9BQU8sQ0FBQyxNQUFNLEVBQUUsQ0FBQztRQUVqQixPQUFPLENBQUMsV0FBVyxHQUFHLE9BQU8sQ0FBQztRQUM5QixPQUFPLENBQUMsU0FBUyxFQUFFLENBQUM7UUFDcEIsT0FBTyxDQUFDLE1BQU0sQ0FBQyxRQUFRLENBQUMsQ0FBQyxFQUFFLFFBQVEsQ0FBQyxDQUFDLENBQUMsQ0FBQztRQUN2QyxPQUFPLENBQUMsTUFBTSxDQUFDLFFBQVEsQ0FBQyxDQUFDLEVBQUUsUUFBUSxDQUFDLENBQUMsR0FBRyxJQUFJLENBQUMsQ0FBQztRQUM5QyxPQUFPLENBQUMsTUFBTSxFQUFFLENBQUM7UUFFakIsT0FBTyxDQUFDLE9BQU8sRUFBRSxDQUFDO0lBQ3BCLENBQUM7SUFFTyxTQUFTLENBQ2YsT0FBaUMsRUFDakMsUUFBYSxFQUNiLE1BQWMsRUFDZCxTQUFpQixFQUNqQixJQUFZO1FBRVosT0FBTyxDQUFDLElBQUksRUFBRSxDQUFDO1FBRWYsT0FBTyxDQUFDLFNBQVMsR0FBRyxTQUFTLENBQUM7UUFFOUIsTUFBTSxRQUFRLEdBQUcsSUFBSSxDQUFDLElBQUksQ0FBQyxJQUFJLEdBQUcsQ0FBQyxDQUFDLENBQUM7UUFDckMsT0FBTyxDQUFDLFdBQVcsR0FBRyxNQUFNLENBQUM7UUFDN0IsT0FBTyxDQUFDLFNBQVMsRUFBRSxDQUFDO1FBQ3BCLE9BQU8sQ0FBQyxNQUFNLENBQUMsUUFBUSxDQUFDLENBQUMsR0FBRyxRQUFRLEVBQUUsUUFBUSxDQUFDLENBQUMsR0FBRyxRQUFRLENBQUMsQ0FBQztRQUM3RCxPQUFPLENBQUMsTUFBTSxDQUFDLFFBQVEsQ0FBQyxDQUFDLEdBQUcsUUFBUSxFQUFFLFFBQVEsQ0FBQyxDQUFDLEdBQUcsUUFBUSxDQUFDLENBQUM7UUFDN0QsT0FBTyxDQUFDLE1BQU0sQ0FBQyxRQUFRLENBQUMsQ0FBQyxHQUFHLFFBQVEsRUFBRSxRQUFRLENBQUMsQ0FBQyxHQUFHLFFBQVEsQ0FBQyxDQUFDO1FBQzdELE9BQU8sQ0FBQyxNQUFNLENBQUMsUUFBUSxDQUFDLENBQUMsR0FBRyxRQUFRLEVBQUUsUUFBUSxDQUFDLENBQUMsR0FBRyxRQUFRLENBQUMsQ0FBQztRQUM3RCxPQUFPLENBQUMsTUFBTSxFQUFFLENBQUM7UUFFakIsT0FBTyxDQUFDLE9BQU8sRUFBRSxDQUFDO0lBQ3BCLENBQUM7O0FBcGRILHdCQXFkQztBQXBkeUIsc0JBQWUsR0FBa0I7SUFDdkQsVUFBVSxFQUFFLENBQUMsU0FBUyxDQUFDO0lBQ3ZCLGdCQUFnQixFQUFFLFNBQVM7SUFDM0IsVUFBVSxFQUFFO1FBQ1YsT0FBTyxFQUFFO1lBQ1AsR0FBRyxFQUFFO2dCQUNILElBQUksRUFBRSxTQUFTO2dCQUNmLFVBQVUsRUFBRSxDQUFDO2dCQUNiLFNBQVMsRUFBRSxDQUFDO2dCQUNaLElBQUksRUFBRSx5QkFBeUIsQ0FBQyxlQUFlO2FBQ2hEO1NBQ0Y7S0FDRjtJQUNELGdCQUFnQixFQUFFLFNBQVM7Q0FDNUIsQ0FBQztBQUVzQixnQ0FBeUIsR0FBRyxPQUFPLENBQUM7QUFDcEMsb0NBQTZCLEdBQUcsQ0FBQyxDQUFDO0FBRWxDLGdDQUF5QixHQUFHLEtBQUssQ0FBQztBQUNsQyxnQ0FBeUIsR0FBRyxRQUFRLENBQUM7QUFDckMsa0NBQTJCLEdBQUcsQ0FBQyxDQUFDO0FBQ2hDLDRCQUFxQixHQUFHLEVBQUUsQ0FBQztBQUUzQixvQ0FBNkIsR0FBRyxNQUFNLENBQUM7QUFDdkMsd0NBQWlDLEdBQUcsQ0FBQyxDQUFDO0FBQ3RDLGtDQUEyQixHQUFHLENBQUMsQ0FBQyJ9

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __webpack_exports__ = __webpack_require__("./index.ts");
/******/ 	
/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,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Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGtEQUFrRCxnQkFBZ0I7QUFDbEU7QUFDQTtBQUNBO0FBQ0Esc0VBQXNFLGdCQUFnQjtBQUN0RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsY0FBYztBQUNkO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYTtBQUNiLFNBQVM7QUFDVCxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLDJDQUEyQzs7Ozs7O1VDdlUzQztVQUNBOztVQUVBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBOztVQUVBO1VBQ0E7O1VBRUE7VUFDQTtVQUNBOzs7O1VFdEJBO1VBQ0E7VUFDQTtVQUNBIiwic291cmNlcyI6WyJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2Uvc3ByaXRlL3dlYnBhY2svdW5pdmVyc2FsTW9kdWxlRGVmaW5pdGlvbiIsIndlYnBhY2s6Ly9AYmFzZW1lbnR1bml2ZXJzZS9zcHJpdGUvLi9ub2RlX21vZHVsZXMvQGJhc2VtZW50dW5pdmVyc2UvdXRpbHMvdXRpbHMuanMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2Uvc3ByaXRlLy4vbm9kZV9tb2R1bGVzL0BiYXNlbWVudHVuaXZlcnNlL3ZlYy92ZWMuanMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2Uvc3ByaXRlLy4vaW5kZXgudHMiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2Uvc3ByaXRlL3dlYnBhY2svYm9vdHN0cmFwIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3Nwcml0ZS93ZWJwYWNrL2JlZm9yZS1zdGFydHVwIiwid2VicGFjazovL0BiYXNlbWVudHVuaXZlcnNlL3Nwcml0ZS93ZWJwYWNrL3N0YXJ0dXAiLCJ3ZWJwYWNrOi8vQGJhc2VtZW50dW5pdmVyc2Uvc3ByaXRlL3dlYnBhY2svYWZ0ZXItc3RhcnR1cCJdLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gd2VicGFja1VuaXZlcnNhbE1vZHVsZURlZmluaXRpb24ocm9vdCwgZmFjdG9yeSkge1xuXHRpZih0eXBlb2YgZXhwb3J0cyA9PT0gJ29iamVjdCcgJiYgdHlwZW9mIG1vZHVsZSA9PT0gJ29iamVjdCcpXG5cdFx0bW9kdWxlLmV4cG9ydHMgPSBmYWN0b3J5KCk7XG5cdGVsc2UgaWYodHlwZW9mIGRlZmluZSA9PT0gJ2Z1bmN0aW9uJyAmJiBkZWZpbmUuYW1kKVxuXHRcdGRlZmluZShbXSwgZmFjdG9yeSk7XG5cdGVsc2Uge1xuXHRcdHZhciBhID0gZmFjdG9yeSgpO1xuXHRcdGZvcih2YXIgaSBpbiBhKSAodHlwZW9mIGV4cG9ydHMgPT09ICdvYmplY3QnID8gZXhwb3J0cyA6IHJvb3QpW2ldID0gYVtpXTtcblx0fVxufSkoc2VsZiwgKCkgPT4ge1xucmV0dXJuICIsIi8qKlxuICogQG92ZXJ2aWV3IEEgbGlicmFyeSBvZiB1c2VmdWwgZnVuY3Rpb25zXG4gKiBAYXV0aG9yIEdvcmRvbiBMYXJyaWdhblxuICovXG5cbi8qKlxuICogQ2hlY2sgaWYgdHdvIG51bWJlcnMgYXJlIGFwcHJveGltYXRlbHkgZXF1YWxcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIE51bWJlciBhXG4gKiBAcGFyYW0ge251bWJlcn0gYiBOdW1iZXIgYlxuICogQHBhcmFtIHtudW1iZXJ9IFtwPU51bWJlci5FUFNJTE9OXSBUaGUgcHJlY2lzaW9uIHZhbHVlXG4gKiBAcmV0dXJuIHtib29sZWFufSBUcnVlIGlmIG51bWJlcnMgYSBhbmQgYiBhcmUgYXBwcm94aW1hdGVseSBlcXVhbFxuICovXG5jb25zdCBmbG9hdEVxdWFscyA9IChhLCBiLCBwID0gTnVtYmVyLkVQU0lMT04pID0+IE1hdGguYWJzKGEgLSBiKSA8IHA7XG5cbi8qKlxuICogQ2xhbXAgYSBudW1iZXIgYmV0d2VlbiBtaW4gYW5kIG1heFxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG51bWJlciB0byBjbGFtcFxuICogQHBhcmFtIHtudW1iZXJ9IFttaW49MF0gVGhlIG1pbmltdW0gdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBbbWF4PTFdIFRoZSBtYXhpbXVtIHZhbHVlXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgY2xhbXBlZCBudW1iZXJcbiAqL1xuY29uc3QgY2xhbXAgPSAoYSwgbWluID0gMCwgbWF4ID0gMSkgPT4gYSA8IG1pbiA/IG1pbiA6IChhID4gbWF4ID8gbWF4IDogYSk7XG5cbi8qKlxuICogR2V0IHRoZSBmcmFjdGlvbmFsIHBhcnQgb2YgYSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBudW1iZXIgZnJvbSB3aGljaCB0byBnZXQgdGhlIGZyYWN0aW9uYWwgcGFydFxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgZnJhY3Rpb25hbCBwYXJ0IG9mIHRoZSBudW1iZXJcbiAqL1xuY29uc3QgZnJhYyA9IGEgPT4gYSA+PSAwID8gYSAtIE1hdGguZmxvb3IoYSkgOiBhIC0gTWF0aC5jZWlsKGEpO1xuXG4vKipcbiAqIFJvdW5kIG4gdG8gZCBkZWNpbWFsIHBsYWNlc1xuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIG51bWJlciB0byByb3VuZFxuICogQHBhcmFtIHtudW1iZXJ9IFtkPTBdIFRoZSBudW1iZXIgb2YgZGVjaW1hbCBwbGFjZXMgdG8gcm91bmQgdG9cbiAqIEByZXR1cm4ge251bWJlcn0gQSByb3VuZGVkIG51bWJlclxuICovXG5jb25zdCByb3VuZCA9IChuLCBkID0gMCkgPT4ge1xuICBjb25zdCBwID0gTWF0aC5wb3coMTAsIGQpO1xuICByZXR1cm4gTWF0aC5yb3VuZChuICogcCArIE51bWJlci5FUFNJTE9OKSAvIHA7XG59XG5cbi8qKlxuICogRG8gYSBsaW5lYXIgaW50ZXJwb2xhdGlvbiBiZXR3ZWVuIGEgYW5kIGJcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBtaW5pbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGIgVGhlIG1heGltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgaW50ZXJwb2xhdGlvbiB2YWx1ZSwgc2hvdWxkIGJlIGluIHRoZSBpbnRlcnZhbCBbMCwgMV1cbiAqIEByZXR1cm4ge251bWJlcn0gQW4gaW50ZXJwb2xhdGVkIHZhbHVlIGluIHRoZSBpbnRlcnZhbCBbYSwgYl1cbiAqL1xuY29uc3QgbGVycCA9IChhLCBiLCBpKSA9PiBhICsgKGIgLSBhKSAqIGk7XG5cbi8qKlxuICogR2V0IHRoZSBwb3NpdGlvbiBvZiBpIGJldHdlZW4gYSBhbmQgYlxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG1pbmltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYiBUaGUgbWF4aW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBpbnRlcnBvbGF0ZWQgdmFsdWUgaW4gdGhlIGludGVydmFsIFthLCBiXVxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgcG9zaXRpb24gb2YgaSBiZXR3ZWVuIGEgYW5kIGJcbiAqL1xuY29uc3QgdW5sZXJwID0gKGEsIGIsIGkpID0+IChpIC0gYSkgLyAoYiAtIGEpO1xuXG4vKipcbiAqIERvIGEgYmlsaW5lYXIgaW50ZXJwb2xhdGlvblxuICogQHBhcmFtIHtudW1iZXJ9IGMwMCBUb3AtbGVmdCB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IGMxMCBUb3AtcmlnaHQgdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBjMDEgQm90dG9tLWxlZnQgdmFsdWVcbiAqIEBwYXJhbSB7bnVtYmVyfSBjMTEgQm90dG9tLXJpZ2h0IHZhbHVlXG4gKiBAcGFyYW0ge251bWJlcn0gaXggSW50ZXJwb2xhdGlvbiB2YWx1ZSBhbG9uZyB4XG4gKiBAcGFyYW0ge251bWJlcn0gaXkgSW50ZXJwb2xhdGlvbiB2YWx1ZSBhbG9uZyB5XG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgYmlsaW5lYXIgaW50ZXJwb2xhdGVkIHZhbHVlXG4gKi9cbmNvbnN0IGJsZXJwID0gKGMwMCwgYzEwLCBjMDEsIGMxMSwgaXgsIGl5KSA9PiBsZXJwKGxlcnAoYzAwLCBjMTAsIGl4KSwgbGVycChjMDEsIGMxMSwgaXgpLCBpeSk7XG5cbi8qKlxuICogUmUtbWFwIGEgbnVtYmVyIGkgZnJvbSByYW5nZSBhMS4uLmEyIHRvIGIxLi4uYjJcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBudW1iZXIgdG8gcmUtbWFwXG4gKiBAcGFyYW0ge251bWJlcn0gYTFcbiAqIEBwYXJhbSB7bnVtYmVyfSBhMlxuICogQHBhcmFtIHtudW1iZXJ9IGIxXG4gKiBAcGFyYW0ge251bWJlcn0gYjJcbiAqIEByZXR1cm4ge251bWJlcn1cbiAqL1xuY29uc3QgcmVtYXAgPSAoaSwgYTEsIGEyLCBiMSwgYjIpID0+IGIxICsgKGkgLSBhMSkgKiAoYjIgLSBiMSkgLyAoYTIgLSBhMSk7XG5cbi8qKlxuICogRG8gYSBzbW9vdGggaW50ZXJwb2xhdGlvbiBiZXR3ZWVuIGEgYW5kIGJcbiAqIEBwYXJhbSB7bnVtYmVyfSBhIFRoZSBtaW5pbXVtIG51bWJlclxuICogQHBhcmFtIHtudW1iZXJ9IGIgVGhlIG1heGltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgaW50ZXJwb2xhdGlvbiB2YWx1ZVxuICogQHJldHVybiB7bnVtYmVyfSBBbiBpbnRlcnBvbGF0ZWQgdmFsdWUgaW4gdGhlIGludGVydmFsIFthLCBiXVxuICovXG5jb25zdCBzbW9vdGhzdGVwID0gKGEsIGIsIGkpID0+IGxlcnAoYSwgYiwgMyAqIE1hdGgucG93KGksIDIpIC0gMiAqIE1hdGgucG93KGksIDMpKTtcblxuLyoqXG4gKiBHZXQgYW4gYW5nbGUgaW4gcmFkaWFuc1xuICogQHBhcmFtIHtudW1iZXJ9IGRlZ3JlZXMgVGhlIGFuZ2xlIGluIGRlZ3JlZXNcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIGFuZ2xlIGluIHJhZGlhbnNcbiAqL1xuY29uc3QgcmFkaWFucyA9IGRlZ3JlZXMgPT4gKE1hdGguUEkgLyAxODApICogZGVncmVlcztcblxuLyoqXG4gKiBHZXQgYW4gYW5nbGUgaW4gZGVncmVlc1xuICogQHBhcmFtIHtudW1iZXJ9IHJhZGlhbnMgVGhlIGFuZ2xlIGluIHJhZGlhbnNcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIGFuZ2xlIGluIGRlZ3JlZXNcbiAqL1xuY29uc3QgZGVncmVlcyA9IHJhZGlhbnMgPT4gKDE4MCAvIE1hdGguUEkpICogcmFkaWFucztcblxuLyoqXG4gKiBHZXQgYSByYW5kb20gZmxvYXQgaW4gdGhlIGludGVydmFsIFttaW4sIG1heClcbiAqIEBwYXJhbSB7bnVtYmVyfSBtaW4gSW5jbHVzaXZlIG1pblxuICogQHBhcmFtIHtudW1iZXJ9IG1heCBFeGNsdXNpdmUgbWF4XG4gKiBAcmV0dXJuIHtudW1iZXJ9IEEgcmFuZG9tIGZsb2F0IGluIHRoZSBpbnRlcnZhbCBbbWluLCBtYXgpXG4gKi9cbmNvbnN0IHJhbmRvbUJldHdlZW4gPSAobWluLCBtYXgpID0+IE1hdGgucmFuZG9tKCkgKiAobWF4IC0gbWluKSArIG1pbjtcblxuLyoqXG4gKiBHZXQgYSByYW5kb20gaW50ZWdlciBpbiB0aGUgaW50ZXJ2YWwgW21pbiwgbWF4XVxuICogQHBhcmFtIHtudW1iZXJ9IG1pbiBJbmNsdXNpdmUgbWluXG4gKiBAcGFyYW0ge251bWJlcn0gbWF4IEluY2x1c2l2ZSBtYXhcbiAqIEByZXR1cm4ge251bWJlcn0gQSByYW5kb20gaW50ZWdlciBpbiB0aGUgaW50ZXJ2YWwgW21pbiwgbWF4XVxuICovXG5jb25zdCByYW5kb21JbnRCZXR3ZWVuID0gKG1pbiwgbWF4KSA9PiBNYXRoLmZsb29yKE1hdGgucmFuZG9tKCkgKiAobWF4IC0gbWluICsgMSkpICsgbWluO1xuXG4vKipcbiAqIEdldCBhIG5vcm1hbGx5LWRpc3RyaWJ1dGVkIHJhbmRvbSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBbbXU9MC41XSBUaGUgbWVhbiB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IFtzaWdtYT0wLjVdIFRoZSBzdGFuZGFyZCBkZXZpYXRpb25cbiAqIEBwYXJhbSB7bnVtYmVyfSBbc2FtcGxlcz0yXSBUaGUgbnVtYmVyIG9mIHNhbXBsZXNcbiAqIEByZXR1cm4ge251bWJlcn0gQSBub3JtYWxseS1kaXN0cmlidXRlZCByYW5kb20gbnVtYmVyXG4gKi9cbmNvbnN0IGNsdFJhbmRvbSA9IChtdSA9IDAuNSwgc2lnbWEgPSAwLjUsIHNhbXBsZXMgPSAyKSA9PiB7XG4gIGxldCB0b3RhbCA9IDA7XG4gIGZvciAobGV0IGkgPSBzYW1wbGVzOyBpLS07KSB7XG4gICAgdG90YWwgKz0gTWF0aC5yYW5kb20oKTtcbiAgfVxuICByZXR1cm4gbXUgKyAodG90YWwgLSBzYW1wbGVzIC8gMikgLyAoc2FtcGxlcyAvIDIpICogc2lnbWE7XG59O1xuXG4vKipcbiAqIEdldCBhIG5vcm1hbGx5LWRpc3RyaWJ1dGVkIHJhbmRvbSBpbnRlZ2VyIGluIHRoZSBpbnRlcnZhbCBbbWluLCBtYXhdXG4gKiBAcGFyYW0ge251bWJlcn0gbWluIEluY2x1c2l2ZSBtaW5cbiAqIEBwYXJhbSB7bnVtYmVyfSBtYXggSW5jbHVzaXZlIG1heFxuICogQHJldHVybiB7bnVtYmVyfSBBIG5vcm1hbGx5LWRpc3RyaWJ1dGVkIHJhbmRvbSBpbnRlZ2VyXG4gKi9cbmNvbnN0IGNsdFJhbmRvbUludCA9IChtaW4sIG1heCkgPT4gTWF0aC5mbG9vcihtaW4gKyBjbHRSYW5kb20oMC41LCAwLjUsIDIpICogKG1heCArIDEgLSBtaW4pKTtcblxuLyoqXG4gKiBSZXR1cm4gYSB3ZWlnaHRlZCByYW5kb20gaW50ZWdlclxuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSB3IEFuIGFycmF5IG9mIHdlaWdodHNcbiAqIEByZXR1cm4ge251bWJlcn0gQW4gaW5kZXggZnJvbSB3XG4gKi9cbmNvbnN0IHdlaWdodGVkUmFuZG9tID0gdyA9PiB7XG4gIGxldCB0b3RhbCA9IHcucmVkdWNlKChhLCBpKSA9PiBhICsgaSwgMCksIG4gPSAwO1xuICBjb25zdCByID0gTWF0aC5yYW5kb20oKSAqIHRvdGFsO1xuICB3aGlsZSAodG90YWwgPiByKSB7XG4gICAgdG90YWwgLT0gd1tuKytdO1xuICB9XG4gIHJldHVybiBuIC0gMTtcbn07XG5cbi8qKlxuICogQW4gaW50ZXJwb2xhdGlvbiBmdW5jdGlvblxuICogQGNhbGxiYWNrIEludGVycG9sYXRpb25GdW5jdGlvblxuICogQHBhcmFtIHtudW1iZXJ9IGEgVGhlIG1pbmltdW0gbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYiBUaGUgbWF4aW11bSBudW1iZXJcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBpbnRlcnBvbGF0aW9uIHZhbHVlLCBzaG91bGQgYmUgaW4gdGhlIGludGVydmFsIFswLCAxXVxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgaW50ZXJwb2xhdGVkIHZhbHVlIGluIHRoZSBpbnRlcnZhbCBbYSwgYl1cbiAqL1xuXG4vKipcbiAqIFJldHVybiBhbiBpbnRlcnBvbGF0ZWQgdmFsdWUgZnJvbSBhbiBhcnJheVxuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSBhIEFuIGFycmF5IG9mIHZhbHVlcyBpbnRlcnBvbGF0ZVxuICogQHBhcmFtIHtudW1iZXJ9IGkgQSBudW1iZXIgaW4gdGhlIGludGVydmFsIFswLCAxXVxuICogQHBhcmFtIHtJbnRlcnBvbGF0aW9uRnVuY3Rpb259IFtmPU1hdGgubGVycF0gVGhlIGludGVycG9sYXRpb24gZnVuY3Rpb24gdG8gdXNlXG4gKiBAcmV0dXJuIHtudW1iZXJ9IEFuIGludGVycG9sYXRlZCB2YWx1ZSBpbiB0aGUgaW50ZXJ2YWwgW21pbihhKSwgbWF4KGEpXVxuICovXG5jb25zdCBsZXJwQXJyYXkgPSAoYSwgaSwgZiA9IGxlcnApID0+IHtcbiAgY29uc3QgcyA9IGkgKiAoYS5sZW5ndGggLSAxKTtcbiAgY29uc3QgcCA9IGNsYW1wKE1hdGgudHJ1bmMocyksIDAsIGEubGVuZ3RoIC0gMSk7XG4gIHJldHVybiBmKGFbcF0gfHwgMCwgYVtwICsgMV0gfHwgMCwgZnJhYyhzKSk7XG59O1xuXG4vKipcbiAqIEdldCB0aGUgZG90IHByb2R1Y3Qgb2YgdHdvIHZlY3RvcnNcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gYSBWZWN0b3IgYVxuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSBiIFZlY3RvciBiXG4gKiBAcmV0dXJuIHtudW1iZXJ9IGEg4oiZIGJcbiAqL1xuY29uc3QgZG90ID0gKGEsIGIpID0+IGEucmVkdWNlKChuLCB2LCBpKSA9PiBuICsgdiAqIGJbaV0sIDApO1xuXG4vKipcbiAqIEdldCB0aGUgZmFjdG9yaWFsIG9mIGEgbnVtYmVyXG4gKiBAcGFyYW0ge251bWJlcn0gYVxuICogQHJldHVybiB7bnVtYmVyfSBhIVxuICovXG5jb25zdCBmYWN0b3JpYWwgPSBhID0+IHtcbiAgbGV0IHJlc3VsdCA9IDE7XG4gIGZvciAobGV0IGkgPSAyOyBpIDw9IGE7IGkrKykge1xuICAgIHJlc3VsdCAqPSBpO1xuICB9XG4gIHJldHVybiByZXN1bHQ7XG59O1xuXG4vKipcbiAqIEdldCB0aGUgbnVtYmVyIG9mIHBlcm11dGF0aW9ucyBvZiByIGVsZW1lbnRzIGZyb20gYSBzZXQgb2YgbiBlbGVtZW50c1xuICogQHBhcmFtIHtudW1iZXJ9IG5cbiAqIEBwYXJhbSB7bnVtYmVyfSByXG4gKiBAcmV0dXJuIHtudW1iZXJ9IG5QclxuICovXG5jb25zdCBucHIgPSAobiwgcikgPT4gZmFjdG9yaWFsKG4pIC8gZmFjdG9yaWFsKG4gLSByKTtcblxuLyoqXG4gKiBHZXQgdGhlIG51bWJlciBvZiBjb21iaW5hdGlvbnMgb2YgciBlbGVtZW50cyBmcm9tIGEgc2V0IG9mIG4gZWxlbWVudHNcbiAqIEBwYXJhbSB7bnVtYmVyfSBuXG4gKiBAcGFyYW0ge251bWJlcn0gclxuICogQHJldHVybiB7bnVtYmVyfSBuQ3JcbiAqL1xuY29uc3QgbmNyID0gKG4sIHIpID0+IGZhY3RvcmlhbChuKSAvIChmYWN0b3JpYWwocikgKiBmYWN0b3JpYWwobiAtIHIpKTtcblxuLyoqXG4gKiBHZW5lcmF0ZSBhbGwgY29tYmluYXRpb25zIG9mIHIgZWxlbWVudHMgZnJvbSBhbiBhcnJheVxuICpcbiAqIEBleGFtcGxlXG4gKiBgYGBqc1xuICogY29tYmluYXRpb25zKFsxLCAyLCAzXSwgMik7XG4gKiBgYGBcbiAqXG4gKiBPdXRwdXQ6XG4gKiBgYGBqc29uXG4gKiBbXG4gKiAgIFsxLCAyXSxcbiAqICAgWzEsIDNdLFxuICogICBbMiwgM11cbiAqIF1cbiAqIGBgYFxuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHBhcmFtIHtudW1iZXJ9IHIgVGhlIG51bWJlciBvZiBlbGVtZW50cyB0byBjaG9vc2UgaW4gZWFjaCBjb21iaW5hdGlvblxuICogQHJldHVybiB7QXJyYXk8QXJyYXk8Kj4+fSBBbiBhcnJheSBvZiBjb21iaW5hdGlvbiBhcnJheXNcbiAqL1xuY29uc3QgY29tYmluYXRpb25zID0gKGEsIHIpID0+IHtcbiAgaWYgKHIgPT09IDEpIHtcbiAgICByZXR1cm4gYS5tYXAoaXRlbSA9PiBbaXRlbV0pO1xuICB9XG5cbiAgcmV0dXJuIGEucmVkdWNlKFxuICAgIChhY2MsIGl0ZW0sIGkpID0+IFtcbiAgICAgIC4uLmFjYyxcbiAgICAgIC4uLmNvbWJpbmF0aW9ucyhhLnNsaWNlKGkgKyAxKSwgciAtIDEpLm1hcChjID0+IFtpdGVtLCAuLi5jXSksXG4gICAgXSxcbiAgICBbXVxuICApO1xufTtcblxuLyoqXG4gKiBHZXQgYSBjYXJ0ZXNpYW4gcHJvZHVjdCBvZiBhcnJheXNcbiAqXG4gKiBAZXhhbXBsZVxuICogYGBganNcbiAqIGNhcnRlc2lhbihbMSwgMiwgM10sIFsnYScsICdiJ10pO1xuICogYGBgXG4gKlxuICogT3V0cHV0OlxuICogYGBganNvblxuICogW1xuICogICBbMSwgXCJhXCJdLFxuICogICBbMSwgXCJiXCJdLFxuICogICBbMiwgXCJhXCJdLFxuICogICBbMiwgXCJiXCJdLFxuICogICBbMywgXCJhXCJdLFxuICogICBbMywgXCJiXCJdXG4gKiBdXG4gKiBgYGBcbiAqL1xuY29uc3QgY2FydGVzaWFuID0gKC4uLmFycikgPT5cbiAgYXJyLnJlZHVjZShcbiAgICAoYSwgYikgPT4gYS5mbGF0TWFwKGMgPT4gYi5tYXAoZCA9PiBbLi4uYywgZF0pKSxcbiAgICBbW11dXG4gICk7XG5cbi8qKlxuICogQSBmdW5jdGlvbiBmb3IgZ2VuZXJhdGluZyBhcnJheSB2YWx1ZXNcbiAqIEBjYWxsYmFjayBUaW1lc0Z1bmN0aW9uXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgYXJyYXkgaW5kZXhcbiAqIEByZXR1cm4geyp9IFRoZSBhcnJheSB2YWx1ZVxuICovXG5cbi8qKlxuICogUmV0dXJuIGEgbmV3IGFycmF5IHdpdGggbGVuZ3RoIG4gYnkgY2FsbGluZyBmdW5jdGlvbiBmKGkpIG9uIGVhY2ggZWxlbWVudFxuICogQHBhcmFtIHtUaW1lc0Z1bmN0aW9ufSBmXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgc2l6ZSBvZiB0aGUgYXJyYXlcbiAqIEByZXR1cm4ge0FycmF5PCo+fVxuICovXG5jb25zdCB0aW1lcyA9IChmLCBuKSA9PiBBcnJheShuKS5maWxsKDApLm1hcCgoXywgaSkgPT4gZihpKSk7XG5cbi8qKlxuICogUmV0dXJuIGFuIGFycmF5IGNvbnRhaW5pbmcgbnVtYmVycyAwLT4obiAtIDEpXG4gKiBAcGFyYW0ge251bWJlcn0gbiBUaGUgc2l6ZSBvZiB0aGUgYXJyYXlcbiAqIEByZXR1cm4ge0FycmF5PG51bWJlcj59IEFuIGFycmF5IG9mIGludGVnZXJzIDAtPihuIC0gMSlcbiAqL1xuY29uc3QgcmFuZ2UgPSBuID0+IHRpbWVzKGkgPT4gaSwgbik7XG5cbi8qKlxuICogWmlwIDIgYXJyYXlzIHRvZ2V0aGVyLCBpLmUuIChbMSwgMiwgM10sIFthLCBiLCBjXSkgPT4gW1sxLCBhXSwgWzIsIGJdLCBbMywgY11dXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBiXG4gKiBAcmV0dXJuIHtBcnJheTxBcnJheTwqPj59XG4gKi9cbmNvbnN0IHppcCA9IChhLCBiKSA9PiBhLm1hcCgoaywgaSkgPT4gW2ssIGJbaV1dKTtcblxuLyoqXG4gKiBSZXR1cm4gYXJyYXlbaV0gd2l0aCBwb3NpdGl2ZSBhbmQgbmVnYXRpdmUgd3JhcHBpbmdcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSBwb3NpdGl2ZWx5L25lZ2F0aXZlbHkgd3JhcHBlZCBhcnJheSBpbmRleFxuICogQHJldHVybiB7Kn0gQW4gZWxlbWVudCBmcm9tIHRoZSBhcnJheVxuICovXG5jb25zdCBhdCA9IChhLCBpKSA9PiBhW2kgPCAwID8gYS5sZW5ndGggLSAoTWF0aC5hYnMoaSArIDEpICUgYS5sZW5ndGgpIC0gMSA6IGkgJSBhLmxlbmd0aF07XG5cbi8qKlxuICogUmV0dXJuIHRoZSBsYXN0IGVsZW1lbnQgb2YgYW4gYXJyYXkgd2l0aG91dCByZW1vdmluZyBpdFxuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHJldHVybiB7Kn0gVGhlIGxhc3QgZWxlbWVudCBmcm9tIHRoZSBhcnJheVxuICovXG5jb25zdCBwZWVrID0gKGEpID0+IHtcbiAgaWYgKCFhLmxlbmd0aCkge1xuICAgIHJldHVybiB1bmRlZmluZWQ7XG4gIH1cblxuICByZXR1cm4gYVthLmxlbmd0aCAtIDFdO1xufTtcblxuLyoqXG4gKiBDaG9wIGFuIGFycmF5IGludG8gY2h1bmtzIG9mIHNpemUgblxuICogQHBhcmFtIHtBcnJheTwqPn0gYVxuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIGNodW5rIHNpemVcbiAqIEByZXR1cm4ge0FycmF5PEFycmF5PCo+Pn0gQW4gYXJyYXkgb2YgYXJyYXkgY2h1bmtzXG4gKi9cbmNvbnN0IGNodW5rID0gKGEsIG4pID0+IHRpbWVzKGkgPT4gYS5zbGljZShpICogbiwgaSAqIG4gKyBuKSwgTWF0aC5jZWlsKGEubGVuZ3RoIC8gbikpO1xuXG4vKipcbiAqIFJhbmRvbWx5IHNodWZmbGUgYSBzaGFsbG93IGNvcHkgb2YgYW4gYXJyYXlcbiAqIEBwYXJhbSB7QXJyYXk8Kj59IGFcbiAqIEByZXR1cm4ge0FycmF5PCo+fSBUaGUgc2h1ZmZsZWQgYXJyYXlcbiAqL1xuY29uc3Qgc2h1ZmZsZSA9IGEgPT4gYS5zbGljZSgpLnNvcnQoKCkgPT4gTWF0aC5yYW5kb20oKSAtIDAuNSk7XG5cbi8qKlxuICogRmxhdHRlbiBhbiBvYmplY3RcbiAqIEBwYXJhbSB7b2JqZWN0fSBvXG4gKiBAcGFyYW0ge3N0cmluZ30gY29uY2F0ZW5hdG9yIFRoZSBzdHJpbmcgdG8gdXNlIGZvciBjb25jYXRlbmF0aW5nIGtleXNcbiAqIEByZXR1cm4ge29iamVjdH0gQSBmbGF0dGVuZWQgb2JqZWN0XG4gKi9cbmNvbnN0IGZsYXQgPSAobywgY29uY2F0ZW5hdG9yID0gJy4nKSA9PiB7XG4gIHJldHVybiBPYmplY3Qua2V5cyhvKS5yZWR1Y2UoKGFjYywga2V5KSA9PiB7XG4gICAgaWYgKG9ba2V5XSBpbnN0YW5jZW9mIERhdGUpIHtcbiAgICAgIHJldHVybiB7XG4gICAgICAgIC4uLmFjYyxcbiAgICAgICAgW2tleV06IG9ba2V5XS50b0lTT1N0cmluZygpLFxuICAgICAgfTtcbiAgICB9XG5cbiAgICBpZiAodHlwZW9mIG9ba2V5XSAhPT0gJ29iamVjdCcgfHwgIW9ba2V5XSkge1xuICAgICAgcmV0dXJuIHtcbiAgICAgICAgLi4uYWNjLFxuICAgICAgICBba2V5XTogb1trZXldLFxuICAgICAgfTtcbiAgICB9XG4gICAgY29uc3QgZmxhdHRlbmVkID0gZmxhdChvW2tleV0sIGNvbmNhdGVuYXRvcik7XG5cbiAgICByZXR1cm4ge1xuICAgICAgLi4uYWNjLFxuICAgICAgLi4uT2JqZWN0LmtleXMoZmxhdHRlbmVkKS5yZWR1Y2UoXG4gICAgICAgIChjaGlsZEFjYywgY2hpbGRLZXkpID0+ICh7XG4gICAgICAgICAgLi4uY2hpbGRBY2MsXG4gICAgICAgICAgW2Ake2tleX0ke2NvbmNhdGVuYXRvcn0ke2NoaWxkS2V5fWBdOiBmbGF0dGVuZWRbY2hpbGRLZXldLFxuICAgICAgICB9KSxcbiAgICAgICAge31cbiAgICAgICksXG4gICAgfTtcbiAgfSwge30pO1xufTtcblxuLyoqXG4gKiBVbmZsYXR0ZW4gYW4gb2JqZWN0XG4gKiBAcGFyYW0ge29iamVjdH0gb1xuICogQHBhcmFtIHtzdHJpbmd9IGNvbmNhdGVuYXRvciBUaGUgc3RyaW5nIHRvIGNoZWNrIGZvciBpbiBjb25jYXRlbmF0ZWQga2V5c1xuICogQHJldHVybiB7b2JqZWN0fSBBbiB1bi1mbGF0dGVuZWQgb2JqZWN0XG4gKi9cbmNvbnN0IHVuZmxhdCA9IChvLCBjb25jYXRlbmF0b3IgPSAnLicpID0+IHtcbiAgbGV0IHJlc3VsdCA9IHt9LCB0ZW1wLCBzdWJzdHJpbmdzLCBwcm9wZXJ0eSwgaTtcblxuICBmb3IgKHByb3BlcnR5IGluIG8pIHtcbiAgICBzdWJzdHJpbmdzID0gcHJvcGVydHkuc3BsaXQoY29uY2F0ZW5hdG9yKTtcbiAgICB0ZW1wID0gcmVzdWx0O1xuICAgIGZvciAoaSA9IDA7IGkgPCBzdWJzdHJpbmdzLmxlbmd0aCAtIDE7IGkrKykge1xuICAgICAgaWYgKCEoc3Vic3RyaW5nc1tpXSBpbiB0ZW1wKSkge1xuICAgICAgICBpZiAoaXNGaW5pdGUoc3Vic3RyaW5nc1tpICsgMV0pKSB7XG4gICAgICAgICAgdGVtcFtzdWJzdHJpbmdzW2ldXSA9IFtdO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHRlbXBbc3Vic3RyaW5nc1tpXV0gPSB7fTtcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdGVtcCA9IHRlbXBbc3Vic3RyaW5nc1tpXV07XG4gICAgfVxuICAgIHRlbXBbc3Vic3RyaW5nc1tzdWJzdHJpbmdzLmxlbmd0aCAtIDFdXSA9IG9bcHJvcGVydHldO1xuICB9XG5cbiAgcmV0dXJuIHJlc3VsdDtcbn07XG5cbi8qKlxuICogQSBzcGxpdCBwcmVkaWNhdGVcbiAqIEBjYWxsYmFjayBTcGxpdFByZWRpY2F0ZVxuICogQHBhcmFtIHthbnl9IHZhbHVlIFRoZSBjdXJyZW50IHZhbHVlXG4gKiBAcmV0dXJuIHtib29sZWFufSBUcnVlIGlmIHRoZSBhcnJheSBzaG91bGQgc3BsaXQgYXQgdGhpcyBpbmRleFxuICovXG5cbi8qKlxuICogU3BsaXQgYW4gYXJyYXkgaW50byBzdWItYXJyYXlzIGJhc2VkIG9uIGEgcHJlZGljYXRlXG4gKiBAcGFyYW0ge0FycmF5PCo+fSBhcnJheVxuICogQHBhcmFtIHtTcGxpdFByZWRpY2F0ZX0gcHJlZGljYXRlXG4gKiBAcmV0dXJuIHtBcnJheTxBcnJheTwqPj59IEFuIGFycmF5IG9mIGFycmF5c1xuICovXG5jb25zdCBzcGxpdCA9IChhcnJheSwgcHJlZGljYXRlKSA9PiB7XG4gIGNvbnN0IHJlc3VsdCA9IFtdO1xuICBsZXQgY3VycmVudCA9IFtdO1xuICBmb3IgKGNvbnN0IHZhbHVlIG9mIGFycmF5KSB7XG4gICAgaWYgKHByZWRpY2F0ZSh2YWx1ZSkpIHtcbiAgICAgIGlmIChjdXJyZW50Lmxlbmd0aCkge1xuICAgICAgICByZXN1bHQucHVzaChjdXJyZW50KTtcbiAgICAgIH1cbiAgICAgIGN1cnJlbnQgPSBbdmFsdWVdO1xuICAgIH0gZWxzZSB7XG4gICAgICBjdXJyZW50LnB1c2godmFsdWUpO1xuICAgIH1cbiAgfVxuICByZXN1bHQucHVzaChjdXJyZW50KTtcblxuICByZXR1cm4gcmVzdWx0O1xufTtcblxuLyoqXG4gKiBQbHVjayBrZXlzIGZyb20gYW4gb2JqZWN0XG4gKiBAcGFyYW0ge29iamVjdH0gb1xuICogQHBhcmFtIHsuLi5zdHJpbmd9IGtleXMgVGhlIGtleXMgdG8gcGx1Y2sgZnJvbSB0aGUgb2JqZWN0XG4gKiBAcmV0dXJuIHtvYmplY3R9IEFuIG9iamVjdCBjb250YWluaW5nIHRoZSBwbHVja2VkIGtleXNcbiAqL1xuY29uc3QgcGx1Y2sgPSAobywgLi4ua2V5cykgPT4ge1xuICByZXR1cm4ga2V5cy5yZWR1Y2UoXG4gICAgKHJlc3VsdCwga2V5KSA9PiBPYmplY3QuYXNzaWduKHJlc3VsdCwgeyBba2V5XTogb1trZXldIH0pLFxuICAgIHt9XG4gICk7XG59O1xuXG4vKipcbiAqIEV4Y2x1ZGUga2V5cyBmcm9tIGFuIG9iamVjdFxuICogQHBhcmFtIHtvYmplY3R9IG9cbiAqIEBwYXJhbSB7Li4uc3RyaW5nfSBrZXlzIFRoZSBrZXlzIHRvIGV4Y2x1ZGUgZnJvbSB0aGUgb2JqZWN0XG4gKiBAcmV0dXJuIHtvYmplY3R9IEFuIG9iamVjdCBjb250YWluaW5nIGFsbCBrZXlzIGV4Y2VwdCBleGNsdWRlZCBrZXlzXG4gKi9cbmNvbnN0IGV4Y2x1ZGUgPSAobywgLi4ua2V5cykgPT4ge1xuICByZXR1cm4gT2JqZWN0LmZyb21FbnRyaWVzKFxuICAgIE9iamVjdC5lbnRyaWVzKG8pLmZpbHRlcigoW2tleV0pID0+ICFrZXlzLmluY2x1ZGVzKGtleSkpXG4gICk7XG59O1xuXG5pZiAodHlwZW9mIG1vZHVsZSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgZmxvYXRFcXVhbHMsXG4gICAgY2xhbXAsXG4gICAgZnJhYyxcbiAgICByb3VuZCxcbiAgICBsZXJwLFxuICAgIHVubGVycCxcbiAgICBibGVycCxcbiAgICByZW1hcCxcbiAgICBzbW9vdGhzdGVwLFxuICAgIHJhZGlhbnMsXG4gICAgZGVncmVlcyxcbiAgICByYW5kb21CZXR3ZWVuLFxuICAgIHJhbmRvbUludEJldHdlZW4sXG4gICAgY2x0UmFuZG9tLFxuICAgIGNsdFJhbmRvbUludCxcbiAgICB3ZWlnaHRlZFJhbmRvbSxcbiAgICBsZXJwQXJyYXksXG4gICAgZG90LFxuICAgIGZhY3RvcmlhbCxcbiAgICBucHIsXG4gICAgbmNyLFxuICAgIGNvbWJpbmF0aW9ucyxcbiAgICBjYXJ0ZXNpYW4sXG4gICAgdGltZXMsXG4gICAgcmFuZ2UsXG4gICAgemlwLFxuICAgIGF0LFxuICAgIHBlZWssXG4gICAgY2h1bmssXG4gICAgc2h1ZmZsZSxcbiAgICBmbGF0LFxuICAgIHVuZmxhdCxcbiAgICBzcGxpdCxcbiAgICBwbHVjayxcbiAgICBleGNsdWRlLFxuICB9O1xufVxuIiwiY29uc3QgeyB0aW1lcywgY2h1bmssIGRvdCB9ID0gcmVxdWlyZSgnQGJhc2VtZW50dW5pdmVyc2UvdXRpbHMnKTtcblxuLyoqXG4gKiBAb3ZlcnZpZXcgQSBzbWFsbCB2ZWN0b3IgYW5kIG1hdHJpeCBsaWJyYXJ5XG4gKiBAYXV0aG9yIEdvcmRvbiBMYXJyaWdhblxuICovXG5cbi8qKlxuICogQSAyZCB2ZWN0b3JcbiAqIEB0eXBlZGVmIHtPYmplY3R9IHZlY1xuICogQHByb3BlcnR5IHtudW1iZXJ9IHggVGhlIHggY29tcG9uZW50IG9mIHRoZSB2ZWN0b3JcbiAqIEBwcm9wZXJ0eSB7bnVtYmVyfSB5IFRoZSB5IGNvbXBvbmVudCBvZiB0aGUgdmVjdG9yXG4gKi9cblxuLyoqXG4gKiBDcmVhdGUgYSBuZXcgdmVjdG9yXG4gKiBAcGFyYW0ge251bWJlcnx2ZWN9IFt4XSBUaGUgeCBjb21wb25lbnQgb2YgdGhlIHZlY3Rvciwgb3IgYSB2ZWN0b3IgdG8gY29weVxuICogQHBhcmFtIHtudW1iZXJ9IFt5XSBUaGUgeSBjb21wb25lbnQgb2YgdGhlIHZlY3RvclxuICogQHJldHVybiB7dmVjfSBBIG5ldyB2ZWN0b3JcbiAqIEBleGFtcGxlIDxjYXB0aW9uPlZhcmlvdXMgd2F5cyB0byBpbml0aWFsaXNlIGEgdmVjdG9yPC9jYXB0aW9uPlxuICogbGV0IGEgPSB2ZWMoMywgMik7ICAvLyAoMywgMilcbiAqIGxldCBiID0gdmVjKDQpOyAgICAgLy8gKDQsIDQpXG4gKiBsZXQgYyA9IHZlYyhhKTsgICAgIC8vICgzLCAyKVxuICogbGV0IGQgPSB2ZWMoKTsgICAgICAvLyAoMCwgMClcbiAqL1xuY29uc3QgdmVjID0gKHgsIHkpID0+ICgheCAmJiAheSA/XG4gIHsgeDogMCwgeTogMCB9IDogKHR5cGVvZiB4ID09PSAnb2JqZWN0JyA/XG4gICAgeyB4OiB4LnggfHwgMCwgeTogeC55IHx8IDAgfSA6ICh5ID09PSBudWxsIHx8IHkgPT09IHVuZGVmaW5lZCA/XG4gICAgICB7IHg6IHgsIHk6IHggfSA6IHsgeDogeCwgeTogeSB9KVxuICApXG4pO1xuXG4vKipcbiAqIEdldCB0aGUgY29tcG9uZW50cyBvZiBhIHZlY3RvciBhcyBhbiBhcnJheVxuICogQHBhcmFtIHt2ZWN9IGEgVGhlIHZlY3RvciB0byBnZXQgY29tcG9uZW50cyBmcm9tXG4gKiBAcmV0dXJuIHtBcnJheTxudW1iZXI+fSBUaGUgdmVjdG9yIGNvbXBvbmVudHMgYXMgYW4gYXJyYXlcbiAqL1xudmVjLmNvbXBvbmVudHMgPSBhID0+IFthLngsIGEueV07XG5cbi8qKlxuICogUmV0dXJuIGEgdW5pdCB2ZWN0b3IgKDEsIDApXG4gKiBAcmV0dXJuIHt2ZWN9IEEgdW5pdCB2ZWN0b3IgKDEsIDApXG4gKi9cbnZlYy51eCA9ICgpID0+IHZlYygxLCAwKTtcblxuLyoqXG4gKiBSZXR1cm4gYSB1bml0IHZlY3RvciAoMCwgMSlcbiAqIEByZXR1cm4ge3ZlY30gQSB1bml0IHZlY3RvciAoMCwgMSlcbiAqL1xudmVjLnV5ID0gKCkgPT4gdmVjKDAsIDEpO1xuXG4vKipcbiAqIEFkZCB2ZWN0b3JzXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHt2ZWN9IGIgVmVjdG9yIGJcbiAqIEByZXR1cm4ge3ZlY30gYSArIGJcbiAqL1xudmVjLmFkZCA9IChhLCBiKSA9PiAoeyB4OiBhLnggKyBiLngsIHk6IGEueSArIGIueSB9KTtcblxuLyoqXG4gKiBTY2FsZSBhIHZlY3RvclxuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBiIFNjYWxhciBiXG4gKiBAcmV0dXJuIHt2ZWN9IGEgKiBiXG4gKi9cbnZlYy5tdWwgPSAoYSwgYikgPT4gKHsgeDogYS54ICogYiwgeTogYS55ICogYiB9KTtcblxuLyoqXG4gKiBTdWJ0cmFjdCB2ZWN0b3JzXG4gKiBAcGFyYW0ge3ZlY30gYSBWZWN0b3IgYVxuICogQHBhcmFtIHt2ZWN9IGIgVmVjdG9yIGJcbiAqIEByZXR1cm4ge3ZlY30gYSAtIGJcbiAqL1xudmVjLnN1YiA9IChhLCBiKSA9PiAoeyB4OiBhLnggLSBiLngsIHk6IGEueSAtIGIueSB9KTtcblxuLyoqXG4gKiBHZXQgdGhlIGxlbmd0aCBvZiBhIHZlY3RvclxuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEByZXR1cm4ge251bWJlcn0gfGF8XG4gKi9cbnZlYy5sZW4gPSBhID0+IE1hdGguc3FydChhLnggKiBhLnggKyBhLnkgKiBhLnkpO1xuXG4vKipcbiAqIEdldCB0aGUgbGVuZ3RoIG9mIGEgdmVjdG9yIHVzaW5nIHRheGljYWIgZ2VvbWV0cnlcbiAqIEBwYXJhbSB7dmVjfSBhIFZlY3RvciBhXG4gKiBAcmV0dXJuIHtudW1iZXJ9IHxhfFxuICovXG52ZWMubWFuaGF0dGFuID0gYSA9PiBNYXRoLmFicyhhLngpICsgTWF0aC5hYnMoYS55KTtcblxuLyoqXG4gKiBOb3JtYWxpc2UgYSB2ZWN0b3JcbiAqIEBwYXJhbSB7dmVjfSBhIFRoZSB2ZWN0b3IgdG8gbm9ybWFsaXNlXG4gKiBAcmV0dXJuIHt2ZWN9IF5hXG4gKi9cbnZlYy5ub3IgPSBhID0+IHtcbiAgbGV0IGxlbiA9IHZlYy5sZW4oYSk7XG4gIHJldHVybiBsZW4gPyB7IHg6IGEueCAvIGxlbiwgeTogYS55IC8gbGVuIH0gOiB2ZWMoKTtcbn07XG5cbi8qKlxuICogR2V0IGEgZG90IHByb2R1Y3Qgb2YgdmVjdG9yc1xuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7dmVjfSBiIFZlY3RvciBiXG4gKiBAcmV0dXJuIHtudW1iZXJ9IGEg4oiZIGJcbiAqL1xudmVjLmRvdCA9IChhLCBiKSA9PiBhLnggKiBiLnggKyBhLnkgKiBiLnk7XG5cbi8qKlxuICogUm90YXRlIGEgdmVjdG9yIGJ5IHIgcmFkaWFuc1xuICogQHBhcmFtIHt2ZWN9IGEgVGhlIHZlY3RvciB0byByb3RhdGVcbiAqIEBwYXJhbSB7bnVtYmVyfSByIFRoZSBhbmdsZSB0byByb3RhdGUgYnksIG1lYXN1cmVkIGluIHJhZGlhbnNcbiAqIEByZXR1cm4ge3ZlY30gQSByb3RhdGVkIHZlY3RvclxuICovXG52ZWMucm90ID0gKGEsIHIpID0+IHtcbiAgbGV0IHMgPSBNYXRoLnNpbihyKSxcbiAgICBjID0gTWF0aC5jb3Mocik7XG4gIHJldHVybiB7IHg6IGMgKiBhLnggLSBzICogYS55LCB5OiBzICogYS54ICsgYyAqIGEueSB9O1xufVxuXG4vKipcbiAqIENoZWNrIGlmIHR3byB2ZWN0b3JzIGFyZSBlcXVhbFxuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7dmVjfSBiIFZlY3RvciBiXG4gKiBAcmV0dXJuIHtib29sZWFufSBUcnVlIGlmIHZlY3RvcnMgYSBhbmQgYiBhcmUgZXF1YWwsIGZhbHNlIG90aGVyd2lzZVxuICovXG52ZWMuZXEgPSAoYSwgYikgPT4gYS54ID09PSBiLnggJiYgYS55ID09PSBiLnk7XG5cbi8qKlxuICogR2V0IHRoZSBhbmdsZSBvZiBhIHZlY3RvclxuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEByZXR1cm4ge251bWJlcn0gVGhlIGFuZ2xlIG9mIHZlY3RvciBhIGluIHJhZGlhbnNcbiAqL1xudmVjLnJhZCA9IGEgPT4gTWF0aC5hdGFuMihhLnksIGEueCk7XG5cbi8qKlxuICogQ29weSBhIHZlY3RvclxuICogQHBhcmFtIHt2ZWN9IGEgVGhlIHZlY3RvciB0byBjb3B5XG4gKiBAcmV0dXJuIHt2ZWN9IEEgY29weSBvZiB2ZWN0b3IgYVxuICovXG52ZWMuY3B5ID0gYSA9PiB2ZWMoYSk7XG5cbi8qKlxuICogQSBmdW5jdGlvbiB0byBjYWxsIG9uIGVhY2ggY29tcG9uZW50IG9mIGEgdmVjdG9yXG4gKiBAY2FsbGJhY2sgdmVjdG9yTWFwQ2FsbGJhY2tcbiAqIEBwYXJhbSB7bnVtYmVyfSB2YWx1ZSBUaGUgY29tcG9uZW50IHZhbHVlXG4gKiBAcGFyYW0geyd4JyB8ICd5J30gbGFiZWwgVGhlIGNvbXBvbmVudCBsYWJlbCAoeCBvciB5KVxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgbWFwcGVkIGNvbXBvbmVudFxuICovXG5cbi8qKlxuICogQ2FsbCBhIGZ1bmN0aW9uIG9uIGVhY2ggY29tcG9uZW50IG9mIGEgdmVjdG9yIGFuZCBidWlsZCBhIG5ldyB2ZWN0b3IgZnJvbSB0aGUgcmVzdWx0c1xuICogQHBhcmFtIHt2ZWN9IGEgVmVjdG9yIGFcbiAqIEBwYXJhbSB7dmVjdG9yTWFwQ2FsbGJhY2t9IGYgVGhlIGZ1bmN0aW9uIHRvIGNhbGwgb24gZWFjaCBjb21wb25lbnQgb2YgdGhlIHZlY3RvclxuICogQHJldHVybiB7dmVjfSBWZWN0b3IgYSBtYXBwZWQgdGhyb3VnaCBmXG4gKi9cbnZlYy5tYXAgPSAoYSwgZikgPT4gKHsgeDogZihhLngsICd4JyksIHk6IGYoYS55LCAneScpIH0pO1xuXG4vKipcbiAqIENvbnZlcnQgYSB2ZWN0b3IgaW50byBhIHN0cmluZ1xuICogQHBhcmFtIHt2ZWN9IGEgVGhlIHZlY3RvciB0byBjb252ZXJ0XG4gKiBAcGFyYW0ge3N0cmluZ30gW3M9JywgJ10gVGhlIHNlcGFyYXRvciBzdHJpbmdcbiAqIEByZXR1cm4ge3N0cmluZ30gQSBzdHJpbmcgcmVwcmVzZW50YXRpb24gb2YgdGhlIHZlY3RvclxuICovXG52ZWMuc3RyID0gKGEsIHMgPSAnLCAnKSA9PiBgJHthLnh9JHtzfSR7YS55fWA7XG5cbi8qKlxuICogQSBtYXRyaXhcbiAqIEB0eXBlZGVmIHtPYmplY3R9IG1hdFxuICogQHByb3BlcnR5IHtudW1iZXJ9IG0gVGhlIG51bWJlciBvZiByb3dzIGluIHRoZSBtYXRyaXhcbiAqIEBwcm9wZXJ0eSB7bnVtYmVyfSBuIFRoZSBudW1iZXIgb2YgY29sdW1ucyBpbiB0aGUgbWF0cml4XG4gKiBAcHJvcGVydHkge0FycmF5PG51bWJlcj59IGVudHJpZXMgVGhlIG1hdHJpeCB2YWx1ZXNcbiAqL1xuXG4vKipcbiAqIENyZWF0ZSBhIG5ldyBtYXRyaXhcbiAqIEBwYXJhbSB7bnVtYmVyfSBbbT00XSBUaGUgbnVtYmVyIG9mIHJvd3NcbiAqIEBwYXJhbSB7bnVtYmVyfSBbbj00XSBUaGUgbnVtYmVyIG9mIGNvbHVtbnNcbiAqIEBwYXJhbSB7QXJyYXk8bnVtYmVyPn0gW2VudHJpZXM9W11dIE1hdHJpeCB2YWx1ZXMgaW4gcmVhZGluZyBvcmRlclxuICogQHJldHVybiB7bWF0fSBBIG5ldyBtYXRyaXhcbiAqL1xuY29uc3QgbWF0ID0gKG0gPSA0LCBuID0gNCwgZW50cmllcyA9IFtdKSA9PiAoe1xuICBtLCBuLFxuICBlbnRyaWVzOiBlbnRyaWVzLmNvbmNhdChBcnJheShtICogbikuZmlsbCgwKSkuc2xpY2UoMCwgbSAqIG4pXG59KTtcblxuLyoqXG4gKiBHZXQgYW4gaWRlbnRpdHkgbWF0cml4IG9mIHNpemUgblxuICogQHBhcmFtIHtudW1iZXJ9IG4gVGhlIHNpemUgb2YgdGhlIG1hdHJpeFxuICogQHJldHVybiB7bWF0fSBBbiBpZGVudGl0eSBtYXRyaXhcbiAqL1xubWF0LmlkZW50aXR5ID0gbiA9PiBtYXQobiwgbiwgQXJyYXkobiAqIG4pLmZpbGwoMCkubWFwKCh2LCBpKSA9PiArKE1hdGguZmxvb3IoaSAvIG4pID09PSBpICUgbikpKTtcblxuLyoqXG4gKiBHZXQgYW4gZW50cnkgZnJvbSBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBpIFRoZSByb3cgb2Zmc2V0XG4gKiBAcGFyYW0ge251bWJlcn0gaiBUaGUgY29sdW1uIG9mZnNldFxuICogQHJldHVybiB7bnVtYmVyfSBUaGUgdmFsdWUgYXQgcG9zaXRpb24gKGksIGopIGluIG1hdHJpeCBhXG4gKi9cbm1hdC5nZXQgPSAoYSwgaSwgaikgPT4gYS5lbnRyaWVzWyhqIC0gMSkgKyAoaSAtIDEpICogYS5uXTtcblxuLyoqXG4gKiBTZXQgYW4gZW50cnkgb2YgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge251bWJlcn0gaSBUaGUgcm93IG9mZnNldFxuICogQHBhcmFtIHtudW1iZXJ9IGogVGhlIGNvbHVtbiBvZmZzZXRcbiAqIEBwYXJhbSB7bnVtYmVyfSB2IFRoZSB2YWx1ZSB0byBzZXQgaW4gbWF0cml4IGFcbiAqL1xubWF0LnNldCA9IChhLCBpLCBqLCB2KSA9PiB7IGEuZW50cmllc1soaiAtIDEpICsgKGkgLSAxKSAqIGEubl0gPSB2OyB9O1xuXG4vKipcbiAqIEdldCBhIHJvdyBmcm9tIGEgbWF0cml4IGFzIGFuIGFycmF5XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IG0gVGhlIHJvdyBvZmZzZXRcbiAqIEByZXR1cm4ge0FycmF5PG51bWJlcj59IFJvdyBtIGZyb20gbWF0cml4IGFcbiAqL1xubWF0LnJvdyA9IChhLCBtKSA9PiB7XG4gIGNvbnN0IHMgPSAobSAtIDEpICogYS5uO1xuICByZXR1cm4gYS5lbnRyaWVzLnNsaWNlKHMsIHMgKyBhLm4pO1xufTtcblxuLyoqXG4gKiBHZXQgYSBjb2x1bW4gZnJvbSBhIG1hdHJpeCBhcyBhbiBhcnJheVxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBuIFRoZSBjb2x1bW4gb2Zmc2V0XG4gKiBAcmV0dXJuIHtBcnJheTxudW1iZXI+fSBDb2x1bW4gbiBmcm9tIG1hdHJpeCBhXG4gKi9cbm1hdC5jb2wgPSAoYSwgbikgPT4gdGltZXMoaSA9PiBtYXQuZ2V0KGEsIChpICsgMSksIG4pLCBhLm0pO1xuXG4vKipcbiAqIEFkZCBtYXRyaWNlc1xuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bWF0fSBiIE1hdHJpeCBiXG4gKiBAcmV0dXJuIHttYXR9IGEgKyBiXG4gKi9cbm1hdC5hZGQgPSAoYSwgYikgPT4gYS5tID09PSBiLm0gJiYgYS5uID09PSBiLm4gJiYgbWF0Lm1hcChhLCAodiwgaSkgPT4gdiArIGIuZW50cmllc1tpXSk7XG5cbi8qKlxuICogU3VidHJhY3QgbWF0cmljZXNcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge21hdH0gYiBNYXRyaXggYlxuICogQHJldHVybiB7bWF0fSBhIC0gYlxuICovXG5tYXQuc3ViID0gKGEsIGIpID0+IGEubSA9PT0gYi5tICYmIGEubiA9PT0gYi5uICYmIG1hdC5tYXAoYSwgKHYsIGkpID0+IHYgLSBiLmVudHJpZXNbaV0pO1xuXG4vKipcbiAqIE11bHRpcGx5IG1hdHJpY2VzXG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHttYXR9IGIgTWF0cml4IGJcbiAqIEByZXR1cm4ge21hdHxib29sZWFufSBhYiBvciBmYWxzZSBpZiB0aGUgbWF0cmljZXMgY2Fubm90IGJlIG11bHRpcGxpZWRcbiAqL1xubWF0Lm11bCA9IChhLCBiKSA9PiB7XG4gIGlmIChhLm4gIT09IGIubSkgeyByZXR1cm4gZmFsc2U7IH1cbiAgY29uc3QgcmVzdWx0ID0gbWF0KGEubSwgYi5uKTtcbiAgZm9yIChsZXQgaSA9IDE7IGkgPD0gYS5tOyBpKyspIHtcbiAgICBmb3IgKGxldCBqID0gMTsgaiA8PSBiLm47IGorKykge1xuICAgICAgbWF0LnNldChyZXN1bHQsIGksIGosIGRvdChtYXQucm93KGEsIGkpLCBtYXQuY29sKGIsIGopKSk7XG4gICAgfVxuICB9XG4gIHJldHVybiByZXN1bHQ7XG59O1xuXG4vKipcbiAqIFNjYWxlIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IGIgU2NhbGFyIGJcbiAqIEByZXR1cm4ge21hdH0gYSAqIGJcbiAqL1xubWF0LnNjYWxlID0gKGEsIGIpID0+IG1hdC5tYXAoYSwgdiA9PiB2ICogYik7XG5cbi8qKlxuICogVHJhbnNwb3NlIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IHRvIHRyYW5zcG9zZVxuICogQHJldHVybiB7bWF0fSBBIHRyYW5zcG9zZWQgbWF0cml4XG4gKi9cbm1hdC50cmFucyA9IGEgPT4gbWF0KGEubiwgYS5tLCB0aW1lcyhpID0+IG1hdC5jb2woYSwgKGkgKyAxKSksIGEubikuZmxhdCgpKTtcblxuLyoqXG4gKiBHZXQgdGhlIG1pbm9yIG9mIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBNYXRyaXggYVxuICogQHBhcmFtIHtudW1iZXJ9IGkgVGhlIHJvdyBvZmZzZXRcbiAqIEBwYXJhbSB7bnVtYmVyfSBqIFRoZSBjb2x1bW4gb2Zmc2V0XG4gKiBAcmV0dXJuIHttYXR8Ym9vbGVhbn0gVGhlIChpLCBqKSBtaW5vciBvZiBtYXRyaXggYSBvciBmYWxzZSBpZiB0aGUgbWF0cml4IGlzIG5vdCBzcXVhcmVcbiAqL1xubWF0Lm1pbm9yID0gKGEsIGksIGopID0+IHtcbiAgaWYgKGEubSAhPT0gYS5uKSB7IHJldHVybiBmYWxzZTsgfVxuICBjb25zdCBlbnRyaWVzID0gW107XG4gIGZvciAobGV0IGlpID0gMTsgaWkgPD0gYS5tOyBpaSsrKSB7XG4gICAgaWYgKGlpID09PSBpKSB7IGNvbnRpbnVlOyB9XG4gICAgZm9yIChsZXQgamogPSAxOyBqaiA8PSBhLm47IGpqKyspIHtcbiAgICAgIGlmIChqaiA9PT0gaikgeyBjb250aW51ZTsgfVxuICAgICAgZW50cmllcy5wdXNoKG1hdC5nZXQoYSwgaWksIGpqKSk7XG4gICAgfVxuICB9XG4gIHJldHVybiBtYXQoYS5tIC0gMSwgYS5uIC0gMSwgZW50cmllcyk7XG59O1xuXG4vKipcbiAqIEdldCB0aGUgZGV0ZXJtaW5hbnQgb2YgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcmV0dXJuIHtudW1iZXJ8Ym9vbGVhbn0gfGF8IG9yIGZhbHNlIGlmIHRoZSBtYXRyaXggaXMgbm90IHNxdWFyZVxuICovXG5tYXQuZGV0ID0gYSA9PiB7XG4gIGlmIChhLm0gIT09IGEubikgeyByZXR1cm4gZmFsc2U7IH1cbiAgaWYgKGEubSA9PT0gMSkge1xuICAgIHJldHVybiBhLmVudHJpZXNbMF07XG4gIH1cbiAgaWYgKGEubSA9PT0gMikge1xuICAgIHJldHVybiBhLmVudHJpZXNbMF0gKiBhLmVudHJpZXNbM10gLSBhLmVudHJpZXNbMV0gKiBhLmVudHJpZXNbMl07XG4gIH1cbiAgbGV0IHRvdGFsID0gMCwgc2lnbiA9IDE7XG4gIGZvciAobGV0IGogPSAxOyBqIDw9IGEubjsgaisrKSB7XG4gICAgdG90YWwgKz0gc2lnbiAqIGEuZW50cmllc1tqIC0gMV0gKiBtYXQuZGV0KG1hdC5taW5vcihhLCAxLCBqKSk7XG4gICAgc2lnbiAqPSAtMTtcbiAgfVxuICByZXR1cm4gdG90YWw7XG59O1xuXG4vKipcbiAqIE5vcm1hbGlzZSBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCB0byBub3JtYWxpc2VcbiAqIEByZXR1cm4ge21hdHxib29sZWFufSBeYSBvciBmYWxzZSBpZiB0aGUgbWF0cml4IGlzIG5vdCBzcXVhcmVcbiAqL1xubWF0Lm5vciA9IGEgPT4ge1xuICBpZiAoYS5tICE9PSBhLm4pIHsgcmV0dXJuIGZhbHNlOyB9XG4gIGNvbnN0IGQgPSBtYXQuZGV0KGEpO1xuICByZXR1cm4gbWF0Lm1hcChhLCBpID0+IGkgKiBkKTtcbn07XG5cbi8qKlxuICogR2V0IHRoZSBhZGp1Z2F0ZSBvZiBhIG1hdHJpeFxuICogQHBhcmFtIHttYXR9IGEgVGhlIG1hdHJpeCBmcm9tIHdoaWNoIHRvIGdldCB0aGUgYWRqdWdhdGVcbiAqIEByZXR1cm4ge21hdH0gVGhlIGFkanVnYXRlIG9mIGFcbiAqL1xubWF0LmFkaiA9IGEgPT4ge1xuICBjb25zdCBtaW5vcnMgPSBtYXQoYS5tLCBhLm4pO1xuICBmb3IgKGxldCBpID0gMTsgaSA8PSBhLm07IGkrKykge1xuICAgIGZvciAobGV0IGogPSAxOyBqIDw9IGEubjsgaisrKSB7XG4gICAgICBtYXQuc2V0KG1pbm9ycywgaSwgaiwgbWF0LmRldChtYXQubWlub3IoYSwgaSwgaikpKTtcbiAgICB9XG4gIH1cbiAgY29uc3QgY29mYWN0b3JzID0gbWF0Lm1hcChtaW5vcnMsICh2LCBpKSA9PiB2ICogKGkgJSAyID8gLTEgOiAxKSk7XG4gIHJldHVybiBtYXQudHJhbnMoY29mYWN0b3JzKTtcbn07XG5cbi8qKlxuICogR2V0IHRoZSBpbnZlcnNlIG9mIGEgbWF0cml4XG4gKiBAcGFyYW0ge21hdH0gYSBUaGUgbWF0cml4IHRvIGludmVydFxuICogQHJldHVybiB7bWF0fGJvb2xlYW59IGFeLTEgb3IgZmFsc2UgaWYgdGhlIG1hdHJpeCBoYXMgbm8gaW52ZXJzZVxuICovXG5tYXQuaW52ID0gYSA9PiB7XG4gIGlmIChhLm0gIT09IGEubikgeyByZXR1cm4gZmFsc2U7IH1cbiAgY29uc3QgZCA9IG1hdC5kZXQoYSk7XG4gIGlmIChkID09PSAwKSB7IHJldHVybiBmYWxzZTsgfVxuICByZXR1cm4gbWF0LnNjYWxlKG1hdC5hZGooYSksIDEgLyBkKTtcbn07XG5cbi8qKlxuICogQ2hlY2sgaWYgdHdvIG1hdHJpY2VzIGFyZSBlcXVhbFxuICogQHBhcmFtIHttYXR9IGEgTWF0cml4IGFcbiAqIEBwYXJhbSB7bWF0fSBiIE1hdHJpeCBiXG4gKiBAcmV0dXJuIHtib29sZWFufSBUcnVlIGlmIG1hdHJpY2VzIGEgYW5kIGIgYXJlIGlkZW50aWNhbCwgZmFsc2Ugb3RoZXJ3aXNlXG4gKi9cbm1hdC5lcSA9IChhLCBiKSA9PiBhLm0gPT09IGIubSAmJiBhLm4gPT09IGIubiAmJiBtYXQuc3RyKGEpID09PSBtYXQuc3RyKGIpO1xuXG4vKipcbiAqIENvcHkgYSBtYXRyaXhcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggdG8gY29weVxuICogQHJldHVybiB7bWF0fSBBIGNvcHkgb2YgbWF0cml4IGFcbiAqL1xubWF0LmNweSA9IGEgPT4gbWF0KGEubSwgYS5uLCBbLi4uYS5lbnRyaWVzXSk7XG5cbi8qKlxuICogQSBmdW5jdGlvbiB0byBjYWxsIG9uIGVhY2ggZW50cnkgb2YgYSBtYXRyaXhcbiAqIEBjYWxsYmFjayBtYXRyaXhNYXBDYWxsYmFja1xuICogQHBhcmFtIHtudW1iZXJ9IHZhbHVlIFRoZSBlbnRyeSB2YWx1ZVxuICogQHBhcmFtIHtudW1iZXJ9IGluZGV4IFRoZSBlbnRyeSBpbmRleFxuICogQHBhcmFtIHtBcnJheTxudW1iZXI+fSBlbnRyaWVzIFRoZSBhcnJheSBvZiBtYXRyaXggZW50cmllc1xuICogQHJldHVybiB7bnVtYmVyfSBUaGUgbWFwcGVkIGVudHJ5XG4gKi9cblxuLyoqXG4gKiBDYWxsIGEgZnVuY3Rpb24gb24gZWFjaCBlbnRyeSBvZiBhIG1hdHJpeCBhbmQgYnVpbGQgYSBuZXcgbWF0cml4IGZyb20gdGhlIHJlc3VsdHNcbiAqIEBwYXJhbSB7bWF0fSBhIE1hdHJpeCBhXG4gKiBAcGFyYW0ge21hdHJpeE1hcENhbGxiYWNrfSBmIFRoZSBmdW5jdGlvbiB0byBjYWxsIG9uIGVhY2ggZW50cnkgb2YgdGhlIG1hdHJpeFxuICogQHJldHVybiB7bWF0fSBNYXRyaXggYSBtYXBwZWQgdGhyb3VnaCBmXG4gKi9cbm1hdC5tYXAgPSAoYSwgZikgPT4gbWF0KGEubSwgYS5uLCBhLmVudHJpZXMubWFwKGYpKTtcblxuLyoqXG4gKiBDb252ZXJ0IGEgbWF0cml4IGludG8gYSBzdHJpbmdcbiAqIEBwYXJhbSB7bWF0fSBhIFRoZSBtYXRyaXggdG8gY29udmVydFxuICogQHBhcmFtIHtzdHJpbmd9IFttcz0nLCAnXSBUaGUgc2VwYXJhdG9yIHN0cmluZyBmb3IgY29sdW1uc1xuICogQHBhcmFtIHtzdHJpbmd9IFtucz0nXFxuJ10gVGhlIHNlcGFyYXRvciBzdHJpbmcgZm9yIHJvd3NcbiAqIEByZXR1cm4ge3N0cmluZ30gQSBzdHJpbmcgcmVwcmVzZW50YXRpb24gb2YgdGhlIG1hdHJpeFxuICovXG5tYXQuc3RyID0gKGEsIG1zID0gJywgJywgbnMgPSAnXFxuJykgPT4gY2h1bmsoYS5lbnRyaWVzLCBhLm4pLm1hcChyID0+IHIuam9pbihtcykpLmpvaW4obnMpO1xuXG5pZiAodHlwZW9mIG1vZHVsZSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgbW9kdWxlLmV4cG9ydHMgPSB7IHZlYywgbWF0IH07XG59XG4iLCJcInVzZSBzdHJpY3RcIjtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcIl9fZXNNb2R1bGVcIiwgeyB2YWx1ZTogdHJ1ZSB9KTtcbmV4cG9ydHMuU3ByaXRlID0gZXhwb3J0cy5TcHJpdGVBbmltYXRpb25SZXBlYXRNb2RlID0gdm9pZCAwO1xuY29uc3QgdmVjXzEgPSByZXF1aXJlKFwiQGJhc2VtZW50dW5pdmVyc2UvdmVjXCIpO1xudmFyIFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGU7XG4oZnVuY3Rpb24gKFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGUpIHtcbiAgICAvKipcbiAgICAgKiBMb29wIHRoaXMgYW5pbWF0aW9uIGluZGVmaW5pdGVseVxuICAgICAqL1xuICAgIFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGVbU3ByaXRlQW5pbWF0aW9uUmVwZWF0TW9kZVtcIlJlcGVhdFwiXSA9IDBdID0gXCJSZXBlYXRcIjtcbiAgICAvKipcbiAgICAgKiBQbGF5IG9uY2UgYW5kIHRoZW4gc3RvcCBvbiB0aGUgbGFzdCBmcmFtZVxuICAgICAqL1xuICAgIFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGVbU3ByaXRlQW5pbWF0aW9uUmVwZWF0TW9kZVtcIlBsYXlPbmNlQW5kU3RvcFwiXSA9IDFdID0gXCJQbGF5T25jZUFuZFN0b3BcIjtcbiAgICAvKipcbiAgICAgKiBQbGF5IG9uY2UgYW5kIHRoZW4gcmVzZXQgYmFjayB0byB0aGUgZmlyc3QgZnJhbWVcbiAgICAgKi9cbiAgICBTcHJpdGVBbmltYXRpb25SZXBlYXRNb2RlW1Nwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGVbXCJQbGF5T25jZUFuZFJlc2V0XCJdID0gMl0gPSBcIlBsYXlPbmNlQW5kUmVzZXRcIjtcbn0pKFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGUgPSBleHBvcnRzLlNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGUgfHwgKGV4cG9ydHMuU3ByaXRlQW5pbWF0aW9uUmVwZWF0TW9kZSA9IHt9KSk7XG5jbGFzcyBTcHJpdGUge1xuICAgIGNvbnN0cnVjdG9yKG9wdGlvbnMpIHtcbiAgICAgICAgdmFyIF9hLCBfYjtcbiAgICAgICAgdGhpcy5wb3NpdGlvbiA9ICgwLCB2ZWNfMS52ZWMpKCk7XG4gICAgICAgIHRoaXMuc2l6ZSA9ICgwLCB2ZWNfMS52ZWMpKCk7XG4gICAgICAgIHRoaXMub3JpZ2luID0gKDAsIHZlY18xLnZlYykoKTtcbiAgICAgICAgdGhpcy5zY2FsZSA9IDE7XG4gICAgICAgIHRoaXMucm90YXRpb24gPSAwO1xuICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zID0gbnVsbDtcbiAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uU3RhdGUgPSBudWxsO1xuICAgICAgICB0aGlzLmN1cnJlbnRJbWFnZSA9IG51bGw7XG4gICAgICAgIHRoaXMuY3VycmVudEF0dGFjaG1lbnRQb2ludHMgPSBudWxsO1xuICAgICAgICBjb25zdCBhY3R1YWxPcHRpb25zID0gT2JqZWN0LmFzc2lnbih7fSwgU3ByaXRlLkRFRkFVTFRfT1BUSU9OUywgb3B0aW9ucyAhPT0gbnVsbCAmJiBvcHRpb25zICE9PSB2b2lkIDAgPyBvcHRpb25zIDoge30pO1xuICAgICAgICBpZiAoIWFjdHVhbE9wdGlvbnMuZGVidWcgfHwgYWN0dWFsT3B0aW9ucy5kZWJ1ZyA9PT0gdHJ1ZSkge1xuICAgICAgICAgICAgYWN0dWFsT3B0aW9ucy5kZWJ1ZyA9IHtcbiAgICAgICAgICAgICAgICBzaG93U3ByaXRlVHJhbnNmb3JtczogISFhY3R1YWxPcHRpb25zLmRlYnVnLFxuICAgICAgICAgICAgICAgIHNob3dTcHJpdGVCb3VuZGluZ0JveDogISFhY3R1YWxPcHRpb25zLmRlYnVnLFxuICAgICAgICAgICAgICAgIHNob3dBdHRhY2htZW50UG9pbnRzOiAhIWFjdHVhbE9wdGlvbnMuZGVidWcsXG4gICAgICAgICAgICB9O1xuICAgICAgICB9XG4gICAgICAgIHRoaXMub3B0aW9ucyA9IGFjdHVhbE9wdGlvbnM7XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMucG9zaXRpb24pIHtcbiAgICAgICAgICAgIHRoaXMucG9zaXRpb24gPSB2ZWNfMS52ZWMuY3B5KHRoaXMub3B0aW9ucy5wb3NpdGlvbik7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5zaXplKSB7XG4gICAgICAgICAgICB0aGlzLnNpemUgPSB2ZWNfMS52ZWMuY3B5KHRoaXMub3B0aW9ucy5zaXplKTtcbiAgICAgICAgfVxuICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgIC8vIERlZmF1bHQgdG8gdGhlIHNpemUgb2YgdGhlIGJhc2UgaW1hZ2UgaWYgb25lIGV4aXN0c1xuICAgICAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5pbWFnZSkge1xuICAgICAgICAgICAgICAgIHRoaXMuc2l6ZSA9ICgwLCB2ZWNfMS52ZWMpKHRoaXMub3B0aW9ucy5pbWFnZS53aWR0aCwgdGhpcy5vcHRpb25zLmltYWdlLmhlaWdodCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICAvLyBGYWxsIGJhY2sgdG8gdGhlIHNpemUgb2YgdGhlIGltYWdlIGluIHRoZSBmaXJzdCBmcmFtZSBvZiB0aGUgZmlyc3RcbiAgICAgICAgICAgICAgICAvLyBhdmFpbGFibGUgZGlyZWN0aW9uIG9mIHRoZSBkZWZhdWx0IGFuaW1hdGlvbiBpZiBvbmUgZXhpc3RzXG4gICAgICAgICAgICAgICAgY29uc3QgZGVmYXVsdEFuaW1hdGlvbkRpcmVjdGlvbnMgPSBPYmplY3QudmFsdWVzKHRoaXMub3B0aW9ucy5hbmltYXRpb25zW3RoaXMub3B0aW9ucy5kZWZhdWx0QW5pbWF0aW9uXSlbMF07XG4gICAgICAgICAgICAgICAgaWYgKGRlZmF1bHRBbmltYXRpb25EaXJlY3Rpb25zICYmXG4gICAgICAgICAgICAgICAgICAgICgoX2IgPSAoX2EgPSBkZWZhdWx0QW5pbWF0aW9uRGlyZWN0aW9ucy5pbWFnZXMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS5sZW5ndGgpICE9PSBudWxsICYmIF9iICE9PSB2b2lkIDAgPyBfYiA6IDApID4gMCkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLnNpemUgPSAoMCwgdmVjXzEudmVjKShkZWZhdWx0QW5pbWF0aW9uRGlyZWN0aW9ucy5pbWFnZXNbMF0ud2lkdGgsIGRlZmF1bHRBbmltYXRpb25EaXJlY3Rpb25zLmltYWdlc1swXS5oZWlnaHQpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIC8vIE90aGVyd2lzZSBsZWF2ZSB0aGUgc2l6ZSBhcyAoMCwgMClcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLm9yaWdpbikge1xuICAgICAgICAgICAgdGhpcy5vcmlnaW4gPSB2ZWNfMS52ZWMuY3B5KHRoaXMub3B0aW9ucy5vcmlnaW4pO1xuICAgICAgICB9XG4gICAgICAgIGVsc2Uge1xuICAgICAgICAgICAgLy8gRGVmYXVsdCB0byB0aGUgY2VudGVyIG9mIHRoZSBzcHJpdGUgYmFzZWQgb24gc2l6ZVxuICAgICAgICAgICAgdGhpcy5vcmlnaW4gPSB2ZWNfMS52ZWMubXVsKHRoaXMuc2l6ZSwgMC41KTtcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLnNjYWxlKSB7XG4gICAgICAgICAgICB0aGlzLnNjYWxlID0gdGhpcy5vcHRpb25zLnNjYWxlO1xuICAgICAgICB9XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMucm90YXRpb24pIHtcbiAgICAgICAgICAgIHRoaXMucm90YXRpb24gPSB0aGlzLm9wdGlvbnMucm90YXRpb247XG4gICAgICAgIH1cbiAgICAgICAgLy8gQ2hlY2sgYW5kIGluaXRpYWxpc2UgZGlyZWN0aW9uXG4gICAgICAgIHRoaXMuX2RpcmVjdGlvbiA9IHRoaXMub3B0aW9ucy5kZWZhdWx0RGlyZWN0aW9uO1xuICAgICAgICBpZiAodGhpcy5vcHRpb25zLmRpcmVjdGlvbnMubGVuZ3RoID09PSAwIHx8XG4gICAgICAgICAgICAhdGhpcy5vcHRpb25zLmRpcmVjdGlvbnMuaW5jbHVkZXModGhpcy5fZGlyZWN0aW9uKSkge1xuICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKGBJbnZhbGlkIGRpcmVjdGlvbiBcIiR7dGhpcy5fZGlyZWN0aW9ufVwiYCk7XG4gICAgICAgIH1cbiAgICAgICAgLy8gQ2hlY2sgYW5kIGluaXRpYWxpc2UgYW5pbWF0aW9uXG4gICAgICAgIHRoaXMuX2FuaW1hdGlvbiA9IHRoaXMub3B0aW9ucy5kZWZhdWx0QW5pbWF0aW9uO1xuICAgICAgICBjb25zdCBhbmltYXRpb25zID0gT2JqZWN0LmtleXModGhpcy5vcHRpb25zLmFuaW1hdGlvbnMpO1xuICAgICAgICBpZiAoYW5pbWF0aW9ucy5sZW5ndGggPT09IDAgfHxcbiAgICAgICAgICAgICFhbmltYXRpb25zLmluY2x1ZGVzKHRoaXMuX2FuaW1hdGlvbikpIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgSW52YWxpZCBhbmltYXRpb24gXCIke3RoaXMuX2FuaW1hdGlvbn1cImApO1xuICAgICAgICB9XG4gICAgICAgIC8vIE1ha2Ugc3VyZSBhdHRhY2htZW50IHBvaW50IGtleWZyYW1lcyBhcmUgZGVmaW5lZCBpbiBhc2NlbmRpbmdcbiAgICAgICAgLy8gZnJhbWUgb3JkZXIgaW4gYWxsIGFuaW1hdGlvbnNcbiAgICAgICAgZm9yIChjb25zdCBhbmltYXRpb24gb2YgT2JqZWN0LmtleXModGhpcy5vcHRpb25zLmFuaW1hdGlvbnMpKSB7XG4gICAgICAgICAgICBmb3IgKGNvbnN0IGRpcmVjdGlvbiBvZiBPYmplY3Qua2V5cyh0aGlzLm9wdGlvbnMuYW5pbWF0aW9uc1thbmltYXRpb25dKSkge1xuICAgICAgICAgICAgICAgIGlmICh0aGlzLm9wdGlvbnMuYW5pbWF0aW9uc1thbmltYXRpb25dW2RpcmVjdGlvbl0uYXR0YWNobWVudFBvaW50S2V5ZnJhbWVzKSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAoY29uc3QgYXR0YWNobWVudFBvaW50IG9mIE9iamVjdC5rZXlzKHRoaXNcbiAgICAgICAgICAgICAgICAgICAgICAgIC5vcHRpb25zXG4gICAgICAgICAgICAgICAgICAgICAgICAuYW5pbWF0aW9uc1thbmltYXRpb25dW2RpcmVjdGlvbl1cbiAgICAgICAgICAgICAgICAgICAgICAgIC5hdHRhY2htZW50UG9pbnRLZXlmcmFtZXMpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLm9wdGlvbnNcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAuYW5pbWF0aW9uc1thbmltYXRpb25dW2RpcmVjdGlvbl1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAuYXR0YWNobWVudFBvaW50S2V5ZnJhbWVzW2F0dGFjaG1lbnRQb2ludF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAuc29ydCgoYSwgYikgPT4gYS5mcmFtZSAtIGIuZnJhbWUpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuICAgIGdldCBkaXJlY3Rpb24oKSB7XG4gICAgICAgIHJldHVybiB0aGlzLl9kaXJlY3Rpb247XG4gICAgfVxuICAgIHNldCBkaXJlY3Rpb24odmFsdWUpIHtcbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kaXJlY3Rpb25zLmluY2x1ZGVzKHZhbHVlKSkge1xuICAgICAgICAgICAgdGhpcy5fZGlyZWN0aW9uID0gdmFsdWU7XG4gICAgICAgIH1cbiAgICB9XG4gICAgZ2V0IGFuaW1hdGlvbigpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMuX2FuaW1hdGlvbjtcbiAgICB9XG4gICAgc2V0IGFuaW1hdGlvbih2YWx1ZSkge1xuICAgICAgICBpZiAoT2JqZWN0LmtleXModGhpcy5vcHRpb25zLmFuaW1hdGlvbnMpLmluY2x1ZGVzKHZhbHVlKSkge1xuICAgICAgICAgICAgdGhpcy5fYW5pbWF0aW9uID0gdmFsdWU7XG4gICAgICAgIH1cbiAgICB9XG4gICAgcGxheUFuaW1hdGlvbigpIHtcbiAgICAgICAgaWYgKHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlKSB7XG4gICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5wbGF5aW5nID0gdHJ1ZTtcbiAgICAgICAgfVxuICAgIH1cbiAgICBwYXVzZUFuaW1hdGlvbigpIHtcbiAgICAgICAgaWYgKHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlKSB7XG4gICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5wbGF5aW5nID0gZmFsc2U7XG4gICAgICAgIH1cbiAgICB9XG4gICAgcmVzZXRBbmltYXRpb24oKSB7XG4gICAgICAgIGlmICh0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZSkge1xuICAgICAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uU3RhdGUuY3VycmVudEZyYW1lID0gMDtcbiAgICAgICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlLmN1cnJlbnRGcmFtZVRpbWUgPSAwO1xuICAgICAgICB9XG4gICAgfVxuICAgIGdldEF0dGFjaG1lbnRQb2ludChuYW1lKSB7XG4gICAgICAgIHZhciBfYSwgX2I7XG4gICAgICAgIHJldHVybiAoX2IgPSAoX2EgPSB0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzKSA9PT0gbnVsbCB8fCBfYSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2FbbmFtZV0pICE9PSBudWxsICYmIF9iICE9PSB2b2lkIDAgPyBfYiA6IG51bGw7XG4gICAgfVxuICAgIHVwZGF0ZShkdCkge1xuICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zID0gdGhpcy51cGRhdGVBbmltYXRpb25PcHRpb25zKCk7XG4gICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlID0gdGhpcy51cGRhdGVBbmltYXRpb25TdGF0ZShkdCk7XG4gICAgICAgIHRoaXMuY3VycmVudEltYWdlID0gdGhpcy51cGRhdGVJbWFnZSgpO1xuICAgICAgICB0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzID0gdGhpcy51cGRhdGVBdHRhY2htZW50UG9pbnRzKCk7XG4gICAgfVxuICAgIHVwZGF0ZUFuaW1hdGlvbk9wdGlvbnMoKSB7XG4gICAgICAgIGlmICghKHRoaXMuX2FuaW1hdGlvbiBpbiB0aGlzLm9wdGlvbnMuYW5pbWF0aW9ucykpIHtcbiAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihgSW52YWxpZCBhbmltYXRpb24gXCIke3RoaXMuX2FuaW1hdGlvbn1cImApO1xuICAgICAgICB9XG4gICAgICAgIGNvbnN0IGRpcmVjdGlvbnMgPSBPYmplY3Qua2V5cyh0aGlzLm9wdGlvbnMuYW5pbWF0aW9uc1t0aGlzLl9hbmltYXRpb25dKTtcbiAgICAgICAgaWYgKGRpcmVjdGlvbnMubGVuZ3RoID09PSAwKSB7XG4gICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoYE5vIGRpcmVjdGlvbnMgYXZhaWxhYmxlIGZvciBhbmltYXRpb24gXCIke3RoaXMuX2FuaW1hdGlvbn1cImApO1xuICAgICAgICB9XG4gICAgICAgIGlmICh0aGlzLl9kaXJlY3Rpb24gaW4gdGhpcy5vcHRpb25zLmFuaW1hdGlvbnNbdGhpcy5fYW5pbWF0aW9uXSkge1xuICAgICAgICAgICAgcmV0dXJuIHRoaXMub3B0aW9ucy5hbmltYXRpb25zW3RoaXMuX2FuaW1hdGlvbl1bdGhpcy5fZGlyZWN0aW9uXTtcbiAgICAgICAgfVxuICAgICAgICBpZiAoJyonIGluIHRoaXMub3B0aW9ucy5hbmltYXRpb25zW3RoaXMuX2FuaW1hdGlvbl0pIHtcbiAgICAgICAgICAgIHJldHVybiB0aGlzLm9wdGlvbnMuYW5pbWF0aW9uc1t0aGlzLl9hbmltYXRpb25dWycqJ107XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHRoaXMub3B0aW9ucy5hbmltYXRpb25zW3RoaXMuX2FuaW1hdGlvbl1bZGlyZWN0aW9uc1swXV07XG4gICAgfVxuICAgIHVwZGF0ZUFuaW1hdGlvblN0YXRlKGR0KSB7XG4gICAgICAgIGlmICghdGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucyB8fFxuICAgICAgICAgICAgIXRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlKSB7XG4gICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgIHBsYXlpbmc6IHRydWUsXG4gICAgICAgICAgICAgICAgY3VycmVudEZyYW1lOiAwLFxuICAgICAgICAgICAgICAgIGN1cnJlbnRGcmFtZVRpbWU6IDAsXG4gICAgICAgICAgICB9O1xuICAgICAgICB9XG4gICAgICAgIGlmICh0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5wbGF5aW5nKSB7XG4gICAgICAgICAgICBjb25zdCBmcmFtZVRpbWUgPSAxIC8gdGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucy5mcmFtZVJhdGU7XG4gICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5jdXJyZW50RnJhbWVUaW1lICs9IGR0O1xuICAgICAgICAgICAgaWYgKHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlLmN1cnJlbnRGcmFtZVRpbWUgPiBmcmFtZVRpbWUpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBmcmFtZUNvdW50ID0gdGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucy5mcmFtZUNvdW50O1xuICAgICAgICAgICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlLmN1cnJlbnRGcmFtZSsrO1xuICAgICAgICAgICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlLmN1cnJlbnRGcmFtZVRpbWUgPSAwO1xuICAgICAgICAgICAgICAgIGlmICh0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5jdXJyZW50RnJhbWUgPiBmcmFtZUNvdW50KSB7XG4gICAgICAgICAgICAgICAgICAgIHN3aXRjaCAodGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucy5tb2RlKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBjYXNlIFNwcml0ZUFuaW1hdGlvblJlcGVhdE1vZGUuUGxheU9uY2VBbmRSZXNldDpcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5wbGF5aW5nID0gZmFsc2U7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uU3RhdGUuY3VycmVudEZyYW1lID0gMDtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIGNhc2UgU3ByaXRlQW5pbWF0aW9uUmVwZWF0TW9kZS5QbGF5T25jZUFuZFN0b3A6XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uU3RhdGUucGxheWluZyA9IGZhbHNlO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlLmN1cnJlbnRGcmFtZSA9IGZyYW1lQ291bnQgLSAxO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICAgICAgY2FzZSBTcHJpdGVBbmltYXRpb25SZXBlYXRNb2RlLlJlcGVhdDpcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5jdXJyZW50RnJhbWUgPSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZTtcbiAgICB9XG4gICAgdXBkYXRlSW1hZ2UoKSB7XG4gICAgICAgIHZhciBfYSwgX2IsIF9jO1xuICAgICAgICBpZiAoIXRoaXMuY3VycmVudEFuaW1hdGlvbk9wdGlvbnMgfHxcbiAgICAgICAgICAgICF0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZSkge1xuICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKCF0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zLmltYWdlcyB8fFxuICAgICAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucy5pbWFnZXMubGVuZ3RoID09PSAwKSB7XG4gICAgICAgICAgICByZXR1cm4gKF9hID0gdGhpcy5vcHRpb25zLmltYWdlKSAhPT0gbnVsbCAmJiBfYSAhPT0gdm9pZCAwID8gX2EgOiBudWxsO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiAoX2MgPSAoX2IgPSB0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zLmltYWdlc1t0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5jdXJyZW50RnJhbWVdKSAhPT0gbnVsbCAmJiBfYiAhPT0gdm9pZCAwID8gX2IgOiB0aGlzLm9wdGlvbnMuaW1hZ2UpICE9PSBudWxsICYmIF9jICE9PSB2b2lkIDAgPyBfYyA6IG51bGw7XG4gICAgfVxuICAgIHVwZGF0ZUF0dGFjaG1lbnRQb2ludHMoKSB7XG4gICAgICAgIGlmICghdGhpcy5vcHRpb25zLmF0dGFjaG1lbnRQb2ludHMgfHxcbiAgICAgICAgICAgIHRoaXMub3B0aW9ucy5hdHRhY2htZW50UG9pbnRzLmxlbmd0aCA9PT0gMCkge1xuICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKCF0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzKSB7XG4gICAgICAgICAgICB0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzID0gT2JqZWN0LmZyb21FbnRyaWVzKHRoaXMub3B0aW9ucy5hdHRhY2htZW50UG9pbnRzLm1hcChhdHRhY2htZW50UG9pbnQgPT4gW1xuICAgICAgICAgICAgICAgIGF0dGFjaG1lbnRQb2ludC5uYW1lLFxuICAgICAgICAgICAgICAgIGF0dGFjaG1lbnRQb2ludC5vZmZzZXRcbiAgICAgICAgICAgIF0pKTtcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucyAmJlxuICAgICAgICAgICAgdGhpcy5jdXJyZW50QW5pbWF0aW9uT3B0aW9ucy5hdHRhY2htZW50UG9pbnRLZXlmcmFtZXMgJiZcbiAgICAgICAgICAgIHRoaXMuY3VycmVudEFuaW1hdGlvblN0YXRlKSB7XG4gICAgICAgICAgICBmb3IgKGNvbnN0IG5hbWUgb2YgT2JqZWN0LmtleXModGhpcy5jdXJyZW50QXR0YWNobWVudFBvaW50cykpIHtcbiAgICAgICAgICAgICAgICBpZiAobmFtZSBpbiB0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zLmF0dGFjaG1lbnRQb2ludEtleWZyYW1lcyAmJlxuICAgICAgICAgICAgICAgICAgICB0aGlzLmN1cnJlbnRBbmltYXRpb25PcHRpb25zLmF0dGFjaG1lbnRQb2ludEtleWZyYW1lc1tuYW1lXS5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnN0IHByZXZpb3VzS2V5ZnJhbWUgPSB0aGlzLmZpbmRQcmV2aW91c0tleWZyYW1lKHRoaXMuY3VycmVudEFuaW1hdGlvbk9wdGlvbnMuYXR0YWNobWVudFBvaW50S2V5ZnJhbWVzW25hbWVdLCB0aGlzLmN1cnJlbnRBbmltYXRpb25TdGF0ZS5jdXJyZW50RnJhbWUpO1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzW25hbWVdID0gcHJldmlvdXNLZXlmcmFtZS5vZmZzZXQ7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHJldHVybiB0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzO1xuICAgIH1cbiAgICBmaW5kUHJldmlvdXNLZXlmcmFtZShrZXlmcmFtZXMsIGN1cnJlbnRGcmFtZSkge1xuICAgICAgICBjb25zdCBmb3VuZCA9IFsuLi5rZXlmcmFtZXNdLnJldmVyc2UoKS5maW5kKGtleWZyYW1lID0+IGtleWZyYW1lLmZyYW1lIDw9IGN1cnJlbnRGcmFtZSk7XG4gICAgICAgIGlmICghZm91bmQpIHtcbiAgICAgICAgICAgIHJldHVybiBrZXlmcmFtZXNba2V5ZnJhbWVzLmxlbmd0aCAtIDFdO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBmb3VuZDtcbiAgICB9XG4gICAgZHJhdyhjb250ZXh0KSB7XG4gICAgICAgIGNvbnRleHQuc2F2ZSgpO1xuICAgICAgICBjb250ZXh0LnRyYW5zbGF0ZSh0aGlzLnBvc2l0aW9uLngsIHRoaXMucG9zaXRpb24ueSk7XG4gICAgICAgIGNvbnRleHQuc2NhbGUodGhpcy5zY2FsZSwgdGhpcy5zY2FsZSk7XG4gICAgICAgIGNvbnRleHQucm90YXRlKHRoaXMucm90YXRpb24pO1xuICAgICAgICBpZiAodGhpcy5vcHRpb25zLnByZVJlbmRlcikge1xuICAgICAgICAgICAgdGhpcy5vcHRpb25zLnByZVJlbmRlcihjb250ZXh0LCB0aGlzKTtcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5jdXJyZW50SW1hZ2UpIHtcbiAgICAgICAgICAgIGNvbnRleHQuZHJhd0ltYWdlKHRoaXMuY3VycmVudEltYWdlLCAtdGhpcy5vcmlnaW4ueCwgLXRoaXMub3JpZ2luLnksIHRoaXMuY3VycmVudEltYWdlLndpZHRoLCB0aGlzLmN1cnJlbnRJbWFnZS5oZWlnaHQpO1xuICAgICAgICB9XG4gICAgICAgIGlmICh0aGlzLm9wdGlvbnMucG9zdFJlbmRlcikge1xuICAgICAgICAgICAgdGhpcy5vcHRpb25zLnBvc3RSZW5kZXIoY29udGV4dCwgdGhpcyk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93U3ByaXRlQm91bmRpbmdCb3gpIHtcbiAgICAgICAgICAgIGNvbnRleHQuc3Ryb2tlU3R5bGUgPSBTcHJpdGUuREVCVUdfQk9VTkRJTkdfQk9YX0NPTE9VUjtcbiAgICAgICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gU3ByaXRlLkRFQlVHX0JPVU5ESU5HX0JPWF9MSU5FX1dJRFRIO1xuICAgICAgICAgICAgY29udGV4dC5zdHJva2VSZWN0KC10aGlzLm9yaWdpbi54LCAtdGhpcy5vcmlnaW4ueSwgdGhpcy5zaXplLngsIHRoaXMuc2l6ZS55KTtcbiAgICAgICAgfVxuICAgICAgICBpZiAodGhpcy5vcHRpb25zLmRlYnVnLnNob3dTcHJpdGVUcmFuc2Zvcm1zKSB7XG4gICAgICAgICAgICB0aGlzLmRyYXdUcmFuc2Zvcm1zTWFya2VyKGNvbnRleHQsICgwLCB2ZWNfMS52ZWMpKCksIFNwcml0ZS5ERUJVR19UUkFOU0ZPUk1TX0NPTE9VUl9YLCBTcHJpdGUuREVCVUdfVFJBTlNGT1JNU19DT0xPVVJfWSwgU3ByaXRlLkRFQlVHX1RSQU5TRk9STVNfTElORV9XSURUSCwgU3ByaXRlLkRFQlVHX1RSQU5TRk9STVNfU0laRSk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHRoaXMub3B0aW9ucy5kZWJ1Zy5zaG93QXR0YWNobWVudFBvaW50cyAmJlxuICAgICAgICAgICAgdGhpcy5jdXJyZW50QXR0YWNobWVudFBvaW50cykge1xuICAgICAgICAgICAgZm9yIChjb25zdCBhdHRhY2htZW50UG9pbnQgb2YgT2JqZWN0LnZhbHVlcyh0aGlzLmN1cnJlbnRBdHRhY2htZW50UG9pbnRzKSkge1xuICAgICAgICAgICAgICAgIHRoaXMuZHJhd0Nyb3NzKGNvbnRleHQsIGF0dGFjaG1lbnRQb2ludCwgU3ByaXRlLkRFQlVHX0FUVEFDSE1FTlRfUE9JTlRfQ09MT1VSLCBTcHJpdGUuREVCVUdfQVRUQUNITUVOVF9QT0lOVF9MSU5FX1dJRFRILCBTcHJpdGUuREVCVUdfQVRUQUNITUVOVF9QT0lOVF9TSVpFKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjb250ZXh0LnJlc3RvcmUoKTtcbiAgICB9XG4gICAgZHJhd1RyYW5zZm9ybXNNYXJrZXIoY29udGV4dCwgcG9zaXRpb24sIHhDb2xvdXIsIHlDb2xvdXIsIGxpbmVXaWR0aCwgc2l6ZSkge1xuICAgICAgICBjb250ZXh0LnNhdmUoKTtcbiAgICAgICAgY29udGV4dC5saW5lV2lkdGggPSBsaW5lV2lkdGg7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlU3R5bGUgPSB4Q29sb3VyO1xuICAgICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICBjb250ZXh0Lm1vdmVUbyhwb3NpdGlvbi54LCBwb3NpdGlvbi55KTtcbiAgICAgICAgY29udGV4dC5saW5lVG8ocG9zaXRpb24ueCArIHNpemUsIHBvc2l0aW9uLnkpO1xuICAgICAgICBjb250ZXh0LnN0cm9rZSgpO1xuICAgICAgICBjb250ZXh0LnN0cm9rZVN0eWxlID0geUNvbG91cjtcbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKTtcbiAgICAgICAgY29udGV4dC5tb3ZlVG8ocG9zaXRpb24ueCwgcG9zaXRpb24ueSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKHBvc2l0aW9uLngsIHBvc2l0aW9uLnkgKyBzaXplKTtcbiAgICAgICAgY29udGV4dC5zdHJva2UoKTtcbiAgICAgICAgY29udGV4dC5yZXN0b3JlKCk7XG4gICAgfVxuICAgIGRyYXdDcm9zcyhjb250ZXh0LCBwb3NpdGlvbiwgY29sb3VyLCBsaW5lV2lkdGgsIHNpemUpIHtcbiAgICAgICAgY29udGV4dC5zYXZlKCk7XG4gICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gbGluZVdpZHRoO1xuICAgICAgICBjb25zdCBoYWxmU2l6ZSA9IE1hdGguY2VpbChzaXplIC8gMik7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlU3R5bGUgPSBjb2xvdXI7XG4gICAgICAgIGNvbnRleHQuYmVnaW5QYXRoKCk7XG4gICAgICAgIGNvbnRleHQubW92ZVRvKHBvc2l0aW9uLnggLSBoYWxmU2l6ZSwgcG9zaXRpb24ueSAtIGhhbGZTaXplKTtcbiAgICAgICAgY29udGV4dC5saW5lVG8ocG9zaXRpb24ueCArIGhhbGZTaXplLCBwb3NpdGlvbi55ICsgaGFsZlNpemUpO1xuICAgICAgICBjb250ZXh0Lm1vdmVUbyhwb3NpdGlvbi54IC0gaGFsZlNpemUsIHBvc2l0aW9uLnkgKyBoYWxmU2l6ZSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKHBvc2l0aW9uLnggKyBoYWxmU2l6ZSwgcG9zaXRpb24ueSAtIGhhbGZTaXplKTtcbiAgICAgICAgY29udGV4dC5zdHJva2UoKTtcbiAgICAgICAgY29udGV4dC5yZXN0b3JlKCk7XG4gICAgfVxufVxuZXhwb3J0cy5TcHJpdGUgPSBTcHJpdGU7XG5TcHJpdGUuREVGQVVMVF9PUFRJT05TID0ge1xuICAgIGRpcmVjdGlvbnM6IFsnZGVmYXVsdCddLFxuICAgIGRlZmF1bHREaXJlY3Rpb246ICdkZWZhdWx0JyxcbiAgICBhbmltYXRpb25zOiB7XG4gICAgICAgIGRlZmF1bHQ6IHtcbiAgICAgICAgICAgICcqJzoge1xuICAgICAgICAgICAgICAgIG5hbWU6ICdkZWZhdWx0JyxcbiAgICAgICAgICAgICAgICBmcmFtZUNvdW50OiAxLFxuICAgICAgICAgICAgICAgIGZyYW1lUmF0ZTogMSxcbiAgICAgICAgICAgICAgICBtb2RlOiBTcHJpdGVBbmltYXRpb25SZXBlYXRNb2RlLlBsYXlPbmNlQW5kU3RvcCxcbiAgICAgICAgICAgIH0sXG4gICAgICAgIH0sXG4gICAgfSxcbiAgICBkZWZhdWx0QW5pbWF0aW9uOiAnZGVmYXVsdCcsXG59O1xuU3ByaXRlLkRFQlVHX0JPVU5ESU5HX0JPWF9DT0xPVVIgPSAnZ3JlZW4nO1xuU3ByaXRlLkRFQlVHX0JPVU5ESU5HX0JPWF9MSU5FX1dJRFRIID0gMjtcblNwcml0ZS5ERUJVR19UUkFOU0ZPUk1TX0NPTE9VUl9YID0gJ3JlZCc7XG5TcHJpdGUuREVCVUdfVFJBTlNGT1JNU19DT0xPVVJfWSA9ICdvcmFuZ2UnO1xuU3ByaXRlLkRFQlVHX1RSQU5TRk9STVNfTElORV9XSURUSCA9IDE7XG5TcHJpdGUuREVCVUdfVFJBTlNGT1JNU19TSVpFID0gMTA7XG5TcHJpdGUuREVCVUdfQVRUQUNITUVOVF9QT0lOVF9DT0xPVVIgPSAnYmx1ZSc7XG5TcHJpdGUuREVCVUdfQVRUQUNITUVOVF9QT0lOVF9MSU5FX1dJRFRIID0gMjtcblNwcml0ZS5ERUJVR19BVFRBQ0hNRU5UX1BPSU5UX1NJWkUgPSA1O1xuLy8jIHNvdXJjZU1hcHBpbmdVUkw9ZGF0YTphcHBsaWNhdGlvbi9qc29uO2Jhc2U2NCxleUoyWlhKemFXOXVJam96TENKbWFXeGxJam9pYVc1a1pYZ3Vhbk1pTENKemIzVnlZMlZTYjI5MElqb2lJaXdpYzI5MWNtTmxjeUk2V3lJdUxpOXBibVJsZUM1MGN5SmRMQ0p1WVcxbGN5STZXMTBzSW0xaGNIQnBibWR6SWpvaU96czdRVUZCUVN3clEwRkJORU03UVVGblNUVkRMRWxCUVZrc2VVSkJaVmc3UVVGbVJDeFhRVUZaTEhsQ1FVRjVRanRKUVVOdVF6czdUMEZGUnp0SlFVTklMRFpGUVVGVkxFTkJRVUU3U1VGRlZqczdUMEZGUnp0SlFVTklMQ3RHUVVGbExFTkJRVUU3U1VGRlpqczdUMEZGUnp0SlFVTklMR2xIUVVGblFpeERRVUZCTzBGQlEyeENMRU5CUVVNc1JVRm1WeXg1UWtGQmVVSXNSMEZCZWtJc2FVTkJRWGxDTEV0QlFYcENMR2xEUVVGNVFpeFJRV1Z3UXp0QlFUQkdSQ3hOUVVGaExFMUJRVTA3U1VGdlJHcENMRmxCUVcxQ0xFOUJRV2RET3p0UlFXWTFReXhoUVVGUkxFZEJRVkVzU1VGQlFTeFRRVUZITEVkQlFVVXNRMEZCUXp0UlFVTjBRaXhUUVVGSkxFZEJRVkVzU1VGQlFTeFRRVUZITEVkQlFVVXNRMEZCUXp0UlFVVnNRaXhYUVVGTkxFZEJRVkVzU1VGQlFTeFRRVUZITEVkQlFVVXNRMEZCUXp0UlFVTndRaXhWUVVGTExFZEJRVmNzUTBGQlF5eERRVUZETzFGQlEyeENMR0ZCUVZFc1IwRkJWeXhEUVVGRExFTkJRVU03VVVGTGNFSXNORUpCUVhWQ0xFZEJRV3RETEVsQlFVa3NRMEZCUXp0UlFVTTVSQ3d3UWtGQmNVSXNSMEZCWjBNc1NVRkJTU3hEUVVGRE8xRkJRekZFTEdsQ1FVRlpMRWRCUVdkRUxFbEJRVWtzUTBGQlF6dFJRVU5xUlN3MFFrRkJkVUlzUjBGQmIwTXNTVUZCU1N4RFFVRkRPMUZCUjNSRkxFMUJRVTBzWVVGQllTeEhRVUZITEUxQlFVMHNRMEZCUXl4TlFVRk5MRU5CUTJwRExFVkJRVVVzUlVGRFJpeE5RVUZOTEVOQlFVTXNaVUZCWlN4RlFVTjBRaXhQUVVGUExHRkJRVkFzVDBGQlR5eGpRVUZRTEU5QlFVOHNSMEZCU1N4RlFVRkZMRU5CUTJRc1EwRkJRenRSUVVWR0xFbEJRVWtzUTBGQlF5eGhRVUZoTEVOQlFVTXNTMEZCU3l4SlFVRkpMR0ZCUVdFc1EwRkJReXhMUVVGTExFdEJRVXNzU1VGQlNTeEZRVUZGTzFsQlEzaEVMR0ZCUVdFc1EwRkJReXhMUVVGTExFZEJRVWM3WjBKQlEzQkNMRzlDUVVGdlFpeEZRVUZGTEVOQlFVTXNRMEZCUXl4aFFVRmhMRU5CUVVNc1MwRkJTenRuUWtGRE0wTXNjVUpCUVhGQ0xFVkJRVVVzUTBGQlF5eERRVUZETEdGQlFXRXNRMEZCUXl4TFFVRkxPMmRDUVVNMVF5eHZRa0ZCYjBJc1JVRkJSU3hEUVVGRExFTkJRVU1zWVVGQllTeERRVUZETEV0QlFVczdZVUZETlVNc1EwRkJRenRUUVVOSU8xRkJSVVFzU1VGQlNTeERRVUZETEU5QlFVOHNSMEZCUnl4aFFVRnZReXhEUVVGRE8xRkJSWEJFTEVsQlFVa3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhSUVVGUkxFVkJRVVU3V1VGRGVrSXNTVUZCU1N4RFFVRkRMRkZCUVZFc1IwRkJSeXhUUVVGSExFTkJRVU1zUjBGQlJ5eERRVUZETEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1VVRkJVU3hEUVVGRExFTkJRVU03VTBGRGFFUTdVVUZGUkN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zU1VGQlNTeEZRVUZGTzFsQlEzSkNMRWxCUVVrc1EwRkJReXhKUVVGSkxFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZCUXl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFbEJRVWtzUTBGQlF5eERRVUZETzFOQlEzaERPMkZCUVUwN1dVRkRUQ3h6UkVGQmMwUTdXVUZEZEVRc1NVRkJTU3hKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEV0QlFVc3NSVUZCUlR0blFrRkRkRUlzU1VGQlNTeERRVUZETEVsQlFVa3NSMEZCUnl4SlFVRkJMRk5CUVVjc1JVRkRZaXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEV0QlFVc3NRMEZCUXl4TFFVRkxMRVZCUTNoQ0xFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNTMEZCU3l4RFFVRkRMRTFCUVUwc1EwRkRNVUlzUTBGQlF6dGhRVU5JTzJsQ1FVRk5PMmRDUVVOTUxIRkZRVUZ4UlR0blFrRkRja1VzTmtSQlFUWkVPMmRDUVVNM1JDeE5RVUZOTERCQ1FVRXdRaXhIUVVGSExFMUJRVTBzUTBGQlF5eE5RVUZOTEVOQlF6bERMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVlVGQlZTeERRVUZETEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1owSkJRV2RDTEVOQlFVTXNRMEZEZGtRc1EwRkJReXhEUVVGRExFTkJRVU1zUTBGQlF6dG5Ra0ZEVEN4SlFVTkZMREJDUVVFd1FqdHZRa0ZETVVJc1EwRkJReXhOUVVGQkxFMUJRVUVzTUVKQlFUQkNMRU5CUVVNc1RVRkJUU3d3UTBGQlJTeE5RVUZOTEcxRFFVRkpMRU5CUVVNc1EwRkJReXhIUVVGSExFTkJRVU1zUlVGRGNFUTdiMEpCUTBFc1NVRkJTU3hEUVVGRExFbEJRVWtzUjBGQlJ5eEpRVUZCTEZOQlFVY3NSVUZEWWl3d1FrRkJNRUlzUTBGQlF5eE5RVUZQTEVOQlFVTXNRMEZCUXl4RFFVRkRMRU5CUVVNc1MwRkJTeXhGUVVNelF5d3dRa0ZCTUVJc1EwRkJReXhOUVVGUExFTkJRVU1zUTBGQlF5eERRVUZETEVOQlFVTXNUVUZCVFN4RFFVTTNReXhEUVVGRE8ybENRVU5JTzJGQlEwWTdXVUZGUkN4eFEwRkJjVU03VTBGRGRFTTdVVUZGUkN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVFVGQlRTeEZRVUZGTzFsQlEzWkNMRWxCUVVrc1EwRkJReXhOUVVGTkxFZEJRVWNzVTBGQlJ5eERRVUZETEVkQlFVY3NRMEZCUXl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eERRVUZETzFOQlF6VkRPMkZCUVUwN1dVRkRUQ3h2UkVGQmIwUTdXVUZEY0VRc1NVRkJTU3hEUVVGRExFMUJRVTBzUjBGQlJ5eFRRVUZITEVOQlFVTXNSMEZCUnl4RFFVRkRMRWxCUVVrc1EwRkJReXhKUVVGSkxFVkJRVVVzUjBGQlJ5eERRVUZETEVOQlFVTTdVMEZEZGtNN1VVRkZSQ3hKUVVGSkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNTMEZCU3l4RlFVRkZPMWxCUTNSQ0xFbEJRVWtzUTBGQlF5eExRVUZMTEVkQlFVY3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhMUVVGTExFTkJRVU03VTBGRGFrTTdVVUZGUkN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVVVGQlVTeEZRVUZGTzFsQlEzcENMRWxCUVVrc1EwRkJReXhSUVVGUkxFZEJRVWNzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4UlFVRlJMRU5CUVVNN1UwRkRka003VVVGRlJDeHBRMEZCYVVNN1VVRkRha01zU1VGQlNTeERRVUZETEZWQlFWVXNSMEZCUnl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExHZENRVUZuUWl4RFFVRkRPMUZCUTJoRUxFbEJRMFVzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4VlFVRlZMRU5CUVVNc1RVRkJUU3hMUVVGTExFTkJRVU03V1VGRGNFTXNRMEZCUXl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGVkJRVlVzUTBGQlF5eFJRVUZSTEVOQlFVTXNTVUZCU1N4RFFVRkRMRlZCUVZVc1EwRkJReXhGUVVOc1JEdFpRVU5CTEUxQlFVMHNTVUZCU1N4TFFVRkxMRU5CUVVNc2MwSkJRWE5DTEVsQlFVa3NRMEZCUXl4VlFVRlZMRWRCUVVjc1EwRkJReXhEUVVGRE8xTkJRek5FTzFGQlJVUXNhVU5CUVdsRE8xRkJRMnBETEVsQlFVa3NRMEZCUXl4VlFVRlZMRWRCUVVjc1NVRkJTU3hEUVVGRExFOUJRVThzUTBGQlF5eG5Ra0ZCWjBJc1EwRkJRenRSUVVOb1JDeE5RVUZOTEZWQlFWVXNSMEZCUnl4TlFVRk5MRU5CUVVNc1NVRkJTU3hEUVVGRExFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNWVUZCVlN4RFFVRkRMRU5CUVVNN1VVRkRlRVFzU1VGRFJTeFZRVUZWTEVOQlFVTXNUVUZCVFN4TFFVRkxMRU5CUVVNN1dVRkRka0lzUTBGQlF5eFZRVUZWTEVOQlFVTXNVVUZCVVN4RFFVRkRMRWxCUVVrc1EwRkJReXhWUVVGVkxFTkJRVU1zUlVGRGNrTTdXVUZEUVN4TlFVRk5MRWxCUVVrc1MwRkJTeXhEUVVGRExITkNRVUZ6UWl4SlFVRkpMRU5CUVVNc1ZVRkJWU3hIUVVGSExFTkJRVU1zUTBGQlF6dFRRVU16UkR0UlFVVkVMR2RGUVVGblJUdFJRVU5vUlN4blEwRkJaME03VVVGRGFFTXNTMEZCU3l4TlFVRk5MRk5CUVZNc1NVRkJTU3hOUVVGTkxFTkJRVU1zU1VGQlNTeERRVUZETEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1ZVRkJWU3hEUVVGRExFVkJRVVU3V1VGRE5VUXNTMEZCU3l4TlFVRk5MRk5CUVZNc1NVRkJTU3hOUVVGTkxFTkJRVU1zU1VGQlNTeERRVUZETEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1ZVRkJWU3hEUVVGRExGTkJRVk1zUTBGQlF5eERRVUZETEVWQlFVVTdaMEpCUTNaRkxFbEJRMFVzU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4VlFVRlZMRU5CUVVNc1UwRkJVeXhEUVVGRExFTkJRVU1zVTBGQlV5eERRVUZETEVOQlFVTXNkMEpCUVhkQ0xFVkJRM1JGTzI5Q1FVTkJMRXRCUVVzc1RVRkJUU3hsUVVGbExFbEJRVWtzVFVGQlRTeERRVUZETEVsQlFVa3NRMEZEZGtNc1NVRkJTVHQ1UWtGRFJDeFBRVUZQTzNsQ1FVTlFMRlZCUVZVc1EwRkJReXhUUVVGVExFTkJRVU1zUTBGQlF5eFRRVUZUTEVOQlFVTTdlVUpCUTJoRExIZENRVUY1UWl4RFFVTTNRaXhGUVVGRk8zZENRVU5FTEVsQlFVazdOa0pCUTBRc1QwRkJUenMyUWtGRFVDeFZRVUZWTEVOQlFVTXNVMEZCVXl4RFFVRkRMRU5CUVVNc1UwRkJVeXhEUVVGRE96WkNRVU5vUXl4M1FrRkJlVUlzUTBGQlF5eGxRVUZsTEVOQlFVTTdOa0pCUXpGRExFbEJRVWtzUTBGRFNDeERRVUZETEVOQlFVTXNSVUZCUlN4RFFVRkRMRVZCUVVVc1JVRkJSU3hEUVVGRExFTkJRVU1zUTBGQlF5eExRVUZMTEVkQlFVY3NRMEZCUXl4RFFVRkRMRXRCUVVzc1EwRkROVUlzUTBGQlF6dHhRa0ZEVER0cFFrRkRSanRoUVVOR08xTkJRMFk3U1VGRFNDeERRVUZETzBsQlJVUXNTVUZCVnl4VFFVRlRPMUZCUTJ4Q0xFOUJRVThzU1VGQlNTeERRVUZETEZWQlFWVXNRMEZCUXp0SlFVTjZRaXhEUVVGRE8wbEJSVVFzU1VGQlZ5eFRRVUZUTEVOQlFVTXNTMEZCWVR0UlFVTm9ReXhKUVVGSkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNWVUZCVlN4RFFVRkRMRkZCUVZFc1EwRkJReXhMUVVGTExFTkJRVU1zUlVGQlJUdFpRVU16UXl4SlFVRkpMRU5CUVVNc1ZVRkJWU3hIUVVGSExFdEJRVXNzUTBGQlF6dFRRVU42UWp0SlFVTklMRU5CUVVNN1NVRkZSQ3hKUVVGWExGTkJRVk03VVVGRGJFSXNUMEZCVHl4SlFVRkpMRU5CUVVNc1ZVRkJWU3hEUVVGRE8wbEJRM3BDTEVOQlFVTTdTVUZGUkN4SlFVRlhMRk5CUVZNc1EwRkJReXhMUVVGaE8xRkJRMmhETEVsQlFVa3NUVUZCVFN4RFFVRkRMRWxCUVVrc1EwRkJReXhKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZWQlFWVXNRMEZCUXl4RFFVRkRMRkZCUVZFc1EwRkJReXhMUVVGTExFTkJRVU1zUlVGQlJUdFpRVU40UkN4SlFVRkpMRU5CUVVNc1ZVRkJWU3hIUVVGSExFdEJRVXNzUTBGQlF6dFRRVU42UWp0SlFVTklMRU5CUVVNN1NVRkZUU3hoUVVGaE8xRkJRMnhDTEVsQlFVa3NTVUZCU1N4RFFVRkRMSEZDUVVGeFFpeEZRVUZGTzFsQlF6bENMRWxCUVVrc1EwRkJReXh4UWtGQmNVSXNRMEZCUXl4UFFVRlBMRWRCUVVjc1NVRkJTU3hEUVVGRE8xTkJRek5ETzBsQlEwZ3NRMEZCUXp0SlFVVk5MR05CUVdNN1VVRkRia0lzU1VGQlNTeEpRVUZKTEVOQlFVTXNjVUpCUVhGQ0xFVkJRVVU3V1VGRE9VSXNTVUZCU1N4RFFVRkRMSEZDUVVGeFFpeERRVUZETEU5QlFVOHNSMEZCUnl4TFFVRkxMRU5CUVVNN1UwRkROVU03U1VGRFNDeERRVUZETzBsQlJVMHNZMEZCWXp0UlFVTnVRaXhKUVVGSkxFbEJRVWtzUTBGQlF5eHhRa0ZCY1VJc1JVRkJSVHRaUVVNNVFpeEpRVUZKTEVOQlFVTXNjVUpCUVhGQ0xFTkJRVU1zV1VGQldTeEhRVUZITEVOQlFVTXNRMEZCUXp0WlFVTTFReXhKUVVGSkxFTkJRVU1zY1VKQlFYRkNMRU5CUVVNc1owSkJRV2RDTEVkQlFVY3NRMEZCUXl4RFFVRkRPMU5CUTJwRU8wbEJRMGdzUTBGQlF6dEpRVVZOTEd0Q1FVRnJRaXhEUVVGRExFbEJRVms3TzFGQlEzQkRMRTlCUVU4c1RVRkJRU3hOUVVGQkxFbEJRVWtzUTBGQlF5eDFRa0ZCZFVJc01FTkJRVWNzU1VGQlNTeERRVUZETEcxRFFVRkpMRWxCUVVrc1EwRkJRenRKUVVOMFJDeERRVUZETzBsQlJVMHNUVUZCVFN4RFFVRkRMRVZCUVZVN1VVRkRkRUlzU1VGQlNTeERRVUZETEhWQ1FVRjFRaXhIUVVGSExFbEJRVWtzUTBGQlF5eHpRa0ZCYzBJc1JVRkJSU3hEUVVGRE8xRkJRemRFTEVsQlFVa3NRMEZCUXl4eFFrRkJjVUlzUjBGQlJ5eEpRVUZKTEVOQlFVTXNiMEpCUVc5Q0xFTkJRVU1zUlVGQlJTeERRVUZETEVOQlFVTTdVVUZETTBRc1NVRkJTU3hEUVVGRExGbEJRVmtzUjBGQlJ5eEpRVUZKTEVOQlFVTXNWMEZCVnl4RlFVRkZMRU5CUVVNN1VVRkRka01zU1VGQlNTeERRVUZETEhWQ1FVRjFRaXhIUVVGSExFbEJRVWtzUTBGQlF5eHpRa0ZCYzBJc1JVRkJSU3hEUVVGRE8wbEJReTlFTEVOQlFVTTdTVUZGVHl4elFrRkJjMEk3VVVGRE5VSXNTVUZCU1N4RFFVRkRMRU5CUVVNc1NVRkJTU3hEUVVGRExGVkJRVlVzU1VGQlNTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRlZCUVZVc1EwRkJReXhGUVVGRk8xbEJRMnBFTEUxQlFVMHNTVUZCU1N4TFFVRkxMRU5CUVVNc2MwSkJRWE5DTEVsQlFVa3NRMEZCUXl4VlFVRlZMRWRCUVVjc1EwRkJReXhEUVVGRE8xTkJRek5FTzFGQlJVUXNUVUZCVFN4VlFVRlZMRWRCUVVjc1RVRkJUU3hEUVVGRExFbEJRVWtzUTBGQlF5eEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRlZCUVZVc1EwRkJReXhKUVVGSkxFTkJRVU1zVlVGQlZTeERRVUZETEVOQlFVTXNRMEZCUXp0UlFVTjZSU3hKUVVGSkxGVkJRVlVzUTBGQlF5eE5RVUZOTEV0QlFVc3NRMEZCUXl4RlFVRkZPMWxCUXpOQ0xFMUJRVTBzU1VGQlNTeExRVUZMTEVOQlEySXNNRU5CUVRCRExFbEJRVWtzUTBGQlF5eFZRVUZWTEVkQlFVY3NRMEZETjBRc1EwRkJRenRUUVVOSU8xRkJSVVFzU1VGQlNTeEpRVUZKTEVOQlFVTXNWVUZCVlN4SlFVRkpMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVlVGQlZTeERRVUZETEVsQlFVa3NRMEZCUXl4VlFVRlZMRU5CUVVNc1JVRkJSVHRaUVVNdlJDeFBRVUZQTEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1ZVRkJWU3hEUVVGRExFbEJRVWtzUTBGQlF5eFZRVUZWTEVOQlFVTXNRMEZCUXl4SlFVRkpMRU5CUVVNc1ZVRkJWU3hEUVVGRExFTkJRVU03VTBGRGJFVTdVVUZGUkN4SlFVRkpMRWRCUVVjc1NVRkJTU3hKUVVGSkxFTkJRVU1zVDBGQlR5eERRVUZETEZWQlFWVXNRMEZCUXl4SlFVRkpMRU5CUVVNc1ZVRkJWU3hEUVVGRExFVkJRVVU3V1VGRGJrUXNUMEZCVHl4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExGVkJRVlVzUTBGQlF5eEpRVUZKTEVOQlFVTXNWVUZCVlN4RFFVRkRMRU5CUVVNc1IwRkJSeXhEUVVGRExFTkJRVU03VTBGRGRFUTdVVUZGUkN4UFFVRlBMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zVlVGQlZTeERRVUZETEVsQlFVa3NRMEZCUXl4VlFVRlZMRU5CUVVNc1EwRkJReXhWUVVGVkxFTkJRVU1zUTBGQlF5eERRVUZETEVOQlFVTXNRMEZCUXp0SlFVTnFSU3hEUVVGRE8wbEJSVThzYjBKQlFXOUNMRU5CUVVNc1JVRkJWVHRSUVVOeVF5eEpRVU5GTEVOQlFVTXNTVUZCU1N4RFFVRkRMSFZDUVVGMVFqdFpRVU0zUWl4RFFVRkRMRWxCUVVrc1EwRkJReXh4UWtGQmNVSXNSVUZETTBJN1dVRkRRU3hQUVVGUE8yZENRVU5NTEU5QlFVOHNSVUZCUlN4SlFVRkpPMmRDUVVOaUxGbEJRVmtzUlVGQlJTeERRVUZETzJkQ1FVTm1MR2RDUVVGblFpeEZRVUZGTEVOQlFVTTdZVUZEY0VJc1EwRkJRenRUUVVOSU8xRkJSVVFzU1VGQlNTeEpRVUZKTEVOQlFVTXNjVUpCUVhGQ0xFTkJRVU1zVDBGQlR5eEZRVUZGTzFsQlEzUkRMRTFCUVUwc1UwRkJVeXhIUVVGSExFTkJRVU1zUjBGQlJ5eEpRVUZKTEVOQlFVTXNkVUpCUVhWQ0xFTkJRVU1zVTBGQlV5eERRVUZETzFsQlF6ZEVMRWxCUVVrc1EwRkJReXh4UWtGQmNVSXNRMEZCUXl4blFrRkJaMElzU1VGQlNTeEZRVUZGTEVOQlFVTTdXVUZGYkVRc1NVRkJTU3hKUVVGSkxFTkJRVU1zY1VKQlFYRkNMRU5CUVVNc1owSkJRV2RDTEVkQlFVY3NVMEZCVXl4RlFVRkZPMmRDUVVNelJDeE5RVUZOTEZWQlFWVXNSMEZCUnl4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVOQlFVTXNWVUZCVlN4RFFVRkRPMmRDUVVNelJDeEpRVUZKTEVOQlFVTXNjVUpCUVhGQ0xFTkJRVU1zV1VGQldTeEZRVUZGTEVOQlFVTTdaMEpCUXpGRExFbEJRVWtzUTBGQlF5eHhRa0ZCY1VJc1EwRkJReXhuUWtGQlowSXNSMEZCUnl4RFFVRkRMRU5CUVVNN1owSkJSV2hFTEVsQlFVa3NTVUZCU1N4RFFVRkRMSEZDUVVGeFFpeERRVUZETEZsQlFWa3NSMEZCUnl4VlFVRlZMRVZCUVVVN2IwSkJRM2hFTEZGQlFWRXNTVUZCU1N4RFFVRkRMSFZDUVVGMVFpeERRVUZETEVsQlFVa3NSVUZCUlR0M1FrRkRla01zUzBGQlN5eDVRa0ZCZVVJc1EwRkJReXhuUWtGQlowSTdORUpCUXpkRExFbEJRVWtzUTBGQlF5eHhRa0ZCY1VJc1EwRkJReXhQUVVGUExFZEJRVWNzUzBGQlN5eERRVUZET3pSQ1FVTXpReXhKUVVGSkxFTkJRVU1zY1VKQlFYRkNMRU5CUVVNc1dVRkJXU3hIUVVGSExFTkJRVU1zUTBGQlF6czBRa0ZETlVNc1RVRkJUVHQzUWtGRlVpeExRVUZMTEhsQ1FVRjVRaXhEUVVGRExHVkJRV1U3TkVKQlF6VkRMRWxCUVVrc1EwRkJReXh4UWtGQmNVSXNRMEZCUXl4UFFVRlBMRWRCUVVjc1MwRkJTeXhEUVVGRE96UkNRVU16UXl4SlFVRkpMRU5CUVVNc2NVSkJRWEZDTEVOQlFVTXNXVUZCV1N4SFFVRkhMRlZCUVZVc1IwRkJSeXhEUVVGRExFTkJRVU03TkVKQlEzcEVMRTFCUVUwN2QwSkJSVklzUzBGQlN5eDVRa0ZCZVVJc1EwRkJReXhOUVVGTk96UkNRVU51UXl4SlFVRkpMRU5CUVVNc2NVSkJRWEZDTEVOQlFVTXNXVUZCV1N4SFFVRkhMRU5CUVVNc1EwRkJRenMwUWtGRE5VTXNUVUZCVFR0eFFrRkRWRHRwUWtGRFJqdGhRVU5HTzFOQlEwWTdVVUZGUkN4UFFVRlBMRWxCUVVrc1EwRkJReXh4UWtGQmNVSXNRMEZCUXp0SlFVTndReXhEUVVGRE8wbEJSVThzVjBGQlZ6czdVVUZEYWtJc1NVRkRSU3hEUVVGRExFbEJRVWtzUTBGQlF5eDFRa0ZCZFVJN1dVRkROMElzUTBGQlF5eEpRVUZKTEVOQlFVTXNjVUpCUVhGQ0xFVkJRek5DTzFsQlEwRXNUMEZCVHl4SlFVRkpMRU5CUVVNN1UwRkRZanRSUVVWRUxFbEJRMFVzUTBGQlF5eEpRVUZKTEVOQlFVTXNkVUpCUVhWQ0xFTkJRVU1zVFVGQlRUdFpRVU53UXl4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVOQlFVTXNUVUZCVFN4RFFVRkRMRTFCUVUwc1MwRkJTeXhEUVVGRExFVkJRMmhFTzFsQlEwRXNUMEZCVHl4TlFVRkJMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zUzBGQlN5eHRRMEZCU1N4SlFVRkpMRU5CUVVNN1UwRkRia003VVVGRlJDeFBRVUZQTEUxQlFVRXNUVUZCUVN4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVOQlFVTXNUVUZCVFN4RFFVTjRReXhKUVVGSkxFTkJRVU1zY1VKQlFYRkNMRU5CUVVNc1dVRkJXU3hEUVVONFF5eHRRMEZCU1N4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFdEJRVXNzYlVOQlFVa3NTVUZCU1N4RFFVRkRPMGxCUTJ4RExFTkJRVU03U1VGRlR5eHpRa0ZCYzBJN1VVRkROVUlzU1VGRFJTeERRVUZETEVsQlFVa3NRMEZCUXl4UFFVRlBMRU5CUVVNc1owSkJRV2RDTzFsQlF6bENMRWxCUVVrc1EwRkJReXhQUVVGUExFTkJRVU1zWjBKQlFXZENMRU5CUVVNc1RVRkJUU3hMUVVGTExFTkJRVU1zUlVGRE1VTTdXVUZEUVN4UFFVRlBMRWxCUVVrc1EwRkJRenRUUVVOaU8xRkJSVVFzU1VGQlNTeERRVUZETEVsQlFVa3NRMEZCUXl4MVFrRkJkVUlzUlVGQlJUdFpRVU5xUXl4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVkQlFVY3NUVUZCVFN4RFFVRkRMRmRCUVZjc1EwRkRMME1zU1VGQlNTeERRVUZETEU5QlFVOHNRMEZCUXl4blFrRkJaMElzUTBGQlF5eEhRVUZITEVOQlFVTXNaVUZCWlN4RFFVRkRMRVZCUVVVc1EwRkJRenRuUWtGRGJrUXNaVUZCWlN4RFFVRkRMRWxCUVVrN1owSkJRM0JDTEdWQlFXVXNRMEZCUXl4TlFVRk5PMkZCUTNaQ0xFTkJRVU1zUTBGRFNDeERRVUZETzFOQlEwZzdVVUZGUkN4SlFVTkZMRWxCUVVrc1EwRkJReXgxUWtGQmRVSTdXVUZETlVJc1NVRkJTU3hEUVVGRExIVkNRVUYxUWl4RFFVRkRMSGRDUVVGM1FqdFpRVU55UkN4SlFVRkpMRU5CUVVNc2NVSkJRWEZDTEVWQlF6RkNPMWxCUTBFc1MwRkJTeXhOUVVGTkxFbEJRVWtzU1VGQlNTeE5RVUZOTEVOQlFVTXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkJReXgxUWtGQmRVSXNRMEZCUXl4RlFVRkZPMmRDUVVNMVJDeEpRVU5GTEVsQlFVa3NTVUZCU1N4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVOQlFVTXNkMEpCUVhkQ08yOUNRVU0zUkN4SlFVRkpMRU5CUVVNc2RVSkJRWFZDTEVOQlFVTXNkMEpCUVhkQ0xFTkJRVU1zU1VGQlNTeERRVUZETEVOQlFVTXNUVUZCVFN4SFFVRkhMRU5CUVVNc1JVRkRkRVU3YjBKQlEwRXNUVUZCVFN4blFrRkJaMElzUjBGQlJ5eEpRVUZKTEVOQlFVTXNiMEpCUVc5Q0xFTkJRMmhFTEVsQlFVa3NRMEZCUXl4MVFrRkJkVUlzUTBGQlF5eDNRa0ZCZDBJc1EwRkJReXhKUVVGSkxFTkJRVU1zUlVGRE0wUXNTVUZCU1N4RFFVRkRMSEZDUVVGeFFpeERRVUZETEZsQlFWa3NRMEZEZUVNc1EwRkJRenR2UWtGRFJpeEpRVUZKTEVOQlFVTXNkVUpCUVhWQ0xFTkJRVU1zU1VGQlNTeERRVUZETEVkQlFVY3NaMEpCUVdkQ0xFTkJRVU1zVFVGQlRTeERRVUZETzJsQ1FVTTVSRHRoUVVOR08xTkJRMFk3VVVGRlJDeFBRVUZQTEVsQlFVa3NRMEZCUXl4MVFrRkJkVUlzUTBGQlF6dEpRVU4wUXl4RFFVRkRPMGxCUlU4c2IwSkJRVzlDTEVOQlF6RkNMRk5CUVRCRExFVkJRekZETEZsQlFXOUNPMUZCUlhCQ0xFMUJRVTBzUzBGQlN5eEhRVUZITEVOQlFVTXNSMEZCUnl4VFFVRlRMRU5CUVVNc1EwRkJReXhQUVVGUExFVkJRVVVzUTBGQlF5eEpRVUZKTEVOQlEzcERMRkZCUVZFc1EwRkJReXhGUVVGRkxFTkJRVU1zVVVGQlVTeERRVUZETEV0QlFVc3NTVUZCU1N4WlFVRlpMRU5CUXpORExFTkJRVU03VVVGRlJpeEpRVUZKTEVOQlFVTXNTMEZCU3l4RlFVRkZPMWxCUTFZc1QwRkJUeXhUUVVGVExFTkJRVU1zVTBGQlV5eERRVUZETEUxQlFVMHNSMEZCUnl4RFFVRkRMRU5CUVVNc1EwRkJRenRUUVVONFF6dFJRVVZFTEU5QlFVOHNTMEZCU3l4RFFVRkRPMGxCUTJZc1EwRkJRenRKUVVWTkxFbEJRVWtzUTBGQlF5eFBRVUZwUXp0UlFVTXpReXhQUVVGUExFTkJRVU1zU1VGQlNTeEZRVUZGTEVOQlFVTTdVVUZEWml4UFFVRlBMRU5CUVVNc1UwRkJVeXhEUVVObUxFbEJRVWtzUTBGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXl4RlFVTm1MRWxCUVVrc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eERRVU5vUWl4RFFVRkRPMUZCUTBZc1QwRkJUeXhEUVVGRExFdEJRVXNzUTBGQlF5eEpRVUZKTEVOQlFVTXNTMEZCU3l4RlFVRkZMRWxCUVVrc1EwRkJReXhMUVVGTExFTkJRVU1zUTBGQlF6dFJRVU4wUXl4UFFVRlBMRU5CUVVNc1RVRkJUU3hEUVVGRExFbEJRVWtzUTBGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXp0UlFVVTVRaXhKUVVGSkxFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVMEZCVXl4RlFVRkZPMWxCUXpGQ0xFbEJRVWtzUTBGQlF5eFBRVUZQTEVOQlFVTXNVMEZCVXl4RFFVRkRMRTlCUVU4c1JVRkJSU3hKUVVGSkxFTkJRVU1zUTBGQlF6dFRRVU4yUXp0UlFVVkVMRWxCUVVrc1NVRkJTU3hEUVVGRExGbEJRVmtzUlVGQlJUdFpRVU55UWl4UFFVRlBMRU5CUVVNc1UwRkJVeXhEUVVObUxFbEJRVWtzUTBGQlF5eFpRVUZaTEVWQlEycENMRU5CUVVNc1NVRkJTU3hEUVVGRExFMUJRVTBzUTBGQlF5eERRVUZETEVWQlEyUXNRMEZCUXl4SlFVRkpMRU5CUVVNc1RVRkJUU3hEUVVGRExFTkJRVU1zUlVGRFpDeEpRVUZKTEVOQlFVTXNXVUZCV1N4RFFVRkRMRXRCUVVzc1JVRkRka0lzU1VGQlNTeERRVUZETEZsQlFWa3NRMEZCUXl4TlFVRk5MRU5CUTNwQ0xFTkJRVU03VTBGRFNEdFJRVVZFTEVsQlFVa3NTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhWUVVGVkxFVkJRVVU3V1VGRE0wSXNTVUZCU1N4RFFVRkRMRTlCUVU4c1EwRkJReXhWUVVGVkxFTkJRVU1zVDBGQlR5eEZRVUZGTEVsQlFVa3NRMEZCUXl4RFFVRkRPMU5CUTNoRE8xRkJSVVFzU1VGQlNTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRXRCUVVzc1EwRkJReXh4UWtGQmNVSXNSVUZCUlR0WlFVTTFReXhQUVVGUExFTkJRVU1zVjBGQlZ5eEhRVUZITEUxQlFVMHNRMEZCUXl4NVFrRkJlVUlzUTBGQlF6dFpRVU4yUkN4UFFVRlBMRU5CUVVNc1UwRkJVeXhIUVVGSExFMUJRVTBzUTBGQlF5dzJRa0ZCTmtJc1EwRkJRenRaUVVONlJDeFBRVUZQTEVOQlFVTXNWVUZCVlN4RFFVTm9RaXhEUVVGRExFbEJRVWtzUTBGQlF5eE5RVUZOTEVOQlFVTXNRMEZCUXl4RlFVTmtMRU5CUVVNc1NVRkJTU3hEUVVGRExFMUJRVTBzUTBGQlF5eERRVUZETEVWQlEyUXNTVUZCU1N4RFFVRkRMRWxCUVVrc1EwRkJReXhEUVVGRExFVkJRMWdzU1VGQlNTeERRVUZETEVsQlFVa3NRMEZCUXl4RFFVRkRMRU5CUTFvc1EwRkJRenRUUVVOSU8xRkJSVVFzU1VGQlNTeEpRVUZKTEVOQlFVTXNUMEZCVHl4RFFVRkRMRXRCUVVzc1EwRkJReXh2UWtGQmIwSXNSVUZCUlR0WlFVTXpReXhKUVVGSkxFTkJRVU1zYjBKQlFXOUNMRU5CUTNaQ0xFOUJRVThzUlVGRFVDeEpRVUZCTEZOQlFVY3NSMEZCUlN4RlFVTk1MRTFCUVUwc1EwRkJReXg1UWtGQmVVSXNSVUZEYUVNc1RVRkJUU3hEUVVGRExIbENRVUY1UWl4RlFVTm9ReXhOUVVGTkxFTkJRVU1zTWtKQlFUSkNMRVZCUTJ4RExFMUJRVTBzUTBGQlF5eHhRa0ZCY1VJc1EwRkROMElzUTBGQlF6dFRRVU5JTzFGQlJVUXNTVUZEUlN4SlFVRkpMRU5CUVVNc1QwRkJUeXhEUVVGRExFdEJRVXNzUTBGQlF5eHZRa0ZCYjBJN1dVRkRka01zU1VGQlNTeERRVUZETEhWQ1FVRjFRaXhGUVVNMVFqdFpRVU5CTEV0QlFVc3NUVUZCVFN4bFFVRmxMRWxCUVVrc1RVRkJUU3hEUVVGRExFMUJRVTBzUTBGQlF5eEpRVUZKTEVOQlFVTXNkVUpCUVhWQ0xFTkJRVU1zUlVGQlJUdG5Ra0ZEZWtVc1NVRkJTU3hEUVVGRExGTkJRVk1zUTBGRFdpeFBRVUZQTEVWQlExQXNaVUZCWlN4RlFVTm1MRTFCUVUwc1EwRkJReXcyUWtGQk5rSXNSVUZEY0VNc1RVRkJUU3hEUVVGRExHbERRVUZwUXl4RlFVTjRReXhOUVVGTkxFTkJRVU1zTWtKQlFUSkNMRU5CUTI1RExFTkJRVU03WVVGRFNEdFRRVU5HTzFGQlJVUXNUMEZCVHl4RFFVRkRMRTlCUVU4c1JVRkJSU3hEUVVGRE8wbEJRM0JDTEVOQlFVTTdTVUZGVHl4dlFrRkJiMElzUTBGRE1VSXNUMEZCYVVNc1JVRkRha01zVVVGQllTeEZRVU5pTEU5QlFXVXNSVUZEWml4UFFVRmxMRVZCUTJZc1UwRkJhVUlzUlVGRGFrSXNTVUZCV1R0UlFVVmFMRTlCUVU4c1EwRkJReXhKUVVGSkxFVkJRVVVzUTBGQlF6dFJRVVZtTEU5QlFVOHNRMEZCUXl4VFFVRlRMRWRCUVVjc1UwRkJVeXhEUVVGRE8xRkJSVGxDTEU5QlFVOHNRMEZCUXl4WFFVRlhMRWRCUVVjc1QwRkJUeXhEUVVGRE8xRkJRemxDTEU5QlFVOHNRMEZCUXl4VFFVRlRMRVZCUVVVc1EwRkJRenRSUVVOd1FpeFBRVUZQTEVOQlFVTXNUVUZCVFN4RFFVRkRMRkZCUVZFc1EwRkJReXhEUVVGRExFVkJRVVVzVVVGQlVTeERRVUZETEVOQlFVTXNRMEZCUXl4RFFVRkRPMUZCUTNaRExFOUJRVThzUTBGQlF5eE5RVUZOTEVOQlFVTXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhKUVVGSkxFVkJRVVVzVVVGQlVTeERRVUZETEVOQlFVTXNRMEZCUXl4RFFVRkRPMUZCUXpsRExFOUJRVThzUTBGQlF5eE5RVUZOTEVWQlFVVXNRMEZCUXp0UlFVVnFRaXhQUVVGUExFTkJRVU1zVjBGQlZ5eEhRVUZITEU5QlFVOHNRMEZCUXp0UlFVTTVRaXhQUVVGUExFTkJRVU1zVTBGQlV5eEZRVUZGTEVOQlFVTTdVVUZEY0VJc1QwRkJUeXhEUVVGRExFMUJRVTBzUTBGQlF5eFJRVUZSTEVOQlFVTXNRMEZCUXl4RlFVRkZMRkZCUVZFc1EwRkJReXhEUVVGRExFTkJRVU1zUTBGQlF6dFJRVU4yUXl4UFFVRlBMRU5CUVVNc1RVRkJUU3hEUVVGRExGRkJRVkVzUTBGQlF5eERRVUZETEVWQlFVVXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhKUVVGSkxFTkJRVU1zUTBGQlF6dFJRVU01UXl4UFFVRlBMRU5CUVVNc1RVRkJUU3hGUVVGRkxFTkJRVU03VVVGRmFrSXNUMEZCVHl4RFFVRkRMRTlCUVU4c1JVRkJSU3hEUVVGRE8wbEJRM0JDTEVOQlFVTTdTVUZGVHl4VFFVRlRMRU5CUTJZc1QwRkJhVU1zUlVGRGFrTXNVVUZCWVN4RlFVTmlMRTFCUVdNc1JVRkRaQ3hUUVVGcFFpeEZRVU5xUWl4SlFVRlpPMUZCUlZvc1QwRkJUeXhEUVVGRExFbEJRVWtzUlVGQlJTeERRVUZETzFGQlJXWXNUMEZCVHl4RFFVRkRMRk5CUVZNc1IwRkJSeXhUUVVGVExFTkJRVU03VVVGRk9VSXNUVUZCVFN4UlFVRlJMRWRCUVVjc1NVRkJTU3hEUVVGRExFbEJRVWtzUTBGQlF5eEpRVUZKTEVkQlFVY3NRMEZCUXl4RFFVRkRMRU5CUVVNN1VVRkRja01zVDBGQlR5eERRVUZETEZkQlFWY3NSMEZCUnl4TlFVRk5MRU5CUVVNN1VVRkROMElzVDBGQlR5eERRVUZETEZOQlFWTXNSVUZCUlN4RFFVRkRPMUZCUTNCQ0xFOUJRVThzUTBGQlF5eE5RVUZOTEVOQlFVTXNVVUZCVVN4RFFVRkRMRU5CUVVNc1IwRkJSeXhSUVVGUkxFVkJRVVVzVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4UlFVRlJMRU5CUVVNc1EwRkJRenRSUVVNM1JDeFBRVUZQTEVOQlFVTXNUVUZCVFN4RFFVRkRMRkZCUVZFc1EwRkJReXhEUVVGRExFZEJRVWNzVVVGQlVTeEZRVUZGTEZGQlFWRXNRMEZCUXl4RFFVRkRMRWRCUVVjc1VVRkJVU3hEUVVGRExFTkJRVU03VVVGRE4wUXNUMEZCVHl4RFFVRkRMRTFCUVUwc1EwRkJReXhSUVVGUkxFTkJRVU1zUTBGQlF5eEhRVUZITEZGQlFWRXNSVUZCUlN4UlFVRlJMRU5CUVVNc1EwRkJReXhIUVVGSExGRkJRVkVzUTBGQlF5eERRVUZETzFGQlF6ZEVMRTlCUVU4c1EwRkJReXhOUVVGTkxFTkJRVU1zVVVGQlVTeERRVUZETEVOQlFVTXNSMEZCUnl4UlFVRlJMRVZCUVVVc1VVRkJVU3hEUVVGRExFTkJRVU1zUjBGQlJ5eFJRVUZSTEVOQlFVTXNRMEZCUXp0UlFVTTNSQ3hQUVVGUExFTkJRVU1zVFVGQlRTeEZRVUZGTEVOQlFVTTdVVUZGYWtJc1QwRkJUeXhEUVVGRExFOUJRVThzUlVGQlJTeERRVUZETzBsQlEzQkNMRU5CUVVNN08wRkJjR1JJTEhkQ1FYRmtRenRCUVhCa2VVSXNjMEpCUVdVc1IwRkJhMEk3U1VGRGRrUXNWVUZCVlN4RlFVRkZMRU5CUVVNc1UwRkJVeXhEUVVGRE8wbEJRM1pDTEdkQ1FVRm5RaXhGUVVGRkxGTkJRVk03U1VGRE0wSXNWVUZCVlN4RlFVRkZPMUZCUTFZc1QwRkJUeXhGUVVGRk8xbEJRMUFzUjBGQlJ5eEZRVUZGTzJkQ1FVTklMRWxCUVVrc1JVRkJSU3hUUVVGVE8yZENRVU5tTEZWQlFWVXNSVUZCUlN4RFFVRkRPMmRDUVVOaUxGTkJRVk1zUlVGQlJTeERRVUZETzJkQ1FVTmFMRWxCUVVrc1JVRkJSU3g1UWtGQmVVSXNRMEZCUXl4bFFVRmxPMkZCUTJoRU8xTkJRMFk3UzBGRFJqdEpRVU5FTEdkQ1FVRm5RaXhGUVVGRkxGTkJRVk03UTBGRE5VSXNRMEZCUXp0QlFVVnpRaXhuUTBGQmVVSXNSMEZCUnl4UFFVRlBMRU5CUVVNN1FVRkRjRU1zYjBOQlFUWkNMRWRCUVVjc1EwRkJReXhEUVVGRE8wRkJSV3hETEdkRFFVRjVRaXhIUVVGSExFdEJRVXNzUTBGQlF6dEJRVU5zUXl4blEwRkJlVUlzUjBGQlJ5eFJRVUZSTEVOQlFVTTdRVUZEY2tNc2EwTkJRVEpDTEVkQlFVY3NRMEZCUXl4RFFVRkRPMEZCUTJoRExEUkNRVUZ4UWl4SFFVRkhMRVZCUVVVc1EwRkJRenRCUVVVelFpeHZRMEZCTmtJc1IwRkJSeXhOUVVGTkxFTkJRVU03UVVGRGRrTXNkME5CUVdsRExFZEJRVWNzUTBGQlF5eERRVUZETzBGQlEzUkRMR3REUVVFeVFpeEhRVUZITEVOQlFVTXNRMEZCUXlKOSIsIi8vIFRoZSBtb2R1bGUgY2FjaGVcbnZhciBfX3dlYnBhY2tfbW9kdWxlX2NhY2hlX18gPSB7fTtcblxuLy8gVGhlIHJlcXVpcmUgZnVuY3Rpb25cbmZ1bmN0aW9uIF9fd2VicGFja19yZXF1aXJlX18obW9kdWxlSWQpIHtcblx0Ly8gQ2hlY2sgaWYgbW9kdWxlIGlzIGluIGNhY2hlXG5cdHZhciBjYWNoZWRNb2R1bGUgPSBfX3dlYnBhY2tfbW9kdWxlX2NhY2hlX19bbW9kdWxlSWRdO1xuXHRpZiAoY2FjaGVkTW9kdWxlICE9PSB1bmRlZmluZWQpIHtcblx0XHRyZXR1cm4gY2FjaGVkTW9kdWxlLmV4cG9ydHM7XG5cdH1cblx0Ly8gQ3JlYXRlIGEgbmV3IG1vZHVsZSAoYW5kIHB1dCBpdCBpbnRvIHRoZSBjYWNoZSlcblx0dmFyIG1vZHVsZSA9IF9fd2VicGFja19tb2R1bGVfY2FjaGVfX1ttb2R1bGVJZF0gPSB7XG5cdFx0Ly8gbm8gbW9kdWxlLmlkIG5lZWRlZFxuXHRcdC8vIG5vIG1vZHVsZS5sb2FkZWQgbmVlZGVkXG5cdFx0ZXhwb3J0czoge31cblx0fTtcblxuXHQvLyBFeGVjdXRlIHRoZSBtb2R1bGUgZnVuY3Rpb25cblx0X193ZWJwYWNrX21vZHVsZXNfX1ttb2R1bGVJZF0obW9kdWxlLCBtb2R1bGUuZXhwb3J0cywgX193ZWJwYWNrX3JlcXVpcmVfXyk7XG5cblx0Ly8gUmV0dXJuIHRoZSBleHBvcnRzIG9mIHRoZSBtb2R1bGVcblx0cmV0dXJuIG1vZHVsZS5leHBvcnRzO1xufVxuXG4iLCIiLCIvLyBzdGFydHVwXG4vLyBMb2FkIGVudHJ5IG1vZHVsZSBhbmQgcmV0dXJuIGV4cG9ydHNcbi8vIFRoaXMgZW50cnkgbW9kdWxlIGlzIHJlZmVyZW5jZWQgYnkgb3RoZXIgbW9kdWxlcyBzbyBpdCBjYW4ndCBiZSBpbmxpbmVkXG52YXIgX193ZWJwYWNrX2V4cG9ydHNfXyA9IF9fd2VicGFja19yZXF1aXJlX18oXCIuL2luZGV4LnRzXCIpO1xuIiwiIl0sIm5hbWVzIjpbXSwic291cmNlUm9vdCI6IiJ9