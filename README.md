# **TridiagonalMatrix**

A high‑performance Swift library for working with **real and complex tridiagonal matrices**, powered by **Accelerate/LAPACK** and designed with a clean, idiomatic Swift API.

This package provides:

- A strongly‑typed `TridiagonalMatrix<T>` type
- Fast LU factorization via LAPACK (`*gttrf`)
- Solving linear systems (`*gttrs`)
- Condition number estimation (`*gtcon`)
- Determinant computation
- Real and complex matrix–vector multiplication
- Optimized AXPY and VDSP‑accelerated operations
- Full support for `Float`, `Double`, `Complex<Float>`, and `Complex<Double>`

The API is fully generic over a `ScalarField` protocol, enabling future extensions to custom numeric types.

---

## **Features**

### ✔ Fully generic over real and complex scalar fields

The `ScalarField` and `RealScalar` protocols abstract over the LAPACK/BLAS backends, allowing:

- `Float`
- `Double`
- `Complex<Float>`
- `Complex<Double>`

to all share the same high‑level API.

### ✔ Fast LU factorization

`TridiagonalMatrix.factorized()` returns a `TridiagonalLUMatrix` that stores:

- L (unit lower bidiagonal)
- U (upper bidiagonal + second superdiagonal)
- Pivot indices
- Determinant
- Reciprocal condition number

### ✔ Solve linear systems

```swift
var lu = A.factorized()
var b = [1, 2, 3]
let x = lu.solve(&b)
```

Supports:

- Single RHS (`[T]`)
- Multiple RHS in column‑major (`[T]`)
- Multiple RHS as `[[T]]`

### ✔ Optimized matrix–vector multiply

```swift
let y = A * x
```

Real types use vDSP; complex types use fused real/imag operations.

### ✔ Optimized AXPY

```swift
let y = aXpY(a: 2.0, x: x, y: y0)
```

### ✔ Determinant and condition number

```swift
let det = lu.determinant
let cond = lu.approximateConditionNumber
```

---

## **Installation**

Add this to your `Package.swift`:

```swift
.package(url: "https://github.com/josephlevy222/TridiagonalMatrix.git", from: "1.0.0")
```

Then add the dependency to your target:

```swift
.target(
    name: "MyApp",
    dependencies: [
        .product(name: "TridiagonalMatrix", package: "TridiagonalMatrix")
    ]
)
```

---

## **Usage**

### **Constructing a tridiagonal matrix**

```swift
let A = TridiagonalMatrix(
    diagonal: [4, 5, 6],
    upper:    [1, 1],
    lower:    [2, 2]
)
```

### **Factorization**

```swift
var lu = A.factorized()
```

### **Solving a system**

```swift
var b = [7, 8, 9]
let x = lu.solve(&b)
```

### **Matrix–vector multiply**

```swift
let y = A * x
```

### **Condition number**

```swift
print(lu.approximateConditionNumber)
```

---

## **Design Overview**

### **Protocols**

#### `ScalarField`

Defines the operations required for a scalar type to participate in tridiagonal algebra:

- LAPACK `gttrf`, `gttrs`, `gtcon`
- BLAS `axpy`
- High‑level `AXpY` and `multiply`

#### `RealScalar`

Extends `ScalarField` with:

- Complex LAPACK backends
- vDSP operations

### **Matrix Types**

#### `TridiagonalMatrix<T>`

Immutable representation of a tridiagonal matrix.

#### `TridiagonalLUMatrix<T>`

Mutable LU factorization with:

- Pivoting
- Determinant
- Condition number
- Solve routines

### **Operations**

- `AXpY` (real & complex)
- `multiply` (real & complex)
- `*` operator overload

---

## **Performance**

This library is built directly on:

- **LAPACK** (`*gttrf`, `*gttrs`, `*gtcon`)
- **BLAS** (`*axpy`)
- **vDSP** (real multiply‑add, multiply, subtract)

All operations are:

- In‑place where possible
- Cache‑friendly
- Zero‑copy when bridging to LAPACK
- Optimized for both real and complex types

---

## **Directory Structure**

```
Sources/
  TridiagonalMatrix/
    LAPACKTypes.swift
    ScalarField.swift
    RealScalar.swift
    LAPACKWrappers.swift
    Float+ScalarField.swift
    Double+ScalarField.swift
    Complex+ScalarField.swift
    TridiagonalMatrix.swift
    TridiagonalLUMatrix.swift
    TridiagonalWorkspace.swift
    AXpY.swift
    Multiply.swift
```

---

## **Testing**

You can test correctness using small known systems:

- Diagonally dominant matrices
- Random matrices with known solutions
- Identity and scalar multiples
- Symmetric positive definite tridiagonal matrices

I can generate a full test suite if you want.

---

## **License**

MIT 

---

## **Contributing**

Pull requests are welcome.  
If you want to extend the library (e.g., eigenvalues, banded matrices, sparse formats), open an issue.

---

