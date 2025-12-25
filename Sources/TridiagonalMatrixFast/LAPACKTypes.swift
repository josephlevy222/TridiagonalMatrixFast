//
//  LAPACKTypes.swift
//  TridiagonalMatrix
//
//  Created by Joseph Levy on 12/24/25.
//

import Accelerate
import Numerics

// MARK: - Core C Typealiases

/// C integer type used by LAPACK.
public typealias CInt = __CLPK_integer

/// Mutable pointer to C memory.
public typealias CMutablePtr<T> = UnsafeMutablePointer<T>

/// Raw pointer to C memory.
public typealias CVoidPtr = UnsafeRawPointer

/// Mutable raw pointer to C memory.
public typealias CMutableVoidPtr = UnsafeMutableRawPointer

// MARK: - LAPACK / BLAS Function Typealiases

/// LAPACK tridiagonal LU factorization (`*gttrf`).
public typealias GttrfFunction<T> = (
	_ N: CMutablePtr<CInt>?, _ DL: CMutablePtr<T>?, _ D: CMutablePtr<T>?,
	_ DU: CMutablePtr<T>?, _ DU2: CMutablePtr<T>?, _ IPIV: CMutablePtr<CInt>?,
	_ INFO: CMutablePtr<CInt>?
) -> CInt

/// LAPACK tridiagonal triangular solve (`*gttrs`).
public typealias GttrsFunction<T> = (
	_ TRANS: CMutablePtr<Int8>?, _ N: CMutablePtr<CInt>?, _ NRHS: CMutablePtr<CInt>?,
	_ DL: CMutablePtr<T>?, _ D: CMutablePtr<T>?, _ DU: CMutablePtr<T>?,
	_ DU2: CMutablePtr<T>?, _ IPIV: CMutablePtr<CInt>?, _ B: CMutablePtr<T>?,
	_ LDB: CMutablePtr<CInt>?, _ INFO: CMutablePtr<CInt>?
) -> CInt

/// LAPACK reciprocal condition number estimator (`*gtcon`) for real types.
public typealias GtconFunction<T, M, W> = (
	_ NORM: CMutablePtr<Int8>?, _ N: CMutablePtr<CInt>?, _ DL: CMutablePtr<T>?,
	_ D: CMutablePtr<T>?, _ DU: CMutablePtr<T>?, _ DU2: CMutablePtr<T>?,
	_ IPIV: CMutablePtr<CInt>?, _ anorm: CMutablePtr<M>?, _ rcond: CMutablePtr<M>?,
	_ WORK: CMutablePtr<W>?, _ IWORK: CMutablePtr<CInt>?, _ info: CMutablePtr<CInt>?
) -> CInt

/// LAPACK reciprocal condition number estimator (`*gtcon`) for complex types.
public typealias GtconComplexFunction<T, M, W> = (
	_ NORM: CMutablePtr<Int8>?, _ N: CMutablePtr<CInt>?, _ DL: CMutablePtr<T>?,
	_ D: CMutablePtr<T>?, _ DU: CMutablePtr<T>?, _ DU2: CMutablePtr<T>?,
	_ IPIV: CMutablePtr<CInt>?, _ anorm: CMutablePtr<M>?, _ rcond: CMutablePtr<M>?,
	_ WORK: CMutablePtr<W>?, _ info: CMutablePtr<CInt>?
) -> CInt

/// BLAS AXPY signature: `y ← a·x + y`.
public typealias BLASAxpyFunction<T> = (
	_ n: Int32, _ a: T, _ x: UnsafePointer<T>, _ incx: Int32,
	_ y: CMutablePtr<T>, _ incy: Int32
) -> Void

/// vDSP vector operation signature: `c = a ∘ b` with strides.
public typealias DSPFunction<T> = (
	_ a: UnsafePointer<T>, _ sa: Int, _ b: UnsafePointer<T>, _ sb: Int,
	_ c: CMutablePtr<T>, _ sc: Int, _ n: Int
) -> Void

/// High-level multiply-add: `y ← A·x + y`.
public typealias MultiplyAdd<T: ScalarField> = (
	_ A: TridiagonalMatrix<T>, _ x: [T], _ y: inout [T]
) -> [T]

/// High-level matrix–vector multiply: `A·x`.
public typealias MatrixVectorMultiply<T: ScalarField> = (
	_ A: TridiagonalMatrix<T>, _ x: [T]
) -> [T]

