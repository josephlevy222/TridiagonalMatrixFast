//
//  RealScalar.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// A real scalar type that also provides complex LAPACK backends.
///
/// Adopted by `Float` and `Double`.
public protocol RealScalar: ScalarField & FloatingPoint & Real {
	
	/// C complex storage type (`__CLPK_complex` or `__CLPK_doublecomplex`).
	associatedtype WType
	
	/// Complex tridiagonal LU factorization.
	static var cgttrf: gttrf<WType> { get }
	
	/// Complex tridiagonal triangular solve.
	static var cgttrs: gttrs<WType> { get }
	
	/// Complex reciprocal condition number estimator.
	static var cgtcon: gtcon<WType, Self, WType> { get }
	
	/// Complex AXPY.
	static var caxpy: axpy<WType> { get }
	
	/// vDSP multiply-add for real scalars.
	static var vma: DSPSignature<Self> { get }
	
	/// vDSP multiply for real scalars.
	static var vmul: DSPSignature<Self> { get }
	
	/// vDSP subtract for real scalars.
	static var vsub: DSPSignature<Self> { get }
	
	/// Complex multiply-add: `y ← A·x + y`.
	static var cAXpY: MultiplyAdd<Complex<Self>> { get }
	
	/// Complex matrix–vector multiply.
	static var cMultiply: MatrixVectorMultiply<Complex<Self>> { get }
}
