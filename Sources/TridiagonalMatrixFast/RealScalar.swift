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
	static var cgttrf: GttrfFunction<WType> { get }
	
	/// Complex tridiagonal triangular solve.
	static var cgttrs: GttrsFunction<WType> { get }
	
	/// Complex reciprocal condition number estimator.
	static var cgtcon: GtconFunction<WType, Self, WType> { get }
	
	/// Complex AXPY.
	static var caxpy: BLASAxpyFunction<WType> { get }
	
	/// vDSP multiply-add for real scalars.
	static var vma: DSPFunction<Self> { get }
	
	/// vDSP multiply for real scalars.
	static var vmul: DSPFunction<Self> { get }
	
	/// vDSP subtract for real scalars.
	static var vsub: DSPFunction<Self> { get }
	
	/// Complex multiply-add: `y ← A·x + y`.
	static var cAXpY: MultiplyAdd<Complex<Self>> { get }
	
	/// Complex matrix–vector multiply.
	static var cMultiply: MatrixVectorMultiply<Complex<Self>> { get }
}
