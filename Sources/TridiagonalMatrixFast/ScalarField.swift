//
//  ScalarField.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// A numeric scalar type that supports tridiagonal matrix operations.
///
/// Conforming types must provide:
/// - LAPACK tridiagonal LU factorization (`gttrf`)
/// - LAPACK triangular solve (`gttrs`)
/// - LAPACK condition estimation (`gtcon`)
/// - BLAS AXPY (`axpy`)
/// - High-level multiply-add (`AXpY`)
/// - High-level matrix–vector multiply (`multiply`)
public protocol ScalarField: AlgebraicField where Magnitude: FloatingPoint {
	
	/// C-compatible storage type for this scalar.
	associatedtype CType = Self
	
	/// Multiplicative identity.
	static var one: Self { get }
	
	/// LAPACK tridiagonal LU factorization.
	static var gttrf: gttrf<Self> { get }
	
	/// LAPACK tridiagonal triangular solve.
	static var gttrs: gttrs<Self> { get }
	
	/// LAPACK reciprocal condition number estimator.
	static var gtcon: gtcon<Self, Magnitude, CType> { get }
	
	/// BLAS AXPY: `y ← a·x + y`.
	static var axpy: axpy<Self> { get }
	
	/// High-level multiply-add: `y ← A·x + y`.
	static var AXpY: MultiplyAdd<Self> { get }
	
	/// High-level matrix–vector multiply: `A·x`.
	static var multiply: MatrixVectorMultiply<Self> { get }
}

extension ScalarField {
	public static var one: Self { 1 }
}
