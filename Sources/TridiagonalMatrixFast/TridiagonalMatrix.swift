//
//  TridiagonalMatrix.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// A tridiagonal matrix with elements of type `T`.
///
/// Storage layout:
/// - `lower[i]`   stores A[i+1, i]
/// - `diagonal[i]` stores A[i, i]
/// - `upper[i]`   stores A[i, i+1]
///
/// The matrix has size `n × n` where `n = diagonal.count`.
@frozen
public struct TridiagonalMatrix<T: ScalarField> {
	
	/// Subdiagonal elements (count = n - 1).
	public var lower: [T]
	
	/// Main diagonal elements (count = n).
	public var diagonal: [T]
	
	/// Superdiagonal elements (count = n - 1).
	public var upper: [T]
	
	/// Matrix dimension.
	public let size: Int
	
	/// Creates a tridiagonal matrix from its three diagonals.
	///
	/// - Parameters:
	///   - diagonal: Main diagonal (count = n)
	///   - upper: Superdiagonal (count = n - 1)
	///   - lower: Subdiagonal (count = n - 1)
	public init(diagonal: [T], upper: [T], lower: [T]) {
		precondition(diagonal.count > 0, "Diagonal must not be empty")
		precondition(diagonal.count == upper.count + 1, "Invalid upper size")
		precondition(diagonal.count == lower.count + 1, "Invalid lower size")
		self.diagonal = diagonal
		self.upper = upper
		self.lower = lower
		self.size = diagonal.count
	}
	
	/// Computes the matrix 1‑norm (maximum absolute column sum).
	@inlinable
	public func oneNorm() -> T.Magnitude {
		let n = size
		guard n > 0 else { return 0 }
		if n == 1 { return diagonal[0].magnitude }
		
		var norm = diagonal[0].magnitude + lower[0].magnitude
		if n > 2 {
			for j in 1..<(n - 1) {
				let col = upper[j - 1].magnitude
				+ diagonal[j].magnitude
				+ lower[j].magnitude
				norm = max(norm, col)
			}
		}
		norm = max(norm, upper[n - 2].magnitude + diagonal[n - 1].magnitude)
		return norm
	}
	
	/// Returns the LU factorization of this matrix.
	@inlinable
	public func factorized() -> TridiagonalLUMatrix<T> {
		TridiagonalLUMatrix(self)
	}
}

/// A column vector type alias for tridiagonal operations.
public typealias ColumnVector<T: ScalarField> = [T]
