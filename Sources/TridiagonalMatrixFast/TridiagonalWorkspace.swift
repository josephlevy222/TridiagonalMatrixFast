//
//  TridiagonalWorkspace.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

/// Workspace for LAPACK tridiagonal routines.
///
/// Manages:
/// - Real/complex work arrays
/// - Integer work arrays (real types only)
public final class TridiagonalWorkspace<T: ScalarField> {
	
	public private(set) var work: CMutablePtr<T.CType>
	public private(set) var iwork: CMutablePtr<CInt>?
	
	private var capacityWork: Int
	private var capacityIWork: Int
	
	/// Creates a workspace with an initial capacity.
	///
	/// - Parameter capacity: Initial matrix dimension estimate.
	public init(capacity: Int = 0) {
		let n = max(1, capacity)
		self.capacityWork = 2 * n
		self.work = CMutablePtr<T.CType>.allocate(capacity: capacityWork)
		
		// Real types (Float/Double) require iwork for gtcon; complex backends do not
		if T.self is any RealScalar.Type {
			self.capacityIWork = n
			self.iwork = CMutablePtr<CInt>.allocate(capacity: capacityIWork)
		} else {
			self.capacityIWork = 0
			self.iwork = nil
		}
	}
	
	deinit {
		work.deallocate()
		iwork?.deallocate()
	}
	
	/// Returns a work buffer of sufficient size for a matrix of dimension `n`.
	///
	/// - Parameter n: Matrix dimension.
	public func workBuffer(for n: Int) -> CMutablePtr<T.CType> {
		let required = max(1, 2 * n)
		if capacityWork < required {
			work.deallocate()
			capacityWork = max(required, capacityWork * 2)
			work = CMutablePtr<T.CType>.allocate(capacity: capacityWork)
		}
		return work
	}
	
	/// Returns an integer work buffer for real types.
	///
	/// - Parameter n: Matrix dimension.
	/// - Returns: A pointer to an integer work buffer, or `nil` for complex types.
	public func iworkBuffer(for n: Int) -> CMutablePtr<CInt>? {
		guard T.self is any RealScalar.Type else { return nil }
		let required = max(1, n)
		if capacityIWork < required {
			iwork?.deallocate()
			capacityIWork = max(required, capacityIWork * 2)
			iwork = CMutablePtr<CInt>.allocate(capacity: capacityIWork)
		}
		return iwork
	}
}
