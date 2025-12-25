//
//  LAPACKWrappers.swift
//  TridiagonalMatrix
//

import Accelerate
import Numerics

// MARK: - gtcon Wrappers

/// Wrapper for real `gtcon` that allocates `iwork` if needed.
///
/// - Parameters mirror the LAPACK `*gtcon` signature, with optional `iwork`.
@inlinable
public func gtcon_<T, R>(
	_ cFunction: @escaping gtcon<T, R, T>,
	_ norm: CMutablePtr<Int8>?, _ n: CMutablePtr<CInt>?,
	_ dl: CMutablePtr<T>?, _ d: CMutablePtr<T>?, _ du: CMutablePtr<T>?,
	_ du2: CMutablePtr<T>?, _ ipiv: CMutablePtr<CInt>?,
	_ anorm: CMutablePtr<R>?, _ rcond: CMutablePtr<R>?,
	_ work: CMutablePtr<T>?, _ iwork: CMutablePtr<CInt>?,
	_ info: CMutablePtr<CInt>?
) -> CInt {
	if let iwork {
		return cFunction(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
	} else {
		let nn = max(1, Int(n!.pointee))
		return withUnsafeTemporaryAllocation(of: CInt.self, capacity: nn) { tempI in
			cFunction(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, tempI.baseAddress!, info)
		}
	}
}

/// Wrapper for complex `gtcon` that allocates `work` if needed.
///
/// - Note: `iwork` is ignored for complex variants.
@inlinable
public func gtcon_<T, R>(
	_ cFunction: @escaping gtconComplex<T, R, T>,
	_ norm: CMutablePtr<Int8>?, _ n: CMutablePtr<CInt>?,
	_ dl: CMutablePtr<T>?, _ d: CMutablePtr<T>?, _ du: CMutablePtr<T>?,
	_ du2: CMutablePtr<T>?, _ ipiv: CMutablePtr<CInt>?,
	_ anorm: CMutablePtr<R>?, _ rcond: CMutablePtr<R>?,
	_ work: CMutablePtr<T>?, _ iwork: CMutablePtr<CInt>?,
	_ info: CMutablePtr<CInt>?
) -> CInt {
	if let work {
		return cFunction(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
	} else {
		let nn = max(1, Int(n!.pointee) * 2)
		return withUnsafeTemporaryAllocation(of: T.self, capacity: nn) { tmp in
			cFunction(norm, n, dl, d, du, du2, ipiv, anorm, rcond, tmp.baseAddress!, info)
		}
	}
}

// MARK: - AXPY Wrappers

/// AXPY wrapper for complex types, bridging through raw pointers.
///
/// - Parameters:
///   - cFunction: CBLAS complex AXPY function
///   - n: Vector length
///   - a: Scalar multiplier
///   - x: Input vector
///   - y: Output vector (updated in-place)
@inlinable
public func axpy_<CType>(
	_ cFunction: @escaping (
		_ N: Int32, _ a: CVoidPtr?, _ X: CVoidPtr?, _ incX: Int32,
		_ Y: CMutableVoidPtr?, _ incY: Int32
	) -> Void,
	_ n: Int32, _ a: CType, _ x: UnsafePointer<CType>, _ incx: Int32,
	_ y: CMutablePtr<CType>, _ incy: Int32
) {
	withUnsafePointer(to: a) { aptr in
		cFunction(n, aptr, x, incx, y, incy)
	}
}

/// Direct AXPY call for real types.
///
/// - Parameters:
///   - cFunction: CBLAS real AXPY function
///   - n: Vector length
///   - a: Scalar multiplier
///   - x: Input vector
///   - y: Output vector (updated in-place)
@inlinable
public func axpy_<T>(
	_ cFunction: @escaping axpy<T>,
	_ n: Int32, _ a: T, _ x: UnsafePointer<T>, _ incx: Int32,
	_ y: CMutablePtr<T>, _ incy: Int32
) {
	cFunction(n, a, x, incx, y, incy)
}


